#' SuSiE_IRLS: SuSiE for non-Gaussian outcomes via IRLS + SuSiE
#'
#' Iterates between IRLS-based pseudo-response construction and SuSiE fine-mapping
#' on covariate-adjusted data. At each outer iteration the algorithm:
#' 1) fits a GLM to obtain pseudo responses and weights, 2) projects out covariates,
#' 3) runs SuSiE on X, 4) refits a GLM with selected variables, and checks convergence.
#'
#' For survival outcomes, pass a `survival::Surv` object as `y`. The Cox path
#' builds score-based sufficient statistics and runs SuSiE-SS without IRLS
#' weights.
#' Credible-set refit terms use their current absolute SuSiE prior variances as
#' fixed Gaussian ridge penalties. A one-column aggregate non-CS term receives
#' the largest positive finite current variance; an error is raised only when
#' no positive finite SuSiE variance is available.
#' Covariates in `Z` are not penalized. mgcv refits use the previous model dispersion to convert
#' the absolute precision `1 / V` to the penalized-deviance scale; Cox refits
#' use `survival::ridge(theta = 1 / V, scale = FALSE)`.
#'
#' @param X An n by p numeric matrix of predictors.
#' @param y Response vector, or a `survival::Surv` object for Cox PH.
#' @param Z An n by q matrix or vector of covariates. If NULL, only an intercept is used.
#' @param family A supported GLM/mgcv family or dispatch string. Ordered-logit
#'   models use fixed-ridge mgcv refits; explicitly non-logit cumulative links
#'   use `ordinal::clm()`. Ignored when `y` is a `Surv` object.
#' @param mgcv_model Either `NULL`, `"gam"`, or `"bam"` for ordinary GLM and
#'   ZIP refits. `NULL` uses `gam` when `n < 50000` and `bam` otherwise.
#' @param L Number of single effects in SuSiE. Default 10.
#' @param L.init Number of predictors used in the low-dimensional fitting step.
#'   Default 1.
#' @param max.iter Maximum outer iterations. Default 10.
#' @param min.iter Minimum outer iterations before convergence can be declared.
#'   Default 2.
#' @param max.eps Convergence threshold on max parameter change. Default 1e-5.
#' @param verbose Logical flag for progress printing. Default TRUE.
#' @param n_threads Integer number of threads for internal parallel blocks. Default 4.
#' @param susie_para Named `susieR::susie_ss()` options. Structural sufficient
#'   statistics and `L` are managed by SuSiEIRLS. The package field
#'   `prior_variance` is an absolute coefficient prior variance, fixed when
#'   `estimate_prior_variance = FALSE` and used as an estimate initializer when
#'   `TRUE`. The legacy `scaled_prior_variance` name is accepted with a warning
#'   and interpreted identically; the two names are mutually exclusive. Through
#'   `min.iter`, the prior-variance value and estimation switch are replaced by
#'   fixed `scaled_prior_variance = 2` and `estimate_prior_variance = FALSE`,
#'   while all other `susie_para` settings remain active.
#' @param weight_cutoff Quantile in (0, 0.05) to clip extreme IRLS weights. Default 0.0025.
#' @param noncs_var Minimum non-CS variance fraction required to add the
#'   non-CS residual summary variable. For example, `noncs_var = 0.1`
#'   adds it when the CS summary explains less than 90% of the posterior mean
#'   linear predictor variance. Default 0.1.
#' @param noncs_max_abs_cor Maximum allowed absolute correlation between a
#'   non-CS term and applicable existing refit covariates, including `Z` and
#'   credible-set summary variables. Correlations equal to or above the
#'   threshold are rejected. The effective threshold is capped at 0.9;
#'   smaller user-supplied values remain effective. Default 0.9.
#' @param scale_data Logical. If TRUE, standardize `X` with
#'   `SuSiE4I::large_scale()` and center and scale non-binary columns of `Z`;
#'   binary columns of `Z` remain on their original scale. If FALSE, `X` and
#'   non-binary columns of `Z` are assumed to already
#'   be standardized. Default TRUE.
#' @param suff_block_size Row block size for weighted sufficient-statistic
#'   crossproducts. Larger values can be faster for small-to-moderate p when
#'   memory is sufficient. Default 10000.
#' @return A list containing the main-effect SuSiE fit, final joint model, and
#'   `discovery_summary` table. `diagnostics` is a one-row data frame
#'   containing the number of outer iterations, final convergence eps, and
#'   runtime in seconds. The effective sample size used by the algorithm is
#'   stored in `fitJoint$n_eff`.
#'
#' @importFrom stats var coef glm binomial quantile sd
#' @importFrom graphics text
#' @importFrom susieR susie_ss coef.susie
#' @importFrom survival coxph Surv
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixCor
#' @importFrom mgcv gam bam nb tw betar scat ziP ocat
#' @importFrom ordinal clm
#' @export
SuSiE_IRLS <- function(X, Z = NULL, y,
                       family = binomial(link = "logit"),
                       mgcv_model = NULL,
                       n_threads = 4, L = 10,
                       susie_para = NULL,
                       max.iter = 10, max.eps = 1e-5, min.iter = 2,
                       weight_cutoff = 0.0025,
                       L.init = 1,
                       noncs_var = 0.1,
                       noncs_max_abs_cor = 0.9,
                       scale_data = TRUE,
                       suff_block_size = 10000L,
                       verbose = TRUE) {

  # ---- basic checks ----
  if (is.null(X)) stop("X must not be NULL.")
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (ncol(X) == 0) stop("X has zero columns.")
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  n <- nrow(X)

  if (is.null(y)) stop("y must be provided.")
  if (NROW(y) != n) stop("Length(y) must equal nrow(X).")

  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z)) stop("Z must be numeric.")
    if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
    if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  }

  if (!is.logical(scale_data) || length(scale_data) != 1L || is.na(scale_data)) {
    stop("scale_data must be TRUE or FALSE.")
  }

  if (!is.numeric(weight_cutoff) || length(weight_cutoff) != 1L || !is.finite(weight_cutoff)) {
    stop("weight_cutoff must be a finite numeric scalar.")
  }
  if (weight_cutoff <= 0) weight_cutoff <- 1e-6
  if (weight_cutoff >= 0.05) weight_cutoff <- 0.049
  noncs_max_abs_cor <- validate_noncs_max_abs_cor(noncs_max_abs_cor)

  # ---- optional standardization ----
  if (isTRUE(scale_data)) {
    x_dimnames <- dimnames(X)
    X <- SuSiE4I:::large_scale(X, n_threads = n_threads, center = TRUE, scale = TRUE)
    X <- as.matrix(X)
    dimnames(X) <- x_dimnames

    if (!is.null(Z)) {
      z_dimnames <- dimnames(Z)
      Z <- apply(Z, 2L, function(x) {
        vals <- unique(stats::na.omit(x))
        if (length(vals) <= 2L) return(as.numeric(x))
        as.numeric(scale(x))
      })
      Z <- as.matrix(Z)
      dimnames(Z) <- z_dimnames
    }
  } else if (isTRUE(verbose)) {
    cat("scale_data = FALSE: assuming X is column-centered and standardized; non-binary columns of Z should also be standardized.\n")
  }

  # ---- helpers ----
  family_string <- if (is.character(family) && length(family) == 1L) {
    tolower(family)
  } else {
    NULL
  }
  if (identical(family_string, "gaussian")) {
    family <- stats::gaussian()
    family_string <- NULL
  }
  is_zip_flag <- .zip_is_family(family)
  clm_links <- c("logit", "probit", "cloglog", "loglog", "cauchit")
  clm_link <- NULL
  if (!is.null(family_string) && startsWith(family_string, "clm_")) {
    clm_link <- sub("^clm_", "", family_string)
    if (!clm_link %in% clm_links) {
      stop("Unsupported CLM link. Use clm_logit or clm_probit, for example.")
    }
  }
  if (identical(family_string, "clm")) {
    stop("family = 'clm' is incomplete. Use clm_logit or clm_probit, for example.")
  }
  is_clm_flag <- !is.null(clm_link) && !identical(clm_link, "logit")
  ordinal_strings <- c("ordinal", "ocat", "clm_logit")
  ordered_default <- is.ordered(y) && nlevels(y) >= 3L &&
    inherits(family, "family") && identical(family$family, "binomial")
  is_ocat_flag <- .ocat_is_family(family) ||
    (!is.null(family_string) && family_string %in% ordinal_strings) ||
    ordered_default
  # Cox is identified by a Surv-typed response; family is then ignored.
  is_cox_flag <- inherits(y, "Surv")
  susie_para <- .resolve_susie_para(susie_para)
  suff_block_size <- validate_suff_block_size(suff_block_size)

  # ---- dispatch ----
  if (is_cox_flag) {
    # Surv carries time and status; split and hand off to Run_Cox.
    surv_time   <- as.numeric(y[, 1])
    surv_status <- as.integer(y[, 2])
    return(
      Run_Cox(
        X = X, y = surv_time, status = surv_status, Z = Z,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie_para = susie_para,
        verbose = verbose,
        n_threads = n_threads,
        L.init = L.init,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size
      )
    )
  }

  if (is_zip_flag) {
    return(
      Run_ZIP(
        X = X, y = y, Z = Z,
        family = family,
        mgcv_model = mgcv_model,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie_para = susie_para,
        verbose = verbose,
        n_threads = n_threads,
        weight_cutoff = weight_cutoff,
        L.init = L.init,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size
      )
    )
  }

  if (is_clm_flag) {
    return(
      Run_CLM(
        X = X, y = y, Z = Z, clm_link = clm_link,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie_para = susie_para,
        verbose = verbose,
        n_threads = n_threads,
        L.init = L.init,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size
      )
    )
  }

  if (is_ocat_flag) {
    ocat_family <- if (.ocat_is_family(family)) {
      family
    } else {
      mgcv::ocat(R = .ocat_prepare_response(y)$ncat)
    }
    return(
      Run_OCAT(
        X = X, y = y, Z = Z, family = ocat_family,
        weight_cutoff = weight_cutoff,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie_para = susie_para,
        verbose = verbose,
        n_threads = n_threads,
        L.init = L.init,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size
      )
    )
  }

  if (is.character(family)) {
    stop(
      "Unsupported family string. Use 'ordinal', 'ocat', 'clm_logit', ",
      "an explicit non-logit clm_<link>, or supply a family object."
    )
  }

  # general GLM (e.g., Poisson, Gaussian with non-default link, etc.)
  return(
    Run_GLM(
      X = X, y = y, Z = Z, family = family,
      mgcv_model = mgcv_model,
      L = L, max.iter = max.iter, min.iter = min.iter, max.eps = max.eps,
      susie_para = susie_para, verbose = verbose, n_threads = n_threads,
      weight_cutoff = weight_cutoff,
      L.init = L.init,
      noncs_var = noncs_var,
      noncs_max_abs_cor = noncs_max_abs_cor,
      suff_block_size = suff_block_size
    )
  )
}
