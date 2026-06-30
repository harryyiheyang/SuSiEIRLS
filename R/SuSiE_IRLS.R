#' SuSiE_IRLS: SuSiE for non-Gaussian outcomes via IRLS + SuSiE
#'
#' Iterates between IRLS-based pseudo-response construction and SuSiE fine-mapping
#' on covariate-adjusted data. At each outer iteration the algorithm:
#' 1) fits a GLM to obtain pseudo responses and weights, 2) projects out covariates,
#' 3) runs SuSiE on X, 4) refits a GLM with selected variables, and checks convergence.
#'
#' For survival outcomes, pass a \code{survival::Surv} object as \code{y}; the
#' algorithm then dispatches to a Cox proportional-hazards fine-mapping routine
#' (\code{Run_Cox}) that builds score-based sufficient statistics and runs
#' SuSiE-SS. No IRLS weights are used in the Cox path.
#'
#' @param X An n by p numeric matrix of predictors.
#' @param y Response vector, or a \code{survival::Surv} object for Cox PH.
#' @param Z An n by q matrix or vector of covariates. If NULL, only an intercept is used.
#' @param family A GLM or mgcv family object (e.g., \code{binomial(link="logit")}),
#'   or the string \code{"negbin"} for \code{mgcv::nb()}. Ignored when
#'   \code{y} is a \code{Surv} object (Cox PH is used).
#' @param logit_method Method for binomial-logit outcomes. \code{"pg"} uses
#'   \code{Run_Binary}; \code{"glm"} uses the general GLM IRLS path.
#' @param L Number of single effects in SuSiE. Default 10.
#' @param L.init Number of SNPs used in the initial low-dimensional warm start.
#'   Default 1.
#' @param max.iter Maximum outer iterations. Default 15.
#' @param min.iter Minimum iterations before convergence check. Default 3.
#' @param max.eps Convergence threshold on max parameter change. Default 1e-5.
#' @param susie.iter Maximum iterations for each SuSiE fit. Default 300.
#' @param verbose Logical flag for progress printing. Default TRUE.
#' @param n_threads Integer number of threads for internal parallel blocks. Default 4.
#' @param coverage Credible set coverage level in SuSiE. Default 0.9.
#' @param estimate_residual_variance Logical for SuSiE residual variance estimation. Default TRUE.
#' @param residual_variance Initial or fixed residual variance. Default 0.5.
#' @param residual_variance_lowerbound Lower bound when estimating residual variance.
#'   Default 0.1; the PG-logit final correction uses 0.5 when this argument is not supplied.
#' @param residual_variance_upperbound Upper bound when estimating residual variance.
#'   If NULL, defaults to 1.
#' @param scaled_prior_variance Prior variance for SuSiE single effects. Default 1.
#' @param weight_cutoff Quantile in (0, 0.05) to clip extreme IRLS weights. Default 0.005.
#' @param theta_init Initial dispersion parameter for negative binomial. Default 10.
#' @param estimate_theta Logical, whether to estimate theta in negative binomial. Default TRUE.
#' @param ridge Diagonal ridge added to the Cox information matrix for positive
#'   definiteness. Used only in the Cox path. Default 1e-6.
#' @param init_cor_method Deprecated and ignored. The greedy warm start now
#'   selects by the largest absolute residual score \code{X' r}.
#' @param refit_noncs Logical. If TRUE, add a one-dimensional non-CS residual
#'   summary variable to the refit model when the current credible-set summary
#'   leaves enough posterior mean variation outside the CS terms. This variable
#'   is used only to improve the next linear predictor estimate and is not
#'   reported as a credible set. Default TRUE.
#' @param noncs_var Minimum non-CS variance fraction required to add the
#'   non-CS residual summary variable. For example, \code{noncs_var = 0.2}
#'   adds it when the CS summary explains less than 80% of the posterior mean
#'   linear predictor variance. Default 0.2.
#' @param suff_block_size Row block size for weighted sufficient-statistic
#'   crossproducts. Larger values can be faster for small-to-moderate p when
#'   memory is sufficient. Default 10000.
#' @param ... Additional arguments passed to the SuSiE fitting routine.
#'
#' @return A list with elements:
#'   \item{iter}{Number of iterations completed}
#'   \item{error}{Max parameter change per iteration}
#'   \item{converged}{Logical}
#'   \item{fitX}{SuSiE results for main effects}
#'   \item{fitJoint}{Final model fit}
#'   \item{main_index}{Summary table of identified effects with PIPs and p-values}
#'   \item{JointCoef}{Coefficient table from the final joint model}
#'
#' @importFrom stats var lm coef glm
#' @importFrom susieR susie_ss coef.susie
#' @importFrom survival coxph Surv
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixCor
#' @importFrom MASS glm.nb negative.binomial
#' @importFrom mgcv gam bam nb tw betar scat
#' @importFrom SuSiE4I blockwise_crossprod large_scale
#' @export
SuSiE_IRLS <- function(X, Z = NULL, y = NULL,
                       family = binomial(link = "logit"),
                       n_threads = 4, L = 10, coverage = 0.9,
                       estimate_residual_variance = TRUE, residual_variance = 0.5,
                       residual_variance_lowerbound = 0.1,
                       residual_variance_upperbound = NULL,
                       scaled_prior_variance = 1,
                       max.iter = 15, max.eps = 1e-5, min.iter = 4,
                       weight_cutoff = 0.005,
                       theta_init = 10, estimate_theta = TRUE,
                       susie.iter = 30,
                       ridge = 1e-6,
                       logit_method = c("pg", "glm"),
                       L.init = 1,
                       init_cor_method = NULL,
                       refit_noncs = TRUE,
                       noncs_var = 0.2,
                       suff_block_size = 10000L,
                       verbose = TRUE, ...) {

  # ---- helpers ----
  is_logit_binomial <- function(fam) {
    inherits(fam, "family") && identical(fam$family, "binomial") && identical(fam$link, "logit")
  }
  is_negbin_flag <- is.character(family) && length(family) == 1 &&
    family %in% c("negbin", "nb", "negative.binomial")
  # Cox is identified by a Surv-typed response; family is then ignored.
  is_cox_flag <- inherits(y, "Surv")
  logit_method <- match.arg(logit_method)
  rv_upper_default <- if (is.null(residual_variance_upperbound)) 1 else residual_variance_upperbound
  suff_block_size <- validate_suff_block_size(suff_block_size)

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
    if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
    if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  }

  if (!is.numeric(weight_cutoff) || length(weight_cutoff) != 1L || !is.finite(weight_cutoff)) {
    stop("weight_cutoff must be a finite numeric scalar.")
  }
  if (weight_cutoff <= 0) weight_cutoff <- 1e-6
  if (weight_cutoff >= 0.05) weight_cutoff <- 0.049

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
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        ridge = ridge,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is_negbin_flag) {
    return(
      Run_GLM(
        X = X, y = y, Z = Z,
        family = if (isTRUE(estimate_theta)) {
          mgcv::nb(theta = NULL)
        } else {
          mgcv::nb(theta = theta_init)
        },
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        weight_cutoff = weight_cutoff,
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is.character(family)) {
    stop("Unsupported family string. Use \"negbin\" or a GLM/mgcv family object.")
  }

  if (is_logit_binomial(family) && identical(logit_method, "pg")) {
    pg_logit_rv_lower <- if (missing(residual_variance_lowerbound)) 0.5 else residual_variance_lowerbound
    pg_logit_rv_upper <- if (is.null(residual_variance_upperbound)) 1 else residual_variance_upperbound
    return(
      Run_Binary(
        X = X, y = y, Z = Z, family = family,
        L = L, max.iter = max.iter, min.iter = min.iter, max.eps = max.eps,
        susie.iter = susie.iter, verbose = verbose, n_threads = n_threads,
        coverage = coverage, weight_cutoff = weight_cutoff,
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = pg_logit_rv_lower,
        residual_variance_upperbound = pg_logit_rv_upper,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  # general GLM (e.g., Poisson, Gaussian with non-default link, etc.)
  return(
    Run_GLM(
      X = X, y = y, Z = Z, family = family,
      L = L, max.iter = max.iter, min.iter = min.iter, max.eps = max.eps,
      susie.iter = susie.iter, verbose = verbose, n_threads = n_threads,
      coverage = coverage, weight_cutoff = weight_cutoff,
      scaled_prior_variance = scaled_prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = rv_upper_default,
      L.init = L.init,
      init_cor_method = init_cor_method,
      refit_noncs = refit_noncs,
      noncs_var = noncs_var,
      suff_block_size = suff_block_size,
      ...
    )
  )
}
