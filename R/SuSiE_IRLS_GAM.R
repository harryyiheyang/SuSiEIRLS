#' SuSiE_IRLS with an mgcv GAM null model
#'
#' Fine-maps `X` for a GLM/mgcv-family outcome whose covariate adjustment is a
#' GAM null model, written as a `gam()`/`bam()` formula over `data` (for
#' example `y ~ s(age, by = sex) + sex`). The null model is fitted once; at
#' every outer iteration its smooths are integrated out of the SuSiE stage by
#' the penalized projection with \eqn{(B^\top W B + S_\lambda)^{-1}}{(B'WB + S)^-1}
#' on the IRLS weight scale, i.e. the GAM's \eqn{V_p/\phi}{Vp/phi}, where `B`
#' is the null-model design and \eqn{S_\lambda}{S} uses the smoothing
#' parameters of the latest joint refit. The joint refit is the null formula
#' plus the credible-set terms with the fixed ridge \eqn{1/V}{1/V} of
#' [SuSiE_IRLS()]; its smoothing parameters are re-estimated. Otherwise the
#' algorithm is the GLM path of [SuSiE_IRLS()].
#'
#' Every `s()` uses the adaptive Matern basis of mgcv.taps whatever `bs` was
#' written (with a warning); it is built once and reused. A factor `by`
#' becomes a common smooth plus one varying-coefficient smooth per centered
#' contrast; a numeric `by` gives a varying coefficient whose constant column
#' is dropped. `te()`, `ti()` and `t2()` are not supported. All covariates
#' (factor contrasts included) are centered.
#'
#' Zero-inflated, ordinal and Cox outcomes are not supported on this path.
#'
#' @param formula Null-model formula with univariate `s()` terms.
#' @param data Data frame with the response and the null-model covariates.
#' @param X An n by p numeric matrix of predictors.
#' @param family A GLM or mgcv family object (default `binomial()`).
#' @param mgcv_model `NULL` or `"gam"` (REML), or `"bam"` (fREML,
#'   `discrete = TRUE`).
#' @param k Basis dimension of each smooth unless set inside `s()`.
#' @param scale_data Logical. If TRUE, standardize `X` with
#'   `SuSiE4I::large_scale()`; otherwise `X` is only centered. Default TRUE.
#' @inheritParams SuSiE_IRLS
#' @return As the GLM path of [SuSiE_IRLS()], plus the null fit `fitNull`.
#' @export
SuSiE_IRLS_GAM <- function(formula, data, X,
                           family = binomial(link = "logit"),
                           mgcv_model = NULL, k = 10L,
                           n_threads = 4, L = 10,
                           susie_para = NULL,
                           max.iter = 10, max.eps = 1e-5, min.iter = 2,
                           weight_cutoff = 0.0025,
                           noncs_var = 0.1,
                           noncs_max_abs_cor = 0.9,
                           scale_data = TRUE,
                           suff_block_size = 10000L,
                           verbose = TRUE) {
  if (!inherits(formula, "formula")) stop("formula must be a formula.")
  if (!is.data.frame(data)) stop("data must be a data frame.")
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (ncol(X) == 0) stop("X has zero columns.")
  if (nrow(X) != nrow(data)) stop("nrow(X) must equal nrow(data).")
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (!inherits(family, "family")) stop("family must be a GLM or mgcv family object.")
  if (.zip_is_family(family) || .ocat_is_family(family)) {
    stop("SuSiE_IRLS_GAM supports GLM and mgcv families only (no ziP or ocat).")
  }
  .mgcv_validate_family(family)
  .mgcv_fit_engine(nrow(X), mgcv_model)
  if (!is.logical(scale_data) || length(scale_data) != 1L || is.na(scale_data)) {
    stop("scale_data must be TRUE or FALSE.")
  }
  if (!is.numeric(weight_cutoff) || length(weight_cutoff) != 1L || !is.finite(weight_cutoff)) {
    stop("weight_cutoff must be a finite numeric scalar.")
  }
  if (weight_cutoff <= 0) weight_cutoff <- 1e-6
  if (weight_cutoff >= 0.05) weight_cutoff <- 0.049
  noncs_max_abs_cor <- validate_noncs_max_abs_cor(noncs_max_abs_cor)

  x_dimnames <- dimnames(X)
  X <- if (scale_data) {
    as.matrix(SuSiE4I::large_scale(X, n_threads = n_threads, center = TRUE, scale = TRUE))
  } else sweep(X, 2L, colMeans(X))
  dimnames(X) <- x_dimnames

  null <- gam_null_setup(formula, data, family, mgcv_model, k)
  Run_GAM(
    X = X, null = null, family = family, mgcv_model = mgcv_model,
    L = L, max.iter = max.iter, min.iter = min.iter, max.eps = max.eps,
    susie_para = .resolve_susie_para(susie_para), verbose = verbose,
    n_threads = n_threads, weight_cutoff = weight_cutoff,
    noncs_var = noncs_var, noncs_max_abs_cor = noncs_max_abs_cor,
    suff_block_size = validate_suff_block_size(suff_block_size)
  )
}

# Rewrite the null formula (every s() to the frozen AMatern basis), center the
# covariates, fit the null model once, and build the covariate design used by
# the non-CS correlation gate (centered z, and f(z) for smoothed z).
gam_null_setup <- function(formula, data, family, mgcv_model, k) {
  ig <- mgcv::interpret.gam(formula)
  specs <- ig$smooth.spec
  if (any(vapply(specs, inherits, TRUE, what = c("tensor.smooth.spec", "t2.smooth.spec")))) {
    stop("te(), ti() and t2() are not supported; use univariate s() terms.")
  }
  if (!all(vapply(specs, inherits, TRUE, what = "AMatern.smooth.spec"))) {
    warning("All smooths use the adaptive Matern (AMatern) basis; other bs choices were replaced.", call. = FALSE)
  }
  par_vars <- attr(stats::terms(ig$pf), "term.labels")
  smooth_vars <- unique(vapply(specs, function(sp) sp$term, ""))
  by_vars <- setdiff(unique(vapply(specs, function(sp) sp$by, "")), "NA")
  all_vars <- unique(c(par_vars, smooth_vars, by_vars))
  missing_vars <- setdiff(c(all_vars, all.vars(formula[[2L]])), names(data))
  if (length(missing_vars)) {
    stop("Null-model terms must be plain columns of data; not found: ",
         paste(missing_vars, collapse = ", "), ".")
  }

  # Center everything; a factor becomes its centered treatment contrasts.
  new <- data.frame(row.names = seq_len(nrow(data)))
  for (v in all.vars(formula[[2L]])) new[[v]] <- data[[v]]
  cols <- list()
  for (v in all_vars) {
    if (is.numeric(data[[v]])) {
      new[[v]] <- data[[v]] - mean(data[[v]])
      cols[[v]] <- v
    } else {
      f <- as.factor(data[[v]])
      D <- stats::model.matrix(~ f)[, -1L, drop = FALSE]
      colnames(D) <- make.names(paste0(v, levels(f)[-1L]), unique = TRUE)
      for (j in colnames(D)) new[[j]] <- D[, j] - mean(D[, j])
      cols[[v]] <- colnames(D)
    }
  }

  env <- new.env(parent = environment(formula))
  env$.sirls_bases <- list()
  sterm <- function(z, by, kz) {
    if (is.null(env$.sirls_bases[[z]])) env$.sirls_bases[[z]] <- gam_amatern_base(new[[z]], k = kz)
    sprintf("s(%s%s, bs = \"sirlsAM\", xt = list(base = .sirls_bases[[\"%s\"]]))", z,
            if (is.null(by)) "" else paste0(", by = ", by), z)
  }
  rhs <- unlist(cols[par_vars], use.names = FALSE)
  for (sp in specs) {
    kz <- if (sp$bs.dim > 0) sp$bs.dim else k
    if (identical(sp$by, "NA") || !is.numeric(data[[sp$by]])) rhs <- c(rhs, sterm(sp$term, NULL, kz))
    if (!identical(sp$by, "NA")) {
      rhs <- c(rhs, cols[[sp$by]], vapply(cols[[sp$by]], function(b) sterm(sp$term, b, kz), ""))
    }
  }
  if (!length(rhs)) rhs <- "1"
  fml <- stats::as.formula(paste(ig$response, "~", paste(unique(rhs), collapse = " + ")), env = env)

  fit <- .mgcv_fit_explicit(ig$response, character(0), new, family,
                            mgcv_model = mgcv_model, formula = fml)

  Zcor <- as.matrix(new[, unlist(cols, use.names = FALSE), drop = FALSE])
  if (length(smooth_vars)) {
    tt <- stats::predict(fit, type = "terms")
    for (z in smooth_vars) {
      lab <- vapply(fit$smooth, function(sm) if (sm$term == z) sm$label else NA_character_, "")
      fz <- rowSums(tt[, stats::na.omit(lab), drop = FALSE])
      Zcor <- cbind(Zcor, fz - mean(fz))
      colnames(Zcor)[ncol(Zcor)] <- paste0("f(", z, ")")
    }
  }
  if (!ncol(Zcor)) Zcor <- NULL
  list(fit = fit, formula = fml, data = new, response = ig$response,
       B = stats::predict(fit, type = "lpmatrix"), Zcor = Zcor)
}

# S_lambda in the null-model coefficient layout, smoothing parameters taken
# from the latest fit (smooth labels are shared with the joint refit).
gam_null_penalty <- function(fit_null, fit) {
  q <- length(stats::coef(fit_null))
  S <- matrix(0, q, q)
  for (sm in fit_null$smooth) {
    ii <- sm$first.para:sm$last.para
    S[ii, ii] <- S[ii, ii] + fit$sp[[sm$label]] * sm$S[[1L]]
  }
  S
}

# Run_GLM with the linear Z replaced by a GAM null model (see SuSiE_IRLS_GAM):
# the null smooths enter the projection as B with precision S_lambda on the
# scale of the IRLS weights.
Run_GAM <- function(X, null, family, mgcv_model = NULL,
                    L, max.iter, min.iter, max.eps, susie_para,
                    verbose = TRUE, n_threads = 1,
                    weight_cutoff = 0.0025,
                    noncs_var = 0.1,
                    noncs_max_abs_cor = 0.9,
                    suff_block_size = 10000L) {

  run_start <- proc.time()[["elapsed"]]
  n <- nrow(X)
  p <- ncol(X)
  is_gaussian <- identical(family$family, "gaussian") &&
    identical(family$link, "identity")
  B <- null$B
  Zcor <- null$Zcor

  fit_final <- null$fit
  g <- numeric(0)
  beta <- rep(0, p)
  fitX <- NULL
  XCS <- NULL
  XCS_refit <- NULL
  cs_indices <- integer(0)

  fitX_no_cs_streak <- 0L
  for (iter in seq_len(max.iter)) {
    beta_prev <- beta

    work <- .mgcv_extract_working(fit_final, weight_cutoff = weight_cutoff)
    # Gaussian keeps W on the data scale (SuSiE gets phi); otherwise W / phi.
    weight_scale <- if (is_gaussian) 1 else 1 / work$phi0
    S_null <- gam_null_penalty(null$fit, fit_final) * weight_scale
    suff <- weighted_residual_suffstats(
      X = X,
      y = work$pseudo_response,
      ZI = B,
      weights = work$W_diag * weight_scale,
      n_threads = n_threads,
      block_size = suff_block_size,
      nuisance_precision = S_null
    )

    n_ss <- max(0.95 * n, work$n_eff)
    ss_args <- .susie_iteration_args(
      susie_para,
      list(XtX = suff$XtX, Xty = suff$Xty, yty = suff$yty,
           n = n_ss, L = L),
      iter, min.iter
    )
    if (is_gaussian) {
      ss_args$estimate_residual_variance <- FALSE
      ss_args$residual_variance <- work$phi0
    }
    fitX <- do.call(susieR::susie_ss, ss_args)
    rm(suff)

    beta <- susie_main_coef(fitX, p = p)
    CSdt <- summary(fitX)$vars
    cs_list <- susie_cs_list(fitX)
    cs_indices <- cs_list$index
    fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L

    if (!length(cs_indices)) {
      noncs_res <- build_no_cs_noncs_refit_term(
        X, fitX, cor_design = Zcor,
        noncs_max_abs_cor = noncs_max_abs_cor
      )
      if (is.null(noncs_res)) {
        XCS <- NULL
        XCS_refit <- NULL
        if (verbose) {
          cat("No credible set detected; continuing the outer refit without an X term.\n")
        }
      } else {
        XCS <- matrix(noncs_res, ncol = 1)
        colnames(XCS) <- "Main_noncs_res"
        XCS_refit <- XCS
      }
    } else {
      Alpha_filtered <- fitX$alpha * 0
      for (i in cs_indices) {
        vars_in_cs_i <- cs_list$vars[[match(i, cs_list$index)]]
        Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i] / sum(fitX$alpha[i, vars_in_cs_i])
      }
      Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
      XCS <- CppMatrix::matrixMultiply(X, as.matrix(Alpha_filtered), transB = TRUE)
      XCS <- XCS[, cs_indices, drop = FALSE]
      if (is.null(dim(XCS))) XCS <- matrix(XCS, ncol = 1)
      colnames(XCS) <- paste0("Main_CS", cs_indices)
      XCS_refit <- XCS

      noncs_term <- build_noncs_refit_term(
        X = X, fitX = fitX, CSdt = CSdt, cs_indices = cs_indices,
        XCS = XCS, noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Zcor
      )
      if (!is.null(noncs_term)) {
        XCS_refit <- cbind(XCS_refit, Main_noncs_res = noncs_term)
      }
    }

    pred <- .mgcv_predictor_data(NULL, XCS_refit, n = n)
    Data <- cbind(null$data, pred)
    penalty_names <- colnames(pred)
    penalty_V <- .refit_penalty_variance(fitX, cs_indices, penalty_names)
    fit_final <- .mgcv_fit_fixed_ridge(
      null$response, penalty_names, Data, family, penalty_V,
      dispersion = work$phi0, mgcv_model = mgcv_model, formula = null$formula
    )

    err <- sqrt(mean((beta - beta_prev)^2))
    g[iter] <- err

    if (verbose) {
      theta_now <- .mgcv_theta(fit_final)
      theta_msg <- if (is.null(theta_now)) "" else
        sprintf(", theta=%s", paste(signif(theta_now, 4), collapse = ","))
      cat(sprintf("Iteration %d: err = %.3e, n_eff = %.1f%s\n",
                  iter, err, work$n_eff, theta_msg))
    }

    if (fitX_no_cs_streak >= 3L) {
      if (verbose) cat("No main credible set detected in 3 consecutive iterations; stopping.\n")
      break
    }
    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  MainIndex <- if (is.null(fitX)) NULL else Identifying_MainEffect(fitX, colnames(X))
  if (!is.null(XCS_refit)) {
    refit_dispersion <- .mgcv_refit_dispersion(fit_final)
    pred <- .mgcv_predictor_data(NULL, XCS_refit, n = n)
    Data <- cbind(null$data, pred)
    penalty_names <- colnames(pred)
    penalty_V <- .refit_penalty_variance(fitX, cs_indices, penalty_names)
    fit_final <- .mgcv_fit_fixed_ridge(
      null$response, penalty_names, Data, family, penalty_V,
      dispersion = refit_dispersion, mgcv_model = mgcv_model,
      formula = null$formula
    )
  }

  G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
  if (!is.null(G)) MainIndex <- safe_add_p(MainIndex, G)
  fit_final$n_eff <- work$n_eff

  if (verbose && length(g)) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration",
         ylab = "Max Parameter Change",
         main = "Convergence Trace (GAM)")
    for (i in seq_along(g)) {
      graphics::text(x = i, y = g[i],
                     labels = formatC(g[i], format = "e", digits = 1),
                     pos = 3, cex = 0.7, col = "red")
    }
  }

  list(
    diagnostics = make_diagnostics(iter, g, run_start),
    fitNull = null$fit,
    fitX = fitX,
    fitJoint = fit_final,
    discovery_summary = MainIndex
  )
}
