.ocat_extract_working <- function(fit, weight_cutoff = 0.0025) {
  eta <- as.numeric(fit$linear.predictors)
  eta <- pmin(pmax(eta, -50), 50)
  residual <- as.numeric(fit$residuals)
  W_diag <- fit$working.weights
  if (is.null(W_diag)) W_diag <- fit$weights
  W_diag <- as.numeric(W_diag)

  bad <- !is.finite(eta) | !is.finite(residual) |
    !is.finite(W_diag) | W_diag <= 0
  if (mean(bad) > 0.9) stop("Too many invalid ocat working observations.")
  eta[bad] <- 0
  residual[bad] <- 0
  W_diag[bad] <- 0

  pseudo_response <- eta + residual
  W_diag <- robust_weight(W_diag, cutoff = weight_cutoff)
  weight_denom <- sum(W_diag^2)
  if (!is.finite(weight_denom) || weight_denom <= 0) {
    stop("All ocat working weights are zero.")
  }

  list(
    pseudo_response = pseudo_response,
    W_diag = W_diag,
    n_eff = (sum(W_diag)^2) / weight_denom
  )
}

#' Ordered-logit SuSiE runner using mgcv::ocat
#'
#' @inheritParams SuSiE_IRLS
#' @param family An `mgcv::ocat(R = )` family object.
#' @keywords internal
#' @noRd
Run_OCAT <- function(X, y, Z = NULL, weight_cutoff = 0.0025,
                     family,
                     L, max.iter, min.iter, max.eps, susie_para,
                     verbose = TRUE, n_threads = 1,
                     L.init = 1,
                     noncs_var = 0.1,
                     noncs_max_abs_cor = 0.9,
                     suff_block_size = 10000L) {

  run_start <- proc.time()[["elapsed"]]
  n <- NROW(y)
  p <- ncol(X)
  if (!.ocat_is_family(family)) {
    stop("Run_OCAT requires an mgcv::ocat() family.")
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))
  suff_block_size <- validate_suff_block_size(suff_block_size)
  y_info <- .ocat_prepare_response(y, family = family)
  response_info <- list(
    data = data.frame(y = y_info$y_int), response = "y", n = n
  )

  if (is.null(Z)) {
    Z <- matrix(nrow = n, ncol = 0)
    ZI <- matrix(1, nrow = n, ncol = 1)
    colnames(ZI) <- "Intercept"
  } else {
    if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
    if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    ZI <- cbind(Intercept = 1, Z)
  }
  q <- ncol(Z)

  fit_final <- .mgcv_greedy_warm_start(
    X = X, response_info = response_info, Z = Z, family = family,
    L.init = L.init, mgcv_model = "gam"
  )

  alpha <- .ocat_nuisance_coef(fit_final, q = q)
  g <- numeric(0)
  beta <- rep(0, p)
  beta_prev <- beta
  alpha_prev <- alpha * 0
  fitX <- NULL
  XCS <- NULL
  work <- NULL

  fitX_no_cs_streak <- 0L
  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha

    work <- .ocat_extract_working(
      fit_final, weight_cutoff = weight_cutoff
    )
    suff <- weighted_residual_suffstats(
      X = X, y = work$pseudo_response, ZI = ZI,
      weights = work$W_diag,
      n_threads = n_threads, block_size = suff_block_size
    )

    ss_args <- .susie_iteration_args(
      susie_para,
      list(
        XtX = suff$XtX, Xty = suff$Xty, yty = suff$yty,
        n = max(0.95 * n, work$n_eff), L = L
      ),
      iter, min.iter
    )
    fitX <- do.call(susieR::susie_ss, ss_args)
    rm(suff)

    beta <- clean_coef(stats::coef(fitX)[-1L])
    CSdt <- summary(fitX)$vars
    cs_indices <- sort(unique(CSdt$cs[CSdt$cs > 0]))
    fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L

    if (!length(cs_indices)) {
      noncs_res <- build_no_cs_noncs_refit_term(
        X, fitX, cor_design = Z,
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
        vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
        Alpha_filtered[i, vars_in_cs_i] <-
          fitX$alpha[i, vars_in_cs_i] /
          sum(fitX$alpha[i, vars_in_cs_i])
      }
      Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
      XCS <- CppMatrix::matrixMultiply(
        X, as.matrix(Alpha_filtered), transB = TRUE
      )
      XCS <- XCS[, cs_indices, drop = FALSE]
      if (is.null(dim(XCS))) XCS <- matrix(XCS, ncol = 1)
      colnames(XCS) <- paste0("Main_CS", cs_indices)
      XCS_refit <- XCS

      noncs_term <- build_noncs_refit_term(
        X = X, fitX = fitX, CSdt = CSdt, cs_indices = cs_indices,
        XCS = XCS, noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Z
      )
      if (!is.null(noncs_term)) {
        XCS_refit <- cbind(XCS_refit, Main_noncs_res = noncs_term)
      }
    }

    pred <- .mgcv_predictor_data(Z, XCS_refit, n = n)
    Data <- cbind(response_info$data, pred)
    penalty_names <- grep(
      "^(Main_CS[0-9]+|Main_noncs_res)$",
      colnames(XCS_refit), value = TRUE
    )
    penalty_V <- .refit_penalty_variance(
      fitX, cs_indices, penalty_names
    )
    fit_final <- .mgcv_fit_fixed_ridge(
      response_info$response, colnames(pred), Data, family, penalty_V,
      dispersion = 1, mgcv_model = "gam"
    )

    alpha <- .ocat_nuisance_coef(fit_final, q = q)
    err <- max(
      sqrt(mean((beta - beta_prev)^2)),
      sqrt(mean((alpha - alpha_prev)^2))
    )
    g[iter] <- err

    if (verbose) {
      cat(sprintf(
        "Iteration %d: err = %.3e, n_eff = %.1f\n",
        iter, err, work$n_eff
      ))
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

  MainIndex <- if (is.null(fitX)) NULL else
    Identifying_MainEffect(fitX, colnames(X))
  if (!is.null(XCS_refit)) {
    pred <- .mgcv_predictor_data(Z, XCS_refit, n = n)
    Data <- cbind(response_info$data, pred)
    penalty_names <- grep(
      "^(Main_CS[0-9]+|Main_noncs_res)$",
      colnames(XCS_refit), value = TRUE
    )
    penalty_V <- .refit_penalty_variance(
      fitX, cs_indices, penalty_names
    )
    fit_final <- .mgcv_fit_fixed_ridge(
      response_info$response, colnames(pred), Data, family, penalty_V,
      dispersion = 1, mgcv_model = "gam"
    )
  }

  G <- summary(fit_final)$p.table
  MainIndex <- safe_add_p(MainIndex, G)
  fit_final$n_eff <- if (!is.null(work)) work$n_eff else NA_real_
  fit_final <- clean_model_environment(fit_final)

  if (verbose && length(g)) {
    plot(
      g, type = "o", col = "black", pch = 16,
      xlab = "Iteration", ylab = "Max Parameter Change",
      main = "Convergence Trace (mgcv ocat)"
    )
    for (i in seq_along(g)) {
      graphics::text(
        x = i, y = g[i],
        labels = formatC(g[i], format = "e", digits = 1),
        pos = 3, cex = 0.7, col = "red"
      )
    }
  }

  list(
    diagnostics = make_diagnostics(
      if (exists("iter")) iter else 0L, g, run_start
    ),
    fitX = fitX,
    fitJoint = fit_final,
    discovery_summary = MainIndex
  )
}
