#' Cox score IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @param status Event indicator for Cox proportional-hazards outcomes.
#' @keywords internal
#' @noRd
.cox_fit_fixed_ridge <- function(y, status, Z, Xextra = NULL,
                                 penalty_V = numeric(0)) {
  dat <- data.frame(.time = as.numeric(y), .status = as.integer(status))
  if (!is.null(Z) && ncol(as.matrix(Z)) > 0L) {
    Z <- as.data.frame(Z)
    dat <- cbind(dat, Z)
  }
  if (!is.null(Xextra) && ncol(as.matrix(Xextra)) > 0L) {
    Xextra <- as.data.frame(Xextra)
    dat <- cbind(dat, Xextra)
  }

  unpenalized_terms <- attr(penalty_V, "unpenalized_terms")
  penalty_names <- names(penalty_V)
  penalty_V <- as.numeric(penalty_V)
  names(penalty_V) <- penalty_names
  if (length(penalty_V) &&
      (is.null(penalty_names) || any(!nzchar(penalty_names)) ||
       any(!is.finite(penalty_V) | penalty_V <= 0))) {
    stop("Cox penalty_V must be a named vector of positive finite variances.")
  }
  if (!all(penalty_names %in% names(dat))) {
    stop("Every penalized Cox refit term must occur in the model data.")
  }

  ordinary <- setdiff(names(dat), c(".time", ".status", penalty_names))
  rhs <- if (length(ordinary)) .formula_backtick(ordinary) else character(0)
  if (length(penalty_names)) {
    ridge_rhs <- vapply(seq_along(penalty_names), function(i) {
      paste0(
        "survival::ridge(", .formula_backtick(penalty_names[i]),
        ", theta = ", format(1 / penalty_V[i], digits = 17, scientific = TRUE),
        ", scale = FALSE)"
      )
    }, character(1))
    rhs <- c(rhs, ridge_rhs)
  }
  rhs_text <- if (length(rhs)) paste(rhs, collapse = " + ") else "1"
  form <- stats::as.formula(paste(
    "survival::Surv(.time, .status) ~", rhs_text
  ))
  fit <- survival::coxph(form, data = dat, ties = "breslow")

  if (length(penalty_names)) {
    penalized <- utils::tail(
      seq_along(stats::coef(fit)), length(penalty_names)
    )
    names(fit$coefficients)[penalized] <- penalty_names
    attr(fit, "refit_penalty") <- list(
      V = stats::setNames(penalty_V, penalty_names),
      precision = stats::setNames(1 / penalty_V, penalty_names),
      theta = stats::setNames(1 / penalty_V, penalty_names),
      scale = FALSE,
      unpenalized_terms = unpenalized_terms
    )
  } else if (length(unpenalized_terms)) {
    attr(fit, "refit_penalty") <- list(
      V = numeric(0), precision = numeric(0), theta = numeric(0),
      scale = FALSE, unpenalized_terms = unpenalized_terms
    )
  }
  fit
}

.cox_coef_table <- function(fit) {
  G <- summary(fit)$coefficients
  if (is.null(G) || is.null(dim(G))) return(NULL)
  if ("exp(coef)" %in% colnames(G)) {
    G <- G[, colnames(G) != "exp(coef)", drop = FALSE]
  }
  G
}

Run_Cox <- function(X, y, status, Z = NULL,
                    L, max.iter, min.iter, max.eps, susie_para,
                    verbose = TRUE, n_threads = 1,
                    ridge = 1e-6,
                    L.init = 1,
                    noncs_var = 0.1,
                    noncs_max_abs_cor = 0.9,
                    suff_block_size = 10000L) {

  run_start <- proc.time()[["elapsed"]]
  n = length(y)
  p = ncol(X)
  suff_block_size <- validate_suff_block_size(suff_block_size)

  # ============================================
  # Handle Z edge cases
  # ============================================
  if (is.null(Z)) {
    Z = matrix(nrow = n, ncol = 0)
    ZI = matrix(1, nrow = n, ncol = 1)
    colnames(ZI) = "Intercept"

  } else {
    if (is.null(dim(Z))) {
      Z = matrix(Z, ncol = 1)
    }

    if (is.null(colnames(Z))) {
      colnames(Z) = paste0("Z", seq_len(ncol(Z)))
    }

    ZI = cbind(1, Z)
    colnames(ZI)[1] = "Intercept"
  }

  # ============================================
  # Greedy low-dimensional Cox warm start
  # ============================================
  fit_final = greedy_cox_warm_start(
    X = X, y = y, status = status, Z = Z, L.init = L.init
  )
  if (ncol(Z) == 0) {
    alpha = numeric(0)
  } else {
    alpha = clean_coef(coef(fit_final)[seq_len(ncol(Z))])
  }

  # Initialize tracking variables
  g = c()
  beta = rep(0, p)
  beta_prev = beta
  alpha_prev = alpha * 0
  XCS <- NULL

  # ============================================
  # Main iteration loop
  # ============================================
  fitX_no_cs_streak <- 0L
  for (iter in 1:max.iter) {
    beta_prev = beta
    alpha_prev = alpha

    ## ===== Cox score-based sufficient statistics with binary-style projection =====

    # Current linear predictor
    eta = fit_final$linear.predictors

    # Same projection logic as Run_Binary:
    # ZI always contains intercept; ZI is used only for projection.
    q = ncol(ZI)

    XZE = cbind(X, eta, ZI)

    ss  = SuSiE4I:::cox_suffstat(X = XZE, eta = eta, time = y,
                                 status = as.integer(status), n_threads = n_threads)
    a     = as.numeric(ss$a)
    B     = as.matrix(ss$B)
    XZEty = as.numeric(ss$Xty)
    n_eff = ss$d

    XZEa = XZE * sqrt(a)
    A    = SuSiE4I::blockwise_crossprod(
      XZEa, n_threads = n_threads, block_size = suff_block_size
    )
    BtB  = SuSiE4I::blockwise_crossprod(
      B, n_threads = n_threads, block_size = suff_block_size
    )
    XZEtXZE = A - BtB
    XZEtXZE = (XZEtXZE + t(XZEtXZE)) / 2
    idxX = seq_len(p)
    idxE = p + 1L
    idxZ = p + 1L + seq_len(q)

    # Information blocks.
    XtX = XZEtXZE[idxX, idxX, drop = FALSE]
    XtE = XZEtXZE[idxX, idxE, drop = FALSE]
    XtZ = XZEtXZE[idxX, idxZ, drop = FALSE]
    EtZ = XZEtXZE[idxE, idxZ, drop = FALSE]
    ZtZ = XZEtXZE[idxZ, idxZ, drop = FALSE]
    ZtX = XZEtXZE[idxZ, idxX, drop = FALSE]
    ZtE = XZEtXZE[idxZ, idxE, drop = FALSE]

    XtM = XZEty[idxX]

    # Project the (X, eta) block against Z.
    Zinv_ZtX = solve_with_ridge(ZtZ, ZtX, ridge = ridge)
    Zinv_ZtE = solve_with_ridge(ZtZ, ZtE, ridge = ridge)

    XtX_proj = XtX - matrixMultiply(XtZ, Zinv_ZtX)
    XtE_proj = as.vector(XtE - matrixVectorMultiply(XtZ, Zinv_ZtE))

    # Combine projected information with the Cox score.
    Xty = XtE_proj + XtM
    XtX = XtX_proj

    XtX = (XtX + t(XtX)) / 2
    diag(XtX) = diag(XtX) + ridge

    # Run SuSiE-SS on the Cox score sufficient statistics.
    ss_args <- .susie_iteration_args(
      susie_para,
      list(XtX = XtX, Xty = Xty, yty = n - 1, n = n, L = L),
      iter, min.iter
    )
    fitX <- do.call(susieR::susie_ss, ss_args)

    beta = clean_coef(coef(fitX)[-1])

    # Extract credible sets using summary information
    CSdt <- summary(fitX)$vars
    cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
    cs_indices = sort(cs_indices)
    fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L

    if (length(cs_indices) == 0) {
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
        XCS <- as.matrix(XCS)
        XCS_refit <- XCS
      }
    } else {

    Alpha_filtered <- fitX$alpha * 0
    for (i in cs_indices) {
      vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
      Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i] / sum(fitX$alpha[i, vars_in_cs_i])
    }

    # Align within-CS SNP directions while preserving PIP weights.
    Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
    XCS <- matrixMultiply(X, as.matrix(Alpha_filtered), transB = TRUE)
    XCS <- XCS[, cs_indices, drop = FALSE]

    if (is.null(dim(XCS))) {
      XCS <- matrix(XCS, ncol = 1)
    }

    colnames(XCS) <- paste0("Main_CS", cs_indices)
    XCS <- as.matrix(XCS)
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

    # ============================================
    # Refit Cox with selected credible sets
    # ============================================
    penalty_names <- grep("^(Main_CS[0-9]+|Main_noncs_res)$", colnames(XCS_refit), value = TRUE)
    penalty_V <- .refit_penalty_variance(fitX, cs_indices, penalty_names)
    fit_final <- .cox_fit_fixed_ridge(
      y, status, Z, Xextra = XCS_refit, penalty_V = penalty_V
    )


    # Extract covariate coefficients only
    if (ncol(Z) == 0) {
      alpha = numeric(0)
    } else {
      alpha = clean_coef(coef(fit_final)[seq_len(ncol(Z))])
    }

    # Check convergence
    err = max(
      sqrt(mean((beta - beta_prev)^2)),
      if (length(alpha)) sqrt(mean((alpha - alpha_prev)^2)) else 0
    )
    g[iter] = err

    if (verbose) {
      cat(sprintf("Iteration %d: err = %.3e, events = %d\n", iter, err, n_eff))
      cat("iter", iter, "cs:", cs_indices, "eta sd:", sd(eta),
          "beta nonzero:", which(beta != 0), "\n")
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

  # ============================================
  # Post-processing
  # ============================================
  penalty_names <- grep("^(Main_CS[0-9]+|Main_noncs_res)$", colnames(XCS_refit), value = TRUE)
  penalty_V <- .refit_penalty_variance(fitX, cs_indices, penalty_names)
  fit_final <- .cox_fit_fixed_ridge(
    y, status, Z, Xextra = XCS_refit, penalty_V = penalty_V
  )
  MainIndex = Identifying_MainEffect(fitX, colnames(X))
  G = .cox_coef_table(fit_final)
  MainIndex <- safe_add_p(MainIndex, G)
  fit_final$n_eff <- n_eff
  fit_final <- clean_model_environment(fit_final)

  if (verbose) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration",
         ylab = "Max Parameter Change",
         main = "Convergence Trace (Cox PH, Breslow)")
    for (i in seq_along(g)) {
      text(x = i, y = g[i],
           labels = formatC(g[i], format = "e", digits = 1),
           pos = 3, cex = 0.7, col = "red")
    }
  }

  AA = list(
    diagnostics = make_diagnostics(iter, g, run_start),
    fitX = fitX,
    fitJoint = fit_final,
    discovery_summary = MainIndex
  )
  return(AA)
}
