#' Cox score IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @param status Event indicator for Cox proportional-hazards outcomes.
#' @export
Run_Cox <- function(X, y, status, Z = NULL,
                    L, max.iter, min.iter, max.eps, susie.iter,
                    verbose = TRUE, n_threads = 1, coverage = 0.9,
                    estimate_residual_variance = TRUE, scaled_prior_variance = 1,
                    residual_variance = 0.5,
                    residual_variance_lowerbound = 0.1,
                    residual_variance_upperbound = 1,
                    ridge = 1e-6,
                    L.init = 1,
                    init_cor_method = NULL,
                    refit_noncs = TRUE,
                    noncs_var = 0.2,
                    suff_block_size = 10000L, ...) {

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

  # Survival outcome
  surv_y = survival::Surv(y, status)

  # ============================================
  # Greedy low-dimensional Cox warm start
  # ============================================
  fit_final = greedy_cox_warm_start(
    X = X, y = y, status = status, Z = Z, L.init = L.init,
    init_cor_method = init_cor_method
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
  early_no_cs <- FALSE
  XCS <- NULL

  # ============================================
  # Main iteration loop
  # ============================================
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

    ss  = cox_suffstat(X = XZE, eta = eta, time = y,
                       status = as.integer(status), n_threads = n_threads)
    a     = as.numeric(ss$a)
    B     = as.matrix(ss$B)
    XZEty = as.numeric(ss$Xty)
    n_eff = ss$d

    XZEa = XZE * sqrt(a)
    A    = blockwise_crossprod(XZEa, n_threads = n_threads,
                               block_size = suff_block_size)
    BtB  = blockwise_crossprod(B, n_threads = n_threads,
                               block_size = suff_block_size)
    XZEtXZE = A - BtB
    XZEtXZE = (XZEtXZE + t(XZEtXZE)) / 2
    idxX = seq_len(p)
    idxE = p + 1L
    idxZ = p + 1L + seq_len(q)

    # Information blocks.
    XtX = XZEtXZE[idxX, idxX, drop = FALSE]
    XtE = XZEtXZE[idxX, idxE, drop = FALSE]   # X'W eta
    XtZ = XZEtXZE[idxX, idxZ, drop = FALSE]
    EtZ = XZEtXZE[idxE, idxZ, drop = FALSE]
    ZtZ = XZEtXZE[idxZ, idxZ, drop = FALSE]
    ZtX = XZEtXZE[idxZ, idxX, drop = FALSE]
    ZtE = XZEtXZE[idxZ, idxE, drop = FALSE]

    XtM = XZEty[idxX]

    # Project the (X, eta) block against Z.
    diag(ZtZ) = diag(ZtZ) + ridge
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
    fitX <- susieR::susie_ss(
      XtX = XtX, Xty = Xty, yty = n - 1, n = n, L = L,
      residual_variance = residual_variance,
      scaled_prior_variance = scaled_prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = residual_variance_upperbound,
      max_iter = susie.iter,
      estimate_prior_method = "EM",
      coverage = coverage, ...
    )

    beta = clean_coef(coef(fitX)[-1])

    # Extract credible sets using summary information
    CSdt <- summary(fitX)$vars
    cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
    cs_indices = sort(cs_indices)

    if (length(cs_indices) == 0) {
      if (iter <= min.iter) {
        noncs_res <- build_no_cs_noncs_refit_term(X, fitX)
        if (is.null(noncs_res)) {
          early_no_cs <- TRUE
          if (verbose) {
            cat("No credible set detected; returning current no-CS fit.\n")
          }
          break
        }
        XCS <- matrix(noncs_res, ncol = 1)
        colnames(XCS) <- "Main_CS_noncs"
        XCS <- as.matrix(XCS)
        XCS_refit <- XCS
      } else {
        early_no_cs <- TRUE
        if (verbose) {
          cat("No credible set detected; returning current no-CS fit.\n")
        }
        break
      }
    } else {

    Alpha_filtered <- fitX$alpha * 0
    for (i in cs_indices) {
      vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
      Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i]
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
    if (isTRUE(refit_noncs)) {
      noncs_term <- build_noncs_refit_term(
        X = X, fitX = fitX, CSdt = CSdt, cs_indices = cs_indices,
        XCS = XCS, noncs_var = noncs_var
      )
      if (!is.null(noncs_term)) {
        XCS_refit <- cbind(XCS_refit, Main_CS_noncs = noncs_term)
      }
    }
    }

    # ============================================
    # Refit Cox with selected credible sets
    # ============================================
    if (ncol(Z) == 0) {
      Data = data.frame(XCS_refit)
    } else {
      Data = cbind(Z, XCS_refit)
      Data = as.data.frame(Data)
    }

    fit_final = survival::coxph(surv_y ~ ., data = Data, ties = "breslow")

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

    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  # ============================================
  # Post-processing
  # ============================================
  if (early_no_cs) {
    MainIndex <- Identifying_MainEffect(fitX, colnames(X))
    G <- summary(fit_final)$coefficients
    MainIndex <- safe_add_p(MainIndex, G)
    fit_final <- clean_model_environment(fit_final)
    return(list(
      iter = iter,
      error = g,
      converged = FALSE,
      fitX = fitX,
      fitJoint = fit_final,
      main_index = MainIndex,
      JointCoef = G
    ))
  }

  if (ncol(Z) == 0) {
    Data = data.frame(XCS)
  } else {
    Data = cbind(Z, XCS)
    Data = as.data.frame(Data)
  }
  fit_final = survival::coxph(surv_y ~ ., data = Data, ties = "breslow")
  MainIndex = Identifying_MainEffect(fitX, colnames(X))
  G = summary(fit_final)$coefficients[, -2, drop = FALSE]
  MainIndex <- safe_add_p(MainIndex, G)
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
    iter = iter,
    error = g,
    converged = (iter < max.iter && err < max.eps),
    fitX = fitX,
    fitJoint = fit_final,
    main_index = MainIndex,
    JointCoef = G
  )

  return(AA)
}
