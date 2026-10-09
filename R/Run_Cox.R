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

# B' diag(dev) B and B' diag(dev) BN for the risk-set means B of a geno X at the
# event times, streamed over row blocks in descending-time order so the n by p
# matrix B is never formed (as in SuSiE4I's cox_suffstat_block).
cox_riskset_crossprod_geno <- function(X, eta, time, status, dev, BN,
                                       block_size = 10000L) {
  n <- nrow(X)
  p <- ncol(X)
  status <- as.integer(status)
  ord <- order(time, decreasing = TRUE)
  w <- exp(eta - max(eta))
  run <- rle(time[ord])
  run_id <- rep(seq_along(run$lengths), run$lengths)
  has_event <- tabulate(run_id[status[ord] == 1L], length(run$lengths)) > 0
  blk_end <- cumsum(run$lengths)[has_event]
  S0 <- cumsum(w[ord])[blk_end]
  rows <- max(1L, min(as.integer(block_size), 2^26 %/% p))
  BtdB <- matrix(0, p, p)
  BtdN <- matrix(0, p, ncol(BN))
  carry <- numeric(p)
  for (start in seq.int(1L, n, by = rows)) {
    idx <- start:min(n, start + rows - 1L)
    G <- X[ord[idx], , drop = FALSE]
    Gw <- G * w[ord[idx]]
    Gw[1L, ] <- Gw[1L, ] + carry
    S1 <- matrix(apply(Gw, 2L, cumsum), nrow = length(idx))
    carry <- S1[length(idx), ]
    b_idx <- which(blk_end >= start & blk_end <= max(idx))
    if (length(b_idx) == 0L) next
    Bblk <- S1[blk_end[b_idx] - start + 1L, , drop = FALSE] / S0[b_idx]
    BtdB <- BtdB + crossprod(Bblk, Bblk * dev[b_idx])
    BtdN <- BtdN + crossprod(Bblk, BN[b_idx, , drop = FALSE] * dev[b_idx])
  }
  list(XtWX = BtdB, XtM = BtdN)
}

Run_Cox <- function(X, y, status, Z = NULL,
                    L, max.iter, min.iter, max.eps, susie_para,
                    verbose = TRUE, n_threads = 1,
                    ridge = 1e-6,
                    L.init = 1,
                    noncs_var = 0.1,
                    noncs_max_abs_cor = 0.9,
                    lbf_threshold = 1,
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
  XCS_refit <- NULL

  # ============================================
  # Main iteration loop
  # ============================================
  fitX_no_cs_streak <- 0L
  for (iter in seq_len(max.iter)) {
    beta_prev = beta
    alpha_prev = alpha

    ## ===== Cox score-based sufficient statistics with binary-style projection =====

    # Current linear predictor
    eta = fit_final$linear.predictors

    # Same projection logic as Run_Binary:
    # ZI always contains intercept; ZI is used only for projection.
    q = ncol(ZI)

    N = cbind(eta, ZI)
    k = ncol(N)
    rsN = SuSiE4I::cox_riskset(X = N, eta = eta, time = y,
                               status = as.integer(status), n_threads = 1L)
    # a, M, dev and d do not depend on the columns, so a geno X reuses rsN.
    rsX = if (inherits(X, "geno")) rsN else {
      SuSiE4I::cox_riskset(X = X, eta = eta, time = y,
                           status = as.integer(status), n_threads = n_threads)
    }
    a     = as.numeric(rsX$a)
    M     = as.numeric(rsX$M)
    dev   = as.numeric(rsX$dev)
    n_eff = rsX$d

    # [X eta ZI]' diag(a) [X eta ZI] - B' diag(dev) B, split into X and N blocks.
    AX = SuSiE4I::weighted_crossprod(X, a, cbind(N * a, M),
                                     n_threads = n_threads, block_size = suff_block_size)
    BX = if (inherits(X, "geno")) {
      cox_riskset_crossprod_geno(X, eta = eta, time = y, status = status,
                                 dev = dev, BN = rsN$B,
                                 block_size = suff_block_size)
    } else {
      SuSiE4I::weighted_crossprod(rsX$B, dev, rsN$B * dev,
                                  n_threads = n_threads, block_size = suff_block_size)
    }
    XN = AX$XtM[, seq_len(k), drop = FALSE] - BX$XtM
    NN = crossprod(N, N * a) - crossprod(rsN$B, rsN$B * dev)
    NN = (NN + t(NN)) / 2

    # Information blocks.
    XtX = AX$XtWX - BX$XtWX
    XtX = (XtX + t(XtX)) / 2
    dimnames(XtX) = list(colnames(X), colnames(X))
    XtE = XN[, 1L, drop = FALSE]
    XtZ = XN[, 1L + seq_len(q), drop = FALSE]
    EtZ = NN[1L, 1L + seq_len(q), drop = FALSE]
    ZtZ = NN[1L + seq_len(q), 1L + seq_len(q), drop = FALSE]
    ZtX = t(XtZ)
    ZtE = NN[1L + seq_len(q), 1L, drop = FALSE]

    XtM = as.numeric(AX$XtM[, k + 1L])

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
      iter, min.iter, lbf_threshold = lbf_threshold
    )
    fitX <- do.call(susieR::susie_ss, ss_args)

    kill <- !is.null(lbf_threshold) && iter > min.iter
    design <- build_refit_design(
      X, fitX, cor_design = Z, kill = kill,
      lbf_threshold = lbf_threshold,
      noncs_var = noncs_var, noncs_max_abs_cor = noncs_max_abs_cor,
      verbose = verbose
    )
    beta <- design$beta
    cs_indices <- design$cs_indices
    XCS_refit <- design$XCS_refit
    fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L

    # ============================================
    # Refit Cox with selected credible sets
    # ============================================
    penalty_V <- design$penalty_V
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
  penalty_V <- design$penalty_V
  fit_final <- .cox_fit_fixed_ridge(
    y, status, Z, Xextra = XCS_refit, penalty_V = penalty_V
  )
  MainIndex = Identifying_MainEffect(fitX, colnames(X), keep_cs = design$cs_indices)
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
