.mgcv_backtick <- function(x) {
  paste0("`", gsub("`", "``", x, fixed = TRUE), "`")
}

.mgcv_explicit_formula <- function(response, rhs) {
  if (!length(rhs)) return(stats::as.formula(paste(response, "~ 1")))
  stats::as.formula(paste(response, "~", paste(.mgcv_backtick(rhs), collapse = " + ")))
}

.mgcv_validate_family <- function(family) {
  if (!inherits(family, "family")) {
    stop("family must be a GLM or mgcv family object.")
  }

  fam_name <- tolower(paste(family$family, collapse = " "))
  fam_class <- tolower(paste(class(family), collapse = " "))
  blocked <- c("zero inflated", "zip", "ordered", "categorical",
               "cox", "censored")
  if (any(vapply(blocked, grepl, logical(1), x = paste(fam_name, fam_class),
                 fixed = TRUE))) {
    stop("Unsupported family for the mgcv IRLS path: ", family$family)
  }
  if (!is.function(family$variance) || !is.function(family$mu.eta)) {
    stop("family must provide variance() and mu.eta() for working IRLS.")
  }
  invisible(TRUE)
}

.mgcv_fit_engine <- function(n) {
  if (n < 30000L) mgcv::gam else mgcv::bam
}

.mgcv_fit_explicit <- function(response, rhs, data, family) {
  if (nrow(data) < 30000L) {
    mgcv::gam(.mgcv_explicit_formula(response, rhs), data = data,
              family = family, method = "REML")
  } else {
    mgcv::bam(.mgcv_explicit_formula(response, rhs), data = data,
              family = family, method = "fREML")
  }
}

.mgcv_prepare_response <- function(y, family) {
  is_binom <- identical(family$family, "binomial") ||
    identical(family$family, "quasibinomial")
  if (is.matrix(y)) {
    if (!is_binom || ncol(y) != 2L) {
      stop("Matrix y is only supported for two-column binomial responses.")
    }
    return(list(
      data = data.frame(y_success = y[, 1], y_failure = y[, 2]),
      response = "cbind(y_success, y_failure)",
      n = nrow(y)
    ))
  }

  y <- as.numeric(y)
  if (is_binom) {
    ymax <- max(y, na.rm = TRUE)
    is_count <- all(is.finite(y)) && all(y >= 0) &&
      all(abs(y - round(y)) < sqrt(.Machine$double.eps))
    if (is_count && ymax > 1) {
      trials <- as.integer(ymax)
      if (any(y > trials)) stop("Binomial counts exceed inferred trial size.")
      return(list(
        data = data.frame(y_success = y, y_failure = trials - y),
        response = "cbind(y_success, y_failure)",
        n = length(y)
      ))
    }
  }

  list(data = data.frame(y = y), response = "y", n = length(y))
}

.mgcv_predictor_data <- function(Z = NULL, Xextra = NULL, n = NULL) {
  if (is.null(n)) {
    n <- if (!is.null(Z)) nrow(Z) else nrow(Xextra)
  }
  out <- data.frame(row.names = seq_len(n))
  if (!is.null(Z) && ncol(Z) > 0L) {
    Zdf <- as.data.frame(Z)
    colnames(Zdf) <- paste0("Z", seq_len(ncol(Z)))
    out <- cbind(out, Zdf)
  }
  if (!is.null(Xextra) && ncol(as.matrix(Xextra)) > 0L) {
    Xdf <- as.data.frame(Xextra)
    out <- cbind(out, Xdf)
  }
  out
}

.mgcv_fit_init <- function(X, response_info, Z, selected, family) {
  Xinit <- NULL
  if (length(selected) > 0L) {
    Xinit <- X[, selected, drop = FALSE]
    colnames(Xinit) <- paste0("InitX", seq_along(selected))
  }
  pred <- .mgcv_predictor_data(Z, Xinit, n = response_info$n)
  dat <- cbind(response_info$data, pred)
  rhs <- colnames(pred)
  .mgcv_fit_explicit(response_info$response, rhs, dat, family)
}

.mgcv_greedy_warm_start <- function(X, response_info, Z, family, L.init = 1,
                                    init_cor_method = NULL) {
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- .mgcv_fit_init(
    X = X, response_info = response_info, Z = Z,
    selected = selected, family = family
  )

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "response")
    j <- select_by_residual_score(X = X, residual = r, available = available)
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    fit <- .mgcv_fit_init(
      X = X, response_info = response_info, Z = Z,
      selected = selected, family = family
    )
  }

  fit
}

.mgcv_extract_working <- function(fit, weight_cutoff = 0.005) {
  eta <- as.numeric(fit$linear.predictors)
  mu <- as.numeric(fit$fitted.values)
  y_work <- as.numeric(fit$y)
  fam <- fit$family
  g_prime_mu <- 1 / fam$mu.eta(eta)
  var_mu <- fam$variance(mu)
  w0 <- fit$prior.weights
  if (is.null(w0)) w0 <- rep(1, length(mu))

  pseudo_response <- eta + (y_work - mu) * g_prime_mu
  W_diag <- as.numeric(w0) / (var_mu * g_prime_mu^2)
  bad <- !is.finite(pseudo_response) | !is.finite(W_diag) | W_diag <= 0
  if (mean(bad) > 0.9) stop("Too many invalid working observations.")
  if (any(bad)) {
    W_diag[bad] <- 0
    pseudo_response[bad] <- 0
  }
  W_diag <- robust_weight(W_diag, cutoff = weight_cutoff)
  weight_denom <- sum(W_diag^2)
  if (!is.finite(weight_denom) || weight_denom <= 0) {
    stop("All working weights are zero.")
  }

  phi0 <- tryCatch(summary(fit)$dispersion, error = function(e) NA_real_)
  if (!is.finite(phi0) || phi0 <= 0) phi0 <- 1

  list(
    pseudo_response = pseudo_response,
    W_diag = W_diag,
    phi0 = phi0,
    n_eff = (sum(W_diag)^2) / weight_denom
  )
}

.mgcv_theta <- function(fit) {
  if (!is.function(fit$family$getTheta)) return(NULL)
  theta <- tryCatch(fit$family$getTheta(TRUE), error = function(e) NULL)
  if (is.null(theta) || any(!is.finite(theta))) return(NULL)
  as.numeric(theta)
}

.mgcv_family_with_theta <- function(family, theta = NULL) {
  if (is.null(theta)) return(family)
  fam_name <- tolower(paste(family$family, collapse = " "))
  link <- if (!is.null(family$link)) family$link else NULL

  if (grepl("negative binomial", fam_name, fixed = TRUE)) {
    link_val <- if (is.null(link)) "log" else link
    return(do.call(mgcv::nb, list(theta = theta, link = link_val)))
  }
  if (grepl("tweedie", fam_name, fixed = TRUE)) {
    link_val <- if (is.null(link)) "log" else link
    return(do.call(mgcv::tw, list(theta = theta, link = link_val)))
  }
  if (grepl("beta", fam_name, fixed = TRUE)) {
    link_val <- if (is.null(link)) "logit" else link
    return(do.call(mgcv::betar, list(link = link_val, theta = theta)))
  }
  if (grepl("scaled t", fam_name, fixed = TRUE)) {
    return(mgcv::scat(theta = theta))
  }
  family
}

#' General mgcv IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @importFrom mgcv gam bam nb tw betar scat
#' @export
Run_GLM <- function(X, y, Z = NULL, weight_cutoff = 0.005,
                    family = binomial(link = "logit"),
                    L, max.iter, min.iter, max.eps, susie.iter,
                    verbose = TRUE, n_threads = 1, coverage = 0.9,
                    estimate_residual_variance = TRUE,
                    residual_variance = 0.5, scaled_prior_variance = 1,
                    residual_variance_lowerbound = 0.1,
                    residual_variance_upperbound = 1,
                    L.init = 1,
                    init_cor_method = NULL,
                    refit_noncs = TRUE,
                    noncs_var = 0.2,
                    suff_block_size = 10000L, ...) {

  n <- NROW(y)
  p <- ncol(X)
  suff_block_size <- validate_suff_block_size(suff_block_size)
  .mgcv_validate_family(family)
  response_info <- .mgcv_prepare_response(y, family)
  if (response_info$n != nrow(X)) stop("Length(y) must equal nrow(X).")

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

  fit_final <- .mgcv_greedy_warm_start(
    X = X, response_info = response_info, Z = Z, family = family,
    L.init = L.init, init_cor_method = init_cor_method
  )

  alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
  g <- numeric(0)
  beta <- rep(0, p)
  beta_prev <- beta
  alpha_prev <- alpha * 0
  fitX <- NULL
  XCS <- NULL
  theta_lock <- NULL
  early_no_cs <- FALSE

  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha

    work <- .mgcv_extract_working(fit_final, weight_cutoff = weight_cutoff)
    suff <- weighted_residual_suffstats(
      X = X,
      y = work$pseudo_response,
      ZI = ZI,
      weights = work$W_diag / work$phi0,
      n_threads = n_threads,
      block_size = suff_block_size
    )

    fitX <- susieR::susie_ss(
      XtX = suff$XtX, Xty = suff$Xty, yty = suff$yty,
      n = max(n / 2, work$n_eff), L = L,
      scaled_prior_variance = scaled_prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = residual_variance_upperbound,
      max_iter = susie.iter,
      estimate_prior_method = "EM",
      coverage = coverage, ...
    )
    rm(suff)

    beta <- clean_coef(stats::coef(fitX)[-1])
    CSdt <- summary(fitX)$vars
    cs_indices <- sort(unique(CSdt$cs[CSdt$cs > 0]))

    if (!length(cs_indices)) {
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
      Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
      XCS <- CppMatrix::matrixMultiply(X, as.matrix(Alpha_filtered), transB = TRUE)
      XCS <- XCS[, cs_indices, drop = FALSE]
      if (is.null(dim(XCS))) XCS <- matrix(XCS, ncol = 1)
      colnames(XCS) <- paste0("Main_CS", cs_indices)
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

    pred <- .mgcv_predictor_data(Z, XCS_refit)
    Data <- cbind(response_info$data, pred)
    fit_final <- .mgcv_fit_explicit(
      response_info$response, colnames(pred), Data,
      .mgcv_family_with_theta(family, theta_lock)
    )
    if (iter == min.iter) theta_lock <- .mgcv_theta(fit_final)

    alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
    err <- max(sqrt(mean((beta - beta_prev)^2)),
               sqrt(mean((alpha - alpha_prev)^2)))
    g[iter] <- err

    if (verbose) {
      theta_msg <- if (is.null(theta_lock)) "" else
        sprintf(", theta_lock=%s", paste(signif(theta_lock, 4), collapse = ","))
      cat(sprintf("Iteration %d: err = %.3e, n_eff = %.1f%s\n",
                  iter, err, work$n_eff, theta_msg))
    }

    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  MainIndex <- if (is.null(fitX)) NULL else Identifying_MainEffect(fitX, colnames(X))
  if (!is.null(XCS)) {
    pred <- .mgcv_predictor_data(Z, XCS)
    Data <- cbind(response_info$data, pred)
    fit_final <- .mgcv_fit_explicit(
      response_info$response, colnames(pred), Data,
      .mgcv_family_with_theta(family, theta_lock)
    )
  }

  G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
  if (!is.null(G)) MainIndex <- safe_add_p(MainIndex, G)
  fit_final <- clean_model_environment(fit_final)

  if (verbose && length(g)) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration",
         ylab = "Max Parameter Change",
         main = "Convergence Trace (Max |Delta| in alpha and beta)")
    for (i in seq_along(g)) {
      graphics::text(x = i, y = g[i],
                     labels = formatC(g[i], format = "e", digits = 1),
                     pos = 3, cex = 0.7, col = "red")
    }
  }

  last_err <- if (length(g)) utils::tail(g, 1) else Inf
  list(
    iter = if (exists("iter")) iter else 0,
    error = g,
    converged = early_no_cs || (exists("iter") && iter < max.iter && last_err < max.eps),
    fitX = fitX,
    fitJoint = fit_final,
    main_index = MainIndex,
    JointCoef = G,
    n_eff = if (exists("work")) work$n_eff else NA_real_,
    theta = .mgcv_theta(fit_final)
  )
}
