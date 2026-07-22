.formula_backtick <- function(x) {
  paste0("`", gsub("`", "``", x, fixed = TRUE), "`")
}

.mgcv_explicit_formula <- function(response, rhs) {
  if (!length(rhs)) return(stats::as.formula(paste(response, "~ 1")))
  stats::as.formula(paste(response, "~", paste(.formula_backtick(rhs), collapse = " + ")))
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

.mgcv_fit_engine <- function(n, mgcv_model = NULL) {
  if (is.null(mgcv_model)) {
    mgcv_model <- if (n < 50000L) "gam" else "bam"
  } else {
    if (!is.character(mgcv_model) || length(mgcv_model) != 1L ||
        is.na(mgcv_model)) {
      stop("mgcv_model must be NULL, 'gam', or 'bam'.")
    }
    mgcv_model <- tolower(mgcv_model)
    if (!mgcv_model %in% c("gam", "bam")) {
      stop("mgcv_model must be NULL, 'gam', or 'bam'.")
    }
  }
  list(
    fit = if (identical(mgcv_model, "gam")) mgcv::gam else mgcv::bam,
    method = if (identical(mgcv_model, "gam")) "REML" else "fREML",
    model = mgcv_model
  )
}

.mgcv_patch_family_environment <- function(family) {
  fam_name <- tolower(paste(family$family, collapse = " "))
  if (grepl("tweedie", fam_name, fixed = TRUE)) {
    ld <- get("ldTweedie", envir = asNamespace("mgcv"))
    for (nm in names(family)) {
      if (is.function(family[[nm]])) {
        env <- environment(family[[nm]])
        if (!environmentIsLocked(env) &&
            !exists("ldTweedie", envir = env, inherits = TRUE)) {
          assign("ldTweedie", ld, envir = env)
        }
      }
    }
  }
  family
}

.mgcv_fit_explicit <- function(response, rhs, data, family,
                               mgcv_model = NULL) {
  engine <- .mgcv_fit_engine(nrow(data), mgcv_model)
  family <- .mgcv_patch_family_environment(family)
  engine$fit(
    .mgcv_explicit_formula(response, rhs), data = data,
    family = family, method = engine$method
  )
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

.mgcv_fit_init <- function(X, response_info, Z, selected, family,
                           mgcv_model = NULL) {
  Xinit <- NULL
  if (length(selected) > 0L) {
    Xinit <- X[, selected, drop = FALSE]
    colnames(Xinit) <- paste0("InitX", seq_along(selected))
  }
  pred <- .mgcv_predictor_data(Z, Xinit, n = response_info$n)
  dat <- cbind(response_info$data, pred)
  rhs <- colnames(pred)
  .mgcv_fit_explicit(
    response_info$response, rhs, dat, family,
    mgcv_model = mgcv_model
  )
}

.mgcv_greedy_warm_start <- function(X, response_info, Z, family, L.init = 1,
                                    mgcv_model = NULL) {
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- .mgcv_fit_init(
    X = X, response_info = response_info, Z = Z,
    selected = selected, family = family, mgcv_model = mgcv_model
  )

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "response")
    j <- select_by_residual_score(X = X, residual = r, available = available)
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    fit <- .mgcv_fit_init(
      X = X, response_info = response_info, Z = Z,
      selected = selected, family = family, mgcv_model = mgcv_model
    )
  }

  fit
}

.mgcv_extract_working <- function(fit, weight_cutoff = 0.0025) {
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

#' General mgcv IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @param family A GLM or mgcv family object that provides `variance()`
#'   and `mu.eta()` for working IRLS.
#' @param mgcv_model Either `NULL`, `"gam"`, or `"bam"`.
#'   `NULL` uses `gam` when `n < 50000` and `bam` otherwise.
#' @importFrom mgcv gam bam nb tw betar scat
#' @keywords internal
#' @noRd
Run_GLM <- function(X, y, Z = NULL, weight_cutoff = 0.0025,
                    family = binomial(link = "logit"),
                    mgcv_model = NULL,
                    L, max.iter, min.iter, max.eps, susie_para,
                    verbose = TRUE, n_threads = 1,
                    L.init = 1,
                    noncs_var = 0.1,
                    noncs_max_abs_cor = 0.9,
                    suff_block_size = 10000L) {

  run_start <- proc.time()[["elapsed"]]
  n <- NROW(y)
  p <- ncol(X)
  suff_block_size <- validate_suff_block_size(suff_block_size)
  .mgcv_validate_family(family)
  is_gaussian <- identical(family$family, "gaussian") &&
    identical(family$link, "identity")
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
    L.init = L.init, mgcv_model = mgcv_model
  )

  alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
  g <- numeric(0)
  beta <- rep(0, p)
  beta_prev <- beta
  alpha_prev <- alpha * 0
  fitX <- NULL
  XCS <- NULL

  fitX_no_cs_streak <- 0L
  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha

    work <- .mgcv_extract_working(fit_final, weight_cutoff = weight_cutoff)
    suff <- weighted_residual_suffstats(
      X = X,
      y = work$pseudo_response,
      ZI = ZI,
      weights = if (is_gaussian) work$W_diag else work$W_diag / work$phi0,
      n_threads = n_threads,
      block_size = suff_block_size
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

    beta <- clean_coef(stats::coef(fitX)[-1])
    CSdt <- summary(fitX)$vars
    cs_indices <- sort(unique(CSdt$cs[CSdt$cs > 0]))
    fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L

    rm(suff)

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
        noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Z
      )
      if (!is.null(noncs_term)) {
        XCS_refit <- cbind(XCS_refit, Main_noncs_res = noncs_term)
      }
    }

    pred <- .mgcv_predictor_data(Z, XCS_refit)
    Data <- cbind(response_info$data, pred)
    penalty_names <- grep("^(Main_CS[0-9]+|Main_noncs_res)$", colnames(XCS_refit), value = TRUE)
    penalty_V <- .refit_penalty_variance(fitX, cs_indices, penalty_names)
    fit_final <- .mgcv_fit_fixed_ridge(
      response_info$response, colnames(pred), Data, family, penalty_V,
      dispersion = work$phi0, mgcv_model = mgcv_model
    )

    alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
    err <- max(sqrt(mean((beta - beta_prev)^2)),
               sqrt(mean((alpha - alpha_prev)^2)))
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
    pred <- .mgcv_predictor_data(Z, XCS_refit)
    Data <- cbind(response_info$data, pred)
    penalty_names <- grep("^(Main_CS[0-9]+|Main_noncs_res)$", colnames(XCS_refit), value = TRUE)
    penalty_V <- .refit_penalty_variance(fitX, cs_indices, penalty_names)
    fit_final <- .mgcv_fit_fixed_ridge(
      response_info$response, colnames(pred), Data, family, penalty_V,
      dispersion = refit_dispersion, mgcv_model = mgcv_model
    )
  }

  G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
  if (!is.null(G)) MainIndex <- safe_add_p(MainIndex, G)
  fit_final$n_eff <- if (exists("work")) work$n_eff else NA_real_
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

  list(
    diagnostics = make_diagnostics(
      if (exists("iter")) iter else 0L, g, run_start
    ),
    fitX = fitX,
    fitJoint = fit_final,
    discovery_summary = MainIndex
  )
}
