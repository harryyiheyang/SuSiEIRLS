clip_eta <- function(x) {
  pmin(pmax(as.numeric(x), -35), 35)
}

logit_uni_fun <- function(x, y, e, prior_variance,
                          estimate_intercept = 0, ...) {
  v0 <- prior_variance
  off <- attr(y, "hat_etaZ")
  if (is.null(off)) off <- rep(0, length(y))
  e <- as.numeric(e) + as.numeric(off)
  fit <- stats::glm(y ~ x + offset(e), family = stats::binomial())
  co <- summary(fit)$coefficients
  b <- co["x", "Estimate"]
  s <- co["x", "Std. Error"]
  if (!is.finite(b) || !is.finite(s) || s <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  z <- b / s
  lbf <- 0.5 * log(s^2 / (v0 + s^2)) +
    0.5 * z^2 * v0 / (v0 + s^2)
  fit0 <- stats::glm(y ~ 1 + offset(e), family = stats::binomial())
  lrt <- as.numeric(2 * (stats::logLik(fit) - stats::logLik(fit0)))
  lbf <- lbf - 0.5 * z^2 + 0.5 * lrt
  v1 <- 1 / (1 / v0 + 1 / s^2)
  m1 <- v1 * b / s^2
  list(mu = m1, var = v1, lbf = lbf,
       prior_variance = m1^2 + v1, intercept = 0)
}

set_nb_ser_theta <- function(theta) {
  if (!is.finite(theta) || theta <= 0) stop("NB theta must be positive")
  assign(".nb_theta", as.numeric(theta), envir = environment(nb_uni_fun))
  invisible(theta)
}

reset_nb_ser_theta <- function() {
  env <- environment(nb_uni_fun)
  if (exists(".nb_theta", envir = env, inherits = FALSE)) {
    rm(".nb_theta", envir = env)
  }
}

nb_uni_fun <- function(x, y, e, prior_variance,
                       estimate_intercept = 0, ...) {
  v0 <- prior_variance
  theta <- attr(y, "theta")
  if (is.null(theta)) {
    env <- environment(nb_uni_fun)
    if (exists(".nb_theta", envir = env, inherits = FALSE)) {
      theta <- get(".nb_theta", envir = env, inherits = FALSE)
    }
  }
  if (is.null(theta) || !is.finite(theta) || theta <= 0) {
    stop("NB theta attribute is missing or invalid")
  }
  off <- attr(y, "hat_etaZ")
  if (is.null(off)) off <- rep(0, length(y))
  e <- as.numeric(e) + as.numeric(off)
  fam <- MASS::negative.binomial(theta = theta, link = "log")
  fit <- stats::glm(y ~ x + offset(e), family = fam)
  co <- summary(fit)$coefficients
  b <- co["x", "Estimate"]
  s <- co["x", "Std. Error"]
  if (!is.finite(b) || !is.finite(s) || s <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  z <- b / s
  lbf <- 0.5 * log(s^2 / (v0 + s^2)) +
    0.5 * z^2 * v0 / (v0 + s^2)
  fit0 <- stats::glm(y ~ 1 + offset(e), family = fam)
  lrt <- as.numeric(2 * (stats::logLik(fit) - stats::logLik(fit0)))
  lbf <- lbf - 0.5 * z^2 + 0.5 * lrt
  v1 <- 1 / (1 / v0 + 1 / s^2)
  m1 <- v1 * b / s^2
  list(mu = m1, var = v1, lbf = lbf,
       prior_variance = m1^2 + v1, intercept = 0)
}

tweedie_uni_fun <- function(x, y, e, prior_variance,
                            estimate_intercept = 0, ...) {
  v0 <- prior_variance
  tw_p <- attr(y, "tw_p")
  if (is.null(tw_p) || !is.finite(tw_p) || tw_p <= 1 || tw_p >= 2) {
    stop("Tweedie power attribute is missing or invalid")
  }
  tw_phi <- attr(y, "tw_phi")
  if (is.null(tw_phi) || !is.finite(tw_phi) || tw_phi <= 0) {
    stop("Tweedie phi attribute is missing or invalid")
  }
  off <- attr(y, "hat_etaZ")
  if (is.null(off)) off <- rep(0, length(y))
  e <- as.numeric(e) + as.numeric(off)
  fam <- statmod::tweedie(var.power = tw_p, link.power = 0)
  fit <- stats::glm(y ~ x + offset(e), family = fam)
  co <- summary(fit, dispersion = tw_phi)$coefficients
  b <- co["x", "Estimate"]
  s <- co["x", "Std. Error"]
  if (!is.finite(b) || !is.finite(s) || s <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  z <- b / s
  lbf <- 0.5 * log(s^2 / (v0 + s^2)) +
    0.5 * z^2 * v0 / (v0 + s^2)
  lrt <- (fit$null.deviance - fit$deviance) / tw_phi
  lbf <- lbf - 0.5 * z^2 + 0.5 * lrt
  v1 <- 1 / (1 / v0 + 1 / s^2)
  m1 <- v1 * b / s^2
  list(mu = m1, var = v1, lbf = lbf,
       prior_variance = m1^2 + v1, intercept = 0)
}

surv_uni_fun <- function(x, y, e, prior_variance,
                         estimate_intercept = 0, ...) {
  v0 <- prior_variance
  off <- attr(y, "hat_etaZ")
  if (is.null(off)) off <- rep(0, nrow(y))
  e <- as.numeric(e) + as.numeric(off)
  fit <- survival::coxph(y ~ x + offset(e), ties = "breslow")
  co <- summary(fit)$coefficients
  b <- co[1, "coef"]
  s <- co[1, "se(coef)"]
  if (!is.finite(b) || !is.finite(s) || s <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  z <- b / s
  lbf <- 0.5 * log(s^2 / (v0 + s^2)) +
    0.5 * z^2 * v0 / (v0 + s^2)
  lbf <- lbf - 0.5 * z^2 + 0.5 * summary(fit)$logtest["test"]
  v1 <- 1 / (1 / v0 + 1 / s^2)
  m1 <- v1 * b / s^2
  list(mu = m1, var = v1, lbf = lbf,
       prior_variance = m1^2 + v1, intercept = 0)
}

offset_ser <- function(ser_function, offset) {
  offset <- as.numeric(offset)
  force(ser_function)
  function(X, y, o = NULL, ...) {
    if (is.null(o)) o <- rep(0, nrow(X))
    ser_function(X, y, o = as.numeric(o) + offset, ...)
  }
}

fit_logit_hat_eta <- function(y, Z) {
  Zdf <- as.data.frame(Z)
  colnames(Zdf) <- make.names(colnames(Zdf), unique = TRUE)
  fit <- stats::glm(y ~ ., data = data.frame(y = y, Zdf, check.names = FALSE),
                    family = stats::binomial())
  as.numeric(as.matrix(Zdf) %*% stats::coef(fit)[-1L])
}

fit_nb_hat_eta <- function(y, Z) {
  Zdf <- as.data.frame(Z)
  colnames(Zdf) <- make.names(colnames(Zdf), unique = TRUE)
  fit <- MASS::glm.nb(y ~ ., data = data.frame(y = y, Zdf, check.names = FALSE),
                      link = "log")
  list(eta = as.numeric(as.matrix(Zdf) %*% stats::coef(fit)[-1L]),
       theta = as.numeric(fit$theta))
}

fit_cox_hat_eta <- function(y, Z) {
  Zdf <- as.data.frame(Z)
  colnames(Zdf) <- make.names(colnames(Zdf), unique = TRUE)
  dat <- data.frame(
    time = as.numeric(y[, 1]),
    status = as.integer(y[, 2]),
    Zdf,
    check.names = FALSE
  )
  fit <- survival::coxph(survival::Surv(time, status) ~ ., data = dat,
                         ties = "breslow")
  as.numeric(stats::predict(fit, type = "lp"))
}

fit_tweedie_hat_eta <- function(y, Z, tw_p = NULL) {
  if (!is.null(tw_p) && (!is.finite(tw_p) || tw_p <= 1 || tw_p >= 2)) {
    stop("Tweedie power must be in (1, 2)")
  }
  Zdf <- as.data.frame(Z)
  colnames(Zdf) <- make.names(colnames(Zdf), unique = TRUE)
  fam <- mgcv::tw(theta = tw_p, link = "log")
  dat <- data.frame(y = y, Zdf, check.names = FALSE)
  form <- stats::reformulate(colnames(Zdf), response = "y")
  fit <- mgcv::gam(form, data = dat, family = fam, method = "REML")
  p_hat <- as.numeric(fit$family$getTheta(TRUE))
  phi_hat <- as.numeric(summary(fit)$scale)
  list(
    eta = as.numeric(as.matrix(Zdf) %*% stats::coef(fit)[-1L]),
    p = p_hat,
    phi = phi_hat
  )
}

ibss_x_main_index <- function(fit, X, p, coverage = 0.9,
                              min_abs_corr = 0.5) {
  fit0 <- fit
  class(fit0) <- c("susie", "list")
  cs <- susieR::susie_get_cs(
    fit0, X = X, coverage = coverage, min_abs_corr = min_abs_corr
  )
  if (!length(cs$cs)) return(empty_main_index())

  effect_index <- cs$cs_index
  if (is.null(effect_index)) {
    effect_index <- as.integer(gsub("[^0-9]", "", names(cs$cs)))
  }
  if (length(effect_index) != length(cs$cs) || any(!is.finite(effect_index))) {
    effect_index <- seq_along(cs$cs)
  }

  cls <- lapply(cs$cs, function(ii) ii[ii <= p])
  keep <- vapply(cls, length, integer(1)) > 0
  cls <- cls[keep]
  effect_index <- effect_index[keep]
  if (!length(cls)) return(empty_main_index())

  idx <- unlist(cls)
  cs_id <- rep(seq_along(cls), times = vapply(cls, length, integer(1)))
  lbf0 <- if (!is.null(fit$lbf_ser)) fit$lbf_ser else fit$lbf
  lbf <- rep(NA_real_, length(cls))
  if (!is.null(lbf0) && length(lbf0) >= max(effect_index)) {
    lbf <- as.numeric(lbf0[effect_index])
  }

  out <- data.frame(
    Index = idx,
    CS = paste0("Main_CS", cs_id),
    stringsAsFactors = FALSE
  )
  out$PIP <- fit$pip[out$Index]
  out$lbf <- lbf[cs_id]
  out$Pvalue <- NA_real_
  out
}

glmnet_ebic_lambda <- function(fit, n, p, gamma = 1, phi = 1) {
  if (!is.finite(phi) || phi <= 0) stop("EBIC phi must be positive")
  if (!is.finite(gamma) || gamma < 0) stop("EBIC gamma must be nonnegative")
  b <- as.matrix(fit$beta[seq_len(p), , drop = FALSE])
  df <- colSums(is.finite(b) & b != 0)
  dev <- (1 - fit$dev.ratio) * fit$nulldev
  ebic <- dev / phi + (log(n) + gamma * log(p)) * df
  j <- which.min(ebic)
  list(lambda = fit$lambda[j], ebic = ebic[j], df = df[j],
       gamma = gamma, phi = phi)
}

glmnet_qpois_phi <- function(fit_cv, X, y, offset = NULL, s = "lambda.1se",
                             n_extra_df = 0) {
  if (is.null(offset)) offset <- rep(0, length(y))
  mu <- as.numeric(stats::predict(
    fit_cv, newx = X, s = s, newoffset = offset, type = "response"
  ))
  mu <- pmax(mu, 1e-8)
  b <- as.numeric(stats::coef(fit_cv, s = s))[-1L]
  df <- sum(is.finite(b) & b != 0) + 1 + n_extra_df
  phi <- sum((as.numeric(y) - mu)^2 / mu) / max(length(y) - df, 1)
  max(phi, 1e-8)
}

empty_main_index <- function() {
  data.frame(
    Index = integer(0),
    CS = character(0),
    PIP = numeric(0),
    lbf = numeric(0),
    Pvalue = numeric(0)
  )
}

eval_selected_xcs <- function(sel, true_idx, x_cs_id, p) {
  sel <- sort(unique(as.integer(sel)))
  sel <- sel[sel >= 1 & sel <= p]
  x_cs_id <- as.integer(x_cs_id)
  true_idx <- sort(unique(as.integer(true_idx)))
  true_blocks <- sort(unique(x_cs_id[true_idx]))
  sel_blocks <- sort(unique(x_cs_id[sel]))
  block_size <- tabulate(x_cs_id, nbins = max(x_cs_id))
  sel_vars <- if (length(sel_blocks)) which(x_cs_id %in% sel_blocks) else integer(0)

  ok <- sel_blocks %in% true_blocks
  false_cs <- sum(!ok)
  fdr_cs <- if (length(sel_blocks)) mean(!ok) else 0
  power_cs <- if (length(true_blocks)) {
    mean(true_blocks %in% sel_blocks)
  } else {
    NA_real_
  }
  fdr_cs_var <- if (length(sel_vars)) {
    mean(!(x_cs_id[sel_vars] %in% true_blocks))
  } else {
    0
  }
  false_selected <- sum(!(sel %in% true_idx))

  data.frame(
    power_cs = power_cs,
    fdr_cs = fdr_cs,
    false_cs = false_cs,
    type1_cs = as.numeric(false_cs > 0),
    n_cs = length(sel_blocks),
    n_cs_var = length(sel_vars),
    fdr_cs_var = fdr_cs_var,
    type1_cs_var = as.numeric(fdr_cs_var > 0),
    mean_cs_size = if (length(sel_blocks)) mean(block_size[sel_blocks]) else NA_real_,
    min_lbf = NA_real_,
    max_lbf = NA_real_,
    power_pip = NA_real_,
    fdr_pip = NA_real_,
    false_pip = NA_real_,
    type1_pip = NA_real_,
    n_pip = NA_integer_,
    power_selected = if (length(true_idx)) mean(true_idx %in% sel) else NA_real_,
    fdr_selected = if (length(sel)) mean(!(sel %in% true_idx)) else 0,
    false_selected = false_selected,
    type1_selected = as.numeric(false_selected > 0),
    n_selected = length(sel)
  )
}

eval_main_index_xcs <- function(main_index, pip, true_idx, x_cs_id, p,
                                coverage = 0.9) {
  if (is.null(main_index) || !nrow(main_index)) {
    main_index <- empty_main_index()
  }
  main_index <- main_index[main_index$Index >= 1 & main_index$Index <= p, ,
                           drop = FALSE]
  x_cs_id <- as.integer(x_cs_id)
  true_idx <- sort(unique(as.integer(true_idx)))
  true_blocks <- sort(unique(x_cs_id[true_idx]))

  if (nrow(main_index)) {
    cls <- split(main_index$Index, main_index$CS)
    cs_blocks <- lapply(cls, function(ii) sort(unique(x_cs_id[ii])))
    ok <- vapply(cs_blocks, function(ii) any(ii %in% true_blocks), logical(1))
    selected_blocks <- sort(unique(unlist(cs_blocks)))
    selected_vars <- sort(unique(main_index$Index))
    false_cs <- sum(!ok)
    lbf <- numeric(0)
    if ("lbf" %in% names(main_index)) {
      lbf <- unlist(lapply(split(main_index$lbf, main_index$CS), function(x) {
        x <- x[is.finite(x)]
        if (length(x)) x[1] else NA_real_
      }))
      lbf <- lbf[is.finite(lbf)]
    }
    power_cs <- if (length(true_blocks)) {
      mean(true_blocks %in% selected_blocks)
    } else {
      NA_real_
    }
    fdr_cs <- mean(!ok)
    fdr_cs_var <- mean(!(x_cs_id[selected_vars] %in% true_blocks))
    n_cs <- length(cls)
    n_cs_var <- length(selected_vars)
    mean_cs_size <- mean(vapply(cls, length, integer(1)))
    min_lbf <- if (length(lbf)) min(lbf) else NA_real_
    max_lbf <- if (length(lbf)) max(lbf) else NA_real_
  } else {
    power_cs <- if (length(true_blocks)) 0 else NA_real_
    fdr_cs <- 0
    false_cs <- 0
    n_cs <- 0
    n_cs_var <- 0
    fdr_cs_var <- 0
    mean_cs_size <- NA_real_
    min_lbf <- NA_real_
    max_lbf <- NA_real_
  }

  sel <- which(is.finite(pip[seq_len(p)]) & pip[seq_len(p)] > coverage)
  sel_blocks <- sort(unique(x_cs_id[sel]))
  false_pip <- sum(!(sel_blocks %in% true_blocks))
  power_pip <- if (length(true_blocks)) {
    mean(true_blocks %in% sel_blocks)
  } else {
    NA_real_
  }
  fdr_pip <- if (length(sel_blocks)) mean(!(sel_blocks %in% true_blocks)) else 0

  data.frame(
    power_cs = power_cs,
    fdr_cs = fdr_cs,
    false_cs = false_cs,
    type1_cs = as.numeric(false_cs > 0),
    n_cs = n_cs,
    n_cs_var = n_cs_var,
    fdr_cs_var = fdr_cs_var,
    type1_cs_var = as.numeric(fdr_cs_var > 0),
    mean_cs_size = mean_cs_size,
    min_lbf = min_lbf,
    max_lbf = max_lbf,
    power_pip = power_pip,
    fdr_pip = fdr_pip,
    false_pip = false_pip,
    type1_pip = as.numeric(false_pip > 0),
    n_pip = length(sel_blocks),
    power_selected = NA_real_,
    fdr_selected = NA_real_,
    false_selected = NA_real_,
    type1_selected = NA_real_,
    n_selected = NA_integer_
  )
}

mean_finite <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(NA_real_)
  mean(x[is.finite(x)])
}

summarize_sim <- function(sim_df) {
  group_cols <- intersect(c("family", "n", "setting", "case", "structure",
                            "x_structure", "varX_ratio", "varX", "varZ"),
                          names(sim_df))
  rows <- split(sim_df, sim_df[group_cols], drop = TRUE)
  out <- lapply(rows, function(x) {
    ans <- as.data.frame(x[1, group_cols, drop = FALSE])
    ans$n_seed <- nrow(x)
    nm <- setdiff(names(x), c("family", "n", "setting", "case", "structure",
                              "x_structure", "varX_ratio", "varX", "varZ",
                              "seed", "true_signal", "true_cs"))
    for (j in nm) {
      if (is.numeric(x[[j]])) ans[[j]] <- mean_finite(x[[j]])
    }
    ans
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

summarize_result <- function(result_df) {
  group_cols <- intersect(c("family", "n", "setting", "case", "structure",
                            "x_structure", "varX_ratio", "varX", "varZ",
                            "method"),
                          names(result_df))
  rows <- split(result_df, result_df[group_cols], drop = TRUE)
  out <- lapply(rows, function(x) {
    cbind(
      as.data.frame(x[1, group_cols, drop = FALSE]),
      data.frame(
        n_run = nrow(x),
        power_cs = mean_finite(x$power_cs),
        fdr_cs = mean_finite(x$fdr_cs),
        false_cs = mean_finite(x$false_cs),
        type1_cs = mean_finite(x$type1_cs),
        n_cs = mean_finite(x$n_cs),
        n_cs_var = mean_finite(x$n_cs_var),
        fdr_cs_var = mean_finite(x$fdr_cs_var),
        type1_cs_var = mean_finite(x$type1_cs_var),
        mean_cs_size = mean_finite(x$mean_cs_size),
        min_lbf = mean_finite(x$min_lbf),
        max_lbf = mean_finite(x$max_lbf),
        power_pip = mean_finite(x$power_pip),
        fdr_pip = mean_finite(x$fdr_pip),
        false_pip = mean_finite(x$false_pip),
        type1_pip = mean_finite(x$type1_pip),
        n_pip = mean_finite(x$n_pip),
        power_selected = mean_finite(x$power_selected),
        fdr_selected = mean_finite(x$fdr_selected),
        false_selected = mean_finite(x$false_selected),
        type1_selected = mean_finite(x$type1_selected),
        n_selected = mean_finite(x$n_selected),
        time_sec = mean_finite(x$time_sec),
        bic = mean_finite(x$bic),
        qbic = mean_finite(x$qbic),
        qphi = mean_finite(x$qphi),
        ebic = mean_finite(x$ebic),
        ebic_df = mean_finite(x$ebic_df),
        ebic_gamma = mean_finite(x$ebic_gamma),
        theta_hat = mean_finite(x$theta_hat)
      )
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

make_summary <- function(sim_df, result_df) {
  sim_sum <- summarize_sim(sim_df)
  res_sum <- summarize_result(result_df)
  group_cols <- intersect(c("family", "n", "setting", "case", "structure",
                            "x_structure", "varX_ratio", "varX", "varZ"),
                          names(res_sum))
  merge(res_sum, sim_sum, by = group_cols, all.x = TRUE, sort = FALSE)
}

bind_rows_fill <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names)))
  for (i in seq_along(rows)) {
    miss <- setdiff(all_names, names(rows[[i]]))
    for (nm in miss) rows[[i]][[nm]] <- NA
    rows[[i]] <- rows[[i]][, all_names, drop = FALSE]
  }
  do.call(rbind, rows)
}
