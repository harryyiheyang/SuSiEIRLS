require(CppMatrix)
require(MASS)
require(ordinal)
require(susieR)
require(logisticsusie)

source(file.path("R", "basicfunction.R"))
source(file.path("example", "evaluation.R"))

set.seed(97007)

clip_eta <- function(x) {
  pmin(pmax(as.numeric(x), -30), 30)
}

make_clm_data <- function(y, X = NULL, e = NULL) {
  dat <- data.frame(y = y)
  if (!is.null(X) && ncol(as.matrix(X)) > 0L) {
    X <- as.matrix(X)
    Xdf <- as.data.frame(X)
    colnames(Xdf) <- paste0("V", seq_len(ncol(X)))
    dat <- cbind(dat, Xdf)
  }
  if (!is.null(e)) dat$e <- as.numeric(e)
  dat
}

ocat_prob_parts <- function(y_int, eta, alpha, eps = 1e-12) {
  n <- length(y_int)
  K <- length(alpha)
  J <- K + 1L
  a <- c(-Inf, alpha, Inf)

  lo <- a[y_int]
  hi <- a[y_int + 1L]
  tl <- lo - eta
  tu <- hi - eta

  Fl <- stats::plogis(tl)
  Fu <- stats::plogis(tu)
  fl <- stats::dlogis(tl)
  fu <- stats::dlogis(tu)
  fpl <- fl * (1 - 2 * Fl)
  fpu <- fu * (1 - 2 * Fu)

  pr <- pmax(Fu - Fl, eps)
  A <- fl - fu
  B <- fpu - fpl
  u_eta <- A / pr
  h_eta <- A^2 / pr^2 - B / pr

  D <- matrix(0, n, K)
  C <- matrix(0, n, K)
  E <- matrix(0, n, K)

  ii <- which(y_int <= K)
  if (length(ii)) {
    jj <- y_int[ii]
    D[cbind(ii, jj)] <- fu[ii]
    C[cbind(ii, jj)] <- -fpu[ii]
    E[cbind(ii, jj)] <- fpu[ii]
  }

  ii <- which(y_int > 1L)
  if (length(ii)) {
    jj <- y_int[ii] - 1L
    D[cbind(ii, jj)] <- D[cbind(ii, jj)] - fl[ii]
    C[cbind(ii, jj)] <- C[cbind(ii, jj)] + fpl[ii]
    E[cbind(ii, jj)] <- E[cbind(ii, jj)] - fpl[ii]
  }

  u_th <- D / pr
  h_eta_th <- sweep(D, 1L, A / pr^2, "*") - sweep(C, 1L, 1 / pr, "*")

  list(
    u_eta = u_eta,
    h_eta = h_eta,
    u_th = u_th,
    h_eta_th = h_eta_th,
    D = D,
    E = E,
    pr = pr
  )
}

ocat_suffstats <- function(X, y_int, eta, Z, alpha,
                           n_threads = 1, ridge = 1e-6,
                           block_size = 10000L) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  K <- length(alpha)

  pp <- ocat_prob_parts(y_int = y_int, eta = eta, alpha = alpha)
  h <- pmax(as.numeric(pp$h_eta), 1e-8)
  sw <- sqrt(h)

  Xh <- X * sw
  XtX <- crossprod(Xh)
  Eh <- matrix(eta * sw, ncol = 1)
  Zh <- Z * sw
  XtE <- CppMatrix::matrixMultiply(Xh, Eh, transA = TRUE)
  XtZ <- CppMatrix::matrixMultiply(Xh, Zh, transA = TRUE)
  EtZ <- CppMatrix::matrixMultiply(Eh, Zh, transA = TRUE)
  ZtZ <- CppMatrix::matrixMultiply(Zh, Zh, transA = TRUE)

  XtU <- as.numeric(CppMatrix::matrixMultiply(
    X, matrix(pp$u_eta, ncol = 1), transA = TRUE
  ))
  ZtU <- as.numeric(CppMatrix::matrixMultiply(
    Z, matrix(pp$u_eta, ncol = 1), transA = TRUE
  ))
  EtU <- sum(eta * pp$u_eta)

  XtT <- CppMatrix::matrixMultiply(X, pp$h_eta_th, transA = TRUE)
  EtT <- CppMatrix::matrixMultiply(matrix(eta, ncol = 1), pp$h_eta_th,
                                   transA = TRUE)
  ZtT <- CppMatrix::matrixMultiply(Z, pp$h_eta_th, transA = TRUE)

  TtT <- matrix(0, K, K)
  for (i in seq_len(n)) {
    di <- pp$D[i, ]
    TtT <- TtT + tcrossprod(di) / pp$pr[i]^2
    TtT <- TtT - diag(pp$E[i, ] / pp$pr[i], K)
  }
  TtU <- colSums(pp$u_th)

  HNN <- rbind(
    cbind(ZtZ, ZtT),
    cbind(t(ZtT), TtT)
  )
  HXN <- cbind(XtZ, XtT)
  HNE <- rbind(t(EtZ), t(EtT))
  UN <- c(ZtU, TtU)

  HNN_inv_HNX <- solve_with_ridge(HNN, t(HXN), ridge = ridge)
  HNN_inv_HNE <- solve_with_ridge(HNN, HNE, ridge = ridge)
  HNN_inv_UN <- solve_with_ridge(HNN, matrix(UN, ncol = 1), ridge = ridge)

  XtX <- XtX - CppMatrix::matrixMultiply(HXN, HNN_inv_HNX)
  XtE <- as.numeric(XtE - CppMatrix::matrixMultiply(HXN, HNN_inv_HNE))
  XtU <- XtU - as.numeric(CppMatrix::matrixMultiply(HXN, HNN_inv_UN))

  XtX <- (XtX + t(XtX)) / 2
  diag(XtX) <- diag(XtX) + ridge

  list(
    XtX = XtX,
    Xty = XtE + XtU,
    yty = n - 1,
    min_pr = min(pp$pr),
    min_h = min(h),
    med_h = stats::median(h),
    max_h = max(h)
  )
}

fit_ocat_score <- function(X, y, Z, L = 5L, max.iter = 5L,
                           susie.iter = 100L, coverage = 0.9,
                           n_threads = 1L, ridge = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  Zdf <- as.data.frame(Z)
  colnames(Zdf) <- paste0("Z", seq_len(ncol(Z)))
  dat <- data.frame(y = y, Zdf)
  fit <- ordinal::clm(y ~ ., data = dat, link = "logit")
  eta <- as.numeric(as.matrix(Zdf) %*% fit$beta)

  beta <- rep(0, p)
  g <- numeric(0)
  XCS <- NULL
  fitX <- NULL
  stat <- NULL

  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    stat <- ocat_suffstats(
      X = X, y_int = as.integer(y), eta = eta, Z = Z,
      alpha = fit$alpha, n_threads = n_threads, ridge = ridge
    )
    fitX <- susieR::susie_ss(
      XtX = stat$XtX, Xty = stat$Xty, yty = stat$yty,
      n = n, L = L, residual_variance = 1,
      estimate_residual_variance = FALSE,
      scaled_prior_variance = 1,
      max_iter = susie.iter, coverage = coverage,
      estimate_prior_method = "optim"
    )
    beta <- clean_coef(stats::coef(fitX)[-1L])
    CSdt <- summary(fitX)$vars
    cs_indices <- sort(unique(CSdt$cs[CSdt$cs > 0]))
    if (!length(cs_indices)) break

    A <- fitX$alpha * 0
    for (i in cs_indices) {
      jj <- CSdt$variable[CSdt$cs == i]
      A[i, jj] <- fitX$alpha[i, jj]
    }
    A <- A * sign(fitX$mu)
    XCS <- CppMatrix::matrixMultiply(X, as.matrix(A), transB = TRUE)
    XCS <- XCS[, cs_indices, drop = FALSE]
    colnames(XCS) <- paste0("Main_CS", cs_indices)

    Xdf <- as.data.frame(cbind(Z, XCS))
    colnames(Xdf) <- c(paste0("Z", seq_len(ncol(Z))),
                       paste0("Main_CS", cs_indices))
    dat <- data.frame(y = y, Xdf)
    fit <- ordinal::clm(y ~ ., data = dat, link = "logit")
    eta <- as.numeric(as.matrix(Xdf) %*% fit$beta)

    g[iter] <- sqrt(mean((beta - beta_prev)^2))
    if (g[iter] < 1e-4 && iter > 2L) break
  }

  list(fitX = fitX, fitJoint = fit, iter = iter, error = g,
       stat = stat, XCS = XCS)
}

clm_ser <- function(X, y, o = NULL, prior_variance = 1,
                    estimate_prior_variance = TRUE,
                    estimate_intercept = TRUE, num_cores = 1, ...) {
  p <- ncol(X)
  if (is.null(o)) o <- rep(0, length(y))
  off <- attr(y, "hat_etaZ")
  if (is.null(off)) off <- rep(0, length(y))
  e <- as.numeric(o) + as.numeric(off)
  a0 <- attr(y, "clm_alpha")
  v0 <- prior_variance

  dat0 <- make_clm_data(y = y, e = e)
  fit0 <- ordinal::clm(y ~ offset(e), data = dat0, link = "logit",
                       start = a0)
  ll0 <- as.numeric(stats::logLik(fit0))
  a1 <- fit0$alpha

  mu <- var <- lbf <- pv <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    dat <- make_clm_data(y = y, X = X[, j, drop = FALSE], e = e)
    fit <- ordinal::clm(y ~ V1 + offset(e), data = dat, link = "logit",
                        start = c(a1, 0))
    co <- summary(fit)$coefficients
    b <- co["V1", "Estimate"]
    s <- co["V1", "Std. Error"]
    z <- b / s
    lrt <- 2 * (as.numeric(stats::logLik(fit)) - ll0)
    lbf[j] <- 0.5 * log(s^2 / (v0 + s^2)) +
      0.5 * z^2 * v0 / (v0 + s^2) - 0.5 * z^2 + 0.5 * lrt
    v1 <- 1 / (1 / v0 + 1 / s^2)
    mu[j] <- v1 * b / s^2
    var[j] <- v1
    pv[j] <- mu[j]^2 + var[j]
  }

  alpha <- exp(lbf - matrixStats::logSumExp(lbf))
  lbf_model <- sum(alpha * lbf) - logisticsusie:::categorical_kl(
    alpha, rep(1 / p, p)
  )
  if (isTRUE(estimate_prior_variance)) v0 <- sum(alpha * pv)

  res <- list(mu = mu, var = var, alpha = alpha, intercept = rep(0, p),
              lbf = lbf, lbf_model = lbf_model, prior_variance = v0)
  class(res) <- "ser"
  res
}

if (!exists(".ordinal_ocat_probe_no_run", envir = .GlobalEnv,
            inherits = FALSE) ||
    !isTRUE(get(".ordinal_ocat_probe_no_run", envir = .GlobalEnv))) {
  nsim <- 10L
  n <- 500L
  p <- 60L
  q <- 5L
  L <- 5L
  coverage <- 0.9
  true_idx <- c(8L, 25L, 43L)
  rho_x <- 0.5
  rho_z <- 0.2
  varX <- 0.35
  varZ <- 0.25
  cuts <- c(-1.6, -0.6, 0.4, 1.4)

  Sx <- stats::toeplitz(rho_x ^ (0:(p - 1L)))
  Hx <- chol(Sx)
  Sz <- stats::toeplitz(rho_z ^ (0:(q - 1L)))
  Hz <- chol(Sz)

  rows <- list()

  for (r in seq_len(nsim)) {
    X <- matrix(stats::rnorm(n * p), n, p)
    X <- CppMatrix::matrixMultiply(X, Hx)
    X <- as.matrix(scale(X))
    colnames(X) <- paste0("X", seq_len(p))
    Z <- matrix(stats::rnorm(n * q), n, q)
    Z <- CppMatrix::matrixMultiply(Z, Hz)
    Z <- as.matrix(scale(Z))
    colnames(Z) <- paste0("Z", seq_len(q))

    b0 <- rep(0, p)
    b0[true_idx] <- c(0.9, -0.8, 0.75)
    etaX0 <- as.numeric(CppMatrix::matrixMultiply(X, matrix(b0, ncol = 1)))
    b <- b0 * sqrt(varX / stats::var(etaX0))
    etaX <- as.numeric(CppMatrix::matrixMultiply(X, matrix(b, ncol = 1)))

    a0 <- stats::runif(q, -1, 1)
    etaZ0 <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(a0, ncol = 1)))
    a <- a0 * sqrt(varZ / stats::var(etaZ0))
    etaZ <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(a, ncol = 1)))
    eta <- etaX + etaZ
    y_lat <- eta + stats::qlogis(stats::runif(n))
    y_int <- findInterval(y_lat, c(-Inf, cuts, Inf))
    y <- ordered(y_int, levels = seq_len(length(cuts) + 1L))

    Zdf <- as.data.frame(Z)
    fit_null <- ordinal::clm(y ~ ., data = data.frame(y = y, Zdf),
                             link = "logit")
    etaZ_hat <- as.numeric(as.matrix(Zdf) %*% fit_null$beta)
    y_ibss <- y
    attr(y_ibss, "hat_etaZ") <- etaZ_hat
    attr(y_ibss, "clm_alpha") <- fit_null$alpha

    t1 <- proc.time()[["elapsed"]]
    fit_score <- fit_ocat_score(
      X = X, y = y, Z = Z, L = L, max.iter = 5L,
      susie.iter = 100L, coverage = coverage, n_threads = 1L
    )
    tm <- proc.time()[["elapsed"]] - t1
    main_index <- Identifying_MainEffect(fit_score$fitX, colnames(X))
    pip <- fit_score$fitX$pip[seq_len(p)]
    ev <- eval_cs_pip(main_index, pip, true_idx, p, coverage)
    ev$method <- "ocat_score"
    ev$seed <- r
    ev$time_sec <- tm
    ev$iter <- fit_score$iter
    ev$min_pr <- fit_score$stat$min_pr
    ev$min_h <- fit_score$stat$min_h
    ev$med_h <- fit_score$stat$med_h
    ev$max_h <- fit_score$stat$max_h
    ev$cat_min <- min(tabulate(y_int, nbins = length(cuts) + 1L))
    rows[[length(rows) + 1L]] <- ev

    t1 <- proc.time()[["elapsed"]]
    fit_ibss <- logisticsusie::ibss_from_ser(
      X = X, y = y_ibss, L = L, tol = 1e-4, maxit = 5L,
      ser_function = clm_ser, num_cores = 1
    )
    tm <- proc.time()[["elapsed"]] - t1
    main_index <- ibss_x_main_index(fit_ibss, X, p, coverage)
    pip <- fit_ibss$pip[seq_len(p)]
    ev <- eval_cs_pip(main_index, pip, true_idx, p, coverage)
    ev$method <- "ibss_clm"
    ev$seed <- r
    ev$time_sec <- tm
    ev$iter <- fit_ibss$iter
    ev$min_pr <- NA_real_
    ev$min_h <- NA_real_
    ev$med_h <- NA_real_
    ev$max_h <- NA_real_
    ev$cat_min <- min(tabulate(y_int, nbins = length(cuts) + 1L))
    rows[[length(rows) + 1L]] <- ev
  }

  res <- do.call(rbind, rows)
  rr <- split(res, res$method)
  sum_res <- lapply(rr, function(R1) {
    data.frame(
      method = R1$method[1],
      power_cs = mean(R1$power_cs, na.rm = TRUE),
      fdr_cs = mean(R1$fdr_cs, na.rm = TRUE),
      n_cs = mean(R1$n_cs, na.rm = TRUE),
      n_cs_var = mean(R1$n_cs_var, na.rm = TRUE),
      power_pip = mean(R1$power_pip, na.rm = TRUE),
      fdr_pip = mean(R1$fdr_pip, na.rm = TRUE),
      n_pip = mean(R1$n_pip, na.rm = TRUE),
      time_sec = mean(R1$time_sec, na.rm = TRUE),
      iter = mean(R1$iter, na.rm = TRUE),
      min_pr = mean(R1$min_pr, na.rm = TRUE),
      min_h = mean(R1$min_h, na.rm = TRUE),
      med_h = mean(R1$med_h, na.rm = TRUE),
      max_h = mean(R1$max_h, na.rm = TRUE),
      cat_min = mean(R1$cat_min, na.rm = TRUE)
    )
  })
  sum_res <- do.call(rbind, sum_res)
  row.names(sum_res) <- NULL

  print(res[, c("seed", "method", "power_cs", "fdr_cs", "n_cs",
                "power_pip", "fdr_pip", "n_pip", "time_sec",
                "iter", "min_pr", "min_h", "med_h", "max_h", "cat_min")],
        row.names = FALSE)
  print(sum_res, row.names = FALSE)

  saveRDS(list(result = res, summary = sum_res),
          file.path("example", "ordinal_ocat_probe_results.rds"))
}
