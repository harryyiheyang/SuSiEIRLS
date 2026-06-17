suppressPackageStartupMessages({
  library(devtools)
  library(MASS)
  library(survival)
  library(susieR)
  library(logisticsusie)
})

load_all(".", quiet = TRUE)
source("example/otherfunction.R")

ar_cov <- function(p, rho) {
  toeplitz(rho ^ (0:(p - 1)))
}

surv_uni_fun <- function(x, y, e, prior_variance,
                         estimate_intercept = 0, ...) {
  v0 <- prior_variance
  fit <- survival::coxph(y ~ x + offset(e), ties = "breslow")
  out <- summary(fit)$coefficients
  bhat <- out[1, "coef"]
  s <- out[1, "se(coef)"]
  if (!is.finite(bhat) || !is.finite(s) || s <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  z <- bhat / s
  lbf <- 0.5 * log(s^2 / (v0 + s^2)) +
    0.5 * z^2 * v0 / (v0 + s^2)
  lbf <- lbf - 0.5 * z^2 + 0.5 * summary(fit)$logtest["test"]
  v1 <- 1 / (1 / v0 + 1 / s^2)
  mu1 <- v1 * bhat / s^2
  list(mu = mu1, var = v1, lbf = lbf,
       prior_variance = mu1^2 + v1, intercept = 0)
}

simulate_cox_data <- function(seed, n = 1000, p = 10, q = 5,
                              total_h2 = 0.3, target_censor = 0.3,
                              base_hazard = 0.05) {
  set.seed(seed)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = ar_cov(p, 0.5))
  X <- scale(X)
  colnames(X) <- paste0("SNP", seq_len(p))

  Z <- scale(matrix(rnorm(n * q), n, q))
  colnames(Z) <- paste0("Z", seq_len(q))

  true_idx <- c(2, 5, 8)
  beta <- rep(0, p)
  beta[true_idx] <- c(1, -1, 1)
  alpha <- rnorm(q)
  eta_parts <- make_eta_components(
    X = X, Z = Z, beta = beta, alpha = alpha,
    total_h2 = total_h2, z_to_x_ratio = 2
  )
  eta <- scale_cox_liability_h2(eta_parts$eta, h2 = total_h2)

  Ttrue <- -log(runif(n)) / (base_hazard * exp(eta))
  censor_rate <- exp(stats::uniroot(
    function(log_rate) mean(1 - exp(-exp(log_rate) * Ttrue)) - target_censor,
    interval = c(-20, 20)
  )$root)
  Ctime <- rexp(n, rate = censor_rate)
  time <- pmin(Ttrue, Ctime)
  status <- as.integer(Ttrue <= Ctime)

  list(X = X, Z = Z, y = survival::Surv(time, status),
       time = time, status = status, true_idx = true_idx)
}

fit_irls_cox_once <- function(dat, L = 3, L.init = 1) {
  quiet_eval(SuSiE_IRLS(
    X = dat$X, Z = dat$Z, y = dat$y,
    L = L,
    L.init = L.init,
    max.iter = 5,
    min.iter = 1,
    max.eps = 1e-4,
    susie.iter = 100,
    coverage = 0.95,
    n_threads = 2,
    verbose = FALSE
  ))
}

run_cox_benchmark <- function(n_rep = 10, n = 1000, p = 10,
                              seed0 = 1, L = 3, L.init = 1) {
  rows <- vector("list", n_rep * 2L)
  row_id <- 1L

  for (iter in seq_len(n_rep)) {
    dat <- simulate_cox_data(seed = seed0 + iter - 1L, n = n, p = p)

    t1 <- Sys.time()
    fit_irls <- tryCatch(
      fit_irls_cox_once(dat, L = L, L.init = L.init),
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (inherits(fit_irls, "error")) {
      rows[[row_id]] <- data.frame(
        iter = iter, method = "SuSiE_IRLS_cox",
        power = NA_real_, false_cs = NA_real_, n_cs = NA_integer_,
        time_sec = elapsed, error = fit_irls$message
      )
    } else {
      eval <- cs_contains_truth(fit_irls$main_index, dat$true_idx)
      rows[[row_id]] <- data.frame(
        iter = iter, method = "SuSiE_IRLS_cox",
        power = eval$power, false_cs = eval$false_cs, n_cs = eval$n_cs,
        time_sec = elapsed, error = NA_character_
      )
    }
    row_id <- row_id + 1L

    t1 <- Sys.time()
    fit_ibss <- tryCatch(
      quiet_eval(logisticsusie::ibss_from_ser(
        X = cbind(dat$X, dat$Z), y = dat$y, L = L + ncol(dat$Z),
        tol = 1e-4, maxit = 100, num_cores = 1,
        ser_function = logisticsusie::ser_from_univariate(surv_uni_fun)
      )),
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (inherits(fit_ibss, "error")) {
      rows[[row_id]] <- data.frame(
        iter = iter, method = "IBSS_cox_augmented_Z",
        power = NA_real_, false_cs = NA_real_, n_cs = NA_integer_,
        time_sec = elapsed, error = fit_ibss$message
      )
    } else {
      ibss_index <- susie_to_main_index_x_only(
        fit_ibss, X_aug = cbind(dat$X, dat$Z), p = ncol(dat$X),
        coverage = 0.95, min_abs_cor = 0.1
      )
      eval <- cs_contains_truth(ibss_index, dat$true_idx)
      rows[[row_id]] <- data.frame(
        iter = iter, method = "IBSS_cox_augmented_Z",
        power = eval$power, false_cs = eval$false_cs, n_cs = eval$n_cs,
        time_sec = elapsed, error = NA_character_
      )
    }
    row_id <- row_id + 1L
    message("finished replicate ", iter, "/", n_rep)
  }

  per_run <- do.call(rbind, rows)
  list(per_run = per_run, summary = benchmark_summary(per_run))
}

if (sys.nframe() == 0L) {
  bench <- run_cox_benchmark(n_rep = 10, n = 1000, p = 10, L.init = 1)
  print(bench$per_run)
  print(bench$summary)
}
