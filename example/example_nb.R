suppressPackageStartupMessages({
  library(devtools)
  library(MASS)
  library(susieR)
  library(logisticsusie)
})

load_all(".", quiet = TRUE)
source("example/otherfunction.R")

ar_cov <- function(p, rho) {
  toeplitz(rho ^ (0:(p - 1)))
}

nb_uni_fun <- function(x, y, e, prior_variance,
                       estimate_intercept = 0, ...) {
  v0 <- prior_variance
  fit <- MASS::glm.nb(y ~ x + offset(e), link = "log")
  co <- summary(fit)$coefficients
  bhat <- co["x", "Estimate"]
  s <- co["x", "Std. Error"]
  if (!is.finite(bhat) || !is.finite(s) || s <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  z <- bhat / s
  lbf_wake <- 0.5 * log(s^2 / (v0 + s^2)) +
    0.5 * z^2 * v0 / (v0 + s^2)
  fit0 <- MASS::glm.nb(y ~ 1 + offset(e), link = "log")
  lrt <- as.numeric(2 * (logLik(fit) - logLik(fit0)))
  lbf <- lbf_wake - 0.5 * z^2 + 0.5 * lrt
  v1 <- 1 / (1 / v0 + 1 / s^2)
  mu1 <- v1 * bhat / s^2
  list(mu = mu1, var = v1, lbf = lbf,
       prior_variance = mu1^2 + v1, intercept = 0)
}

simulate_nb_data <- function(seed, n = 1000, p = 10, q = 5,
                             total_h2 = 0.3, theta = 10) {
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
  mu <- exp(0.5 + eta_parts$eta)
  y <- rnbinom(n = n, size = theta, mu = mu)

  list(X = X, Z = Z, y = y, true_idx = true_idx, theta = theta)
}

fit_irls_nb_once <- function(dat, L = 3, L.init = 1) {
  quiet_eval(SuSiE_IRLS(
    X = dat$X, Z = dat$Z, y = dat$y,
    family = "negbin",
    theta_init = dat$theta,
    estimate_theta = TRUE,
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

run_nb_benchmark <- function(n_rep = 10, n = 1000, p = 10,
                             seed0 = 1, L = 3, L.init = 1) {
  rows <- vector("list", n_rep * 2L)
  row_id <- 1L

  for (iter in seq_len(n_rep)) {
    dat <- simulate_nb_data(seed = seed0 + iter - 1L, n = n, p = p)

    t1 <- Sys.time()
    fit_irls <- tryCatch(
      fit_irls_nb_once(dat, L = L, L.init = L.init),
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (inherits(fit_irls, "error")) {
      rows[[row_id]] <- data.frame(
        iter = iter, method = "SuSiE_IRLS_nb",
        power = NA_real_, false_cs = NA_real_, n_cs = NA_integer_,
        time_sec = elapsed, error = fit_irls$message
      )
    } else {
      eval <- cs_contains_truth(fit_irls$main_index, dat$true_idx)
      rows[[row_id]] <- data.frame(
        iter = iter, method = "SuSiE_IRLS_nb",
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
        ser_function = logisticsusie::ser_from_univariate(nb_uni_fun)
      )),
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (inherits(fit_ibss, "error")) {
      rows[[row_id]] <- data.frame(
        iter = iter, method = "IBSS_nb_augmented_Z",
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
        iter = iter, method = "IBSS_nb_augmented_Z",
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
  bench <- run_nb_benchmark(n_rep = 10, n = 1000, p = 10, L.init = 1)
  print(bench$per_run)
  print(bench$summary)
}
