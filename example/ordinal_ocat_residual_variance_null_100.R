assign(".ordinal_ocat_probe_no_run", TRUE, envir = .GlobalEnv)
source(file.path("example", "ordinal_ocat_probe.R"))

set.seed(97101)

nsim <- 100L
n <- 500L
p <- 60L
q <- 5L
L <- 5L
rho_x <- 0.5
rho_z <- 0.2
varZ <- 0.25
cuts <- c(-1.6, -0.6, 0.4, 1.4)

Sx <- stats::toeplitz(rho_x ^ (0:(p - 1L)))
Hx <- chol(Sx)
Sz <- stats::toeplitz(rho_z ^ (0:(q - 1L)))
Hz <- chol(Sz)

rows <- vector("list", nsim)

for (r in seq_len(nsim)) {
  X <- matrix(stats::rnorm(n * p), n, p)
  X <- CppMatrix::matrixMultiply(X, Hx)
  X <- as.matrix(scale(X))
  colnames(X) <- paste0("X", seq_len(p))
  Z <- matrix(stats::rnorm(n * q), n, q)
  Z <- CppMatrix::matrixMultiply(Z, Hz)
  Z <- as.matrix(scale(Z))
  colnames(Z) <- paste0("Z", seq_len(q))

  a0 <- stats::runif(q, -1, 1)
  etaZ0 <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(a0, ncol = 1)))
  a <- a0 * sqrt(varZ / stats::var(etaZ0))
  eta <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(a, ncol = 1)))
  y_lat <- eta + stats::qlogis(stats::runif(n))
  y_int <- findInterval(y_lat, c(-Inf, cuts, Inf))
  y <- ordered(y_int, levels = seq_len(length(cuts) + 1L))

  Zdf <- as.data.frame(Z)
  fit_null <- ordinal::clm(y ~ ., data = data.frame(y = y, Zdf),
                           link = "logit")
  etaZ_hat <- as.numeric(as.matrix(Zdf) %*% fit_null$beta)

  stat <- ocat_suffstats(
    X = X, y_int = as.integer(y), eta = etaZ_hat, Z = Z,
    alpha = fit_null$alpha, n_threads = 1L, ridge = 1e-6
  )

  fit_est <- susieR::susie_ss(
    XtX = stat$XtX, Xty = stat$Xty, yty = stat$yty, n = n, L = L,
    residual_variance = 1,
    estimate_residual_variance = TRUE,
    residual_variance_lowerbound = 0.05,
    residual_variance_upperbound = 3,
    scaled_prior_variance = 1,
    max_iter = 100L,
    coverage = 0.9,
    estimate_prior_method = "optim"
  )

  zz <- stat$Xty^2 / diag(stat$XtX)

  rows[[r]] <- data.frame(
    seed = r,
    sigma2_susie = fit_est$sigma2,
    z2_mean = mean(zz),
    z2_median = stats::median(zz) / stats::qchisq(0.5, 1),
    z2_q90 = stats::quantile(zz, 0.9) / stats::qchisq(0.9, 1),
    min_pr = stat$min_pr,
    min_h = stat$min_h,
    med_h = stat$med_h,
    max_h = stat$max_h,
    cat_min = min(tabulate(y_int, nbins = length(cuts) + 1L))
  )
}

res <- do.call(rbind, rows)
summary_df <- data.frame(
  metric = c("sigma2_susie", "z2_mean", "z2_median", "z2_q90",
             "min_pr", "min_h", "med_h", "max_h", "cat_min"),
  mean = c(mean(res$sigma2_susie), mean(res$z2_mean),
           mean(res$z2_median), mean(res$z2_q90),
           mean(res$min_pr), mean(res$min_h), mean(res$med_h),
           mean(res$max_h), mean(res$cat_min)),
  sd = c(stats::sd(res$sigma2_susie), stats::sd(res$z2_mean),
         stats::sd(res$z2_median), stats::sd(res$z2_q90),
         stats::sd(res$min_pr), stats::sd(res$min_h), stats::sd(res$med_h),
         stats::sd(res$max_h), stats::sd(res$cat_min)),
  q10 = c(stats::quantile(res$sigma2_susie, 0.1),
          stats::quantile(res$z2_mean, 0.1),
          stats::quantile(res$z2_median, 0.1),
          stats::quantile(res$z2_q90, 0.1),
          stats::quantile(res$min_pr, 0.1),
          stats::quantile(res$min_h, 0.1),
          stats::quantile(res$med_h, 0.1),
          stats::quantile(res$max_h, 0.1),
          stats::quantile(res$cat_min, 0.1)),
  q50 = c(stats::quantile(res$sigma2_susie, 0.5),
          stats::quantile(res$z2_mean, 0.5),
          stats::quantile(res$z2_median, 0.5),
          stats::quantile(res$z2_q90, 0.5),
          stats::quantile(res$min_pr, 0.5),
          stats::quantile(res$min_h, 0.5),
          stats::quantile(res$med_h, 0.5),
          stats::quantile(res$max_h, 0.5),
          stats::quantile(res$cat_min, 0.5)),
  q90 = c(stats::quantile(res$sigma2_susie, 0.9),
          stats::quantile(res$z2_mean, 0.9),
          stats::quantile(res$z2_median, 0.9),
          stats::quantile(res$z2_q90, 0.9),
          stats::quantile(res$min_pr, 0.9),
          stats::quantile(res$min_h, 0.9),
          stats::quantile(res$med_h, 0.9),
          stats::quantile(res$max_h, 0.9),
          stats::quantile(res$cat_min, 0.9))
)

print(summary_df, row.names = FALSE)
saveRDS(
  list(result = res, summary = summary_df),
  file.path("example", "ordinal_ocat_residual_variance_null_100_results.rds")
)
