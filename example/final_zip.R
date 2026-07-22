require(CppMatrix)
require(logisticsusie)
require(mgcv)
require(SuSiEIRLS)

source(file.path("example", "evaluation.R"))

zip_log1mexp <- function(a) {
  a <- pmax(as.numeric(a), .Machine$double.xmin)
  out <- numeric(length(a))
  ii <- a <= log(2)
  out[ii] <- log(-expm1(-a[ii]))
  out[!ii] <- log1p(-exp(-a[!ii]))
  out
}

zip_ll <- function(y, gamma, theta, b = 0) {
  gamma <- clip_eta(gamma)
  slope <- b + exp(theta[2])
  rho <- clip_eta(theta[1] + slope * gamma)
  et <- exp(rho)
  lam <- exp(gamma)
  z <- y == 0
  ll <- numeric(length(y))
  ll[z] <- -et[z]
  ll[!z] <- zip_log1mexp(et[!z]) + y[!z] * gamma[!z] -
    (lam[!z] + zip_log1mexp(lam[!z])) - lgamma(y[!z] + 1)
  ll
}

zip_rzip <- function(gamma, theta, b = 0) {
  y <- gamma
  n <- length(y)
  lambda <- exp(gamma)
  rho <- theta[1] + (b + exp(theta[2])) * gamma
  pr <- 1 - exp(-exp(rho))
  ind <- pr > stats::runif(n)
  y[!ind] <- 0
  np <- sum(ind)
  y[ind] <- stats::qpois(stats::runif(np, stats::dpois(0, lambda[ind]), 1),
                         lambda[ind])
  y
}

zip_uni_fun <- function(x, y, e, prior_variance,
                        estimate_intercept = 0, ...) {
  v0 <- prior_variance
  theta <- attr(y, "theta")
  bzip <- attr(y, "zip_b")
  if (is.null(theta) || length(theta) != 2L || any(!is.finite(theta))) {
    stop("ziP theta attribute is missing or invalid")
  }
  if (is.null(bzip)) bzip <- 0
  if (!is.finite(bzip) || bzip < 0) stop("ziP b attribute is invalid")
  e <- as.numeric(e)
  x <- as.numeric(x)
  y0 <- as.numeric(y)

  fn1 <- function(par) {
    -sum(zip_ll(y0, e + par[1] + par[2] * x, theta = theta, b = bzip))
  }
  fn0 <- function(a) {
    -sum(zip_ll(y0, e + a[1], theta = theta, b = bzip))
  }

  op0 <- stats::optim(c(0), fn0, method = "BFGS")
  op1 <- stats::optim(c(0, 0), fn1, method = "BFGS", hessian = TRUE)
  if (op0$convergence != 0 || op1$convergence != 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  H <- op1$hessian
  if (any(!is.finite(H))) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  evals <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(evals) <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  V <- CppMatrix::matrixInverse(H)
  b1 <- op1$par[2]
  s2 <- V[2, 2]
  if (!is.finite(b1) || !is.finite(s2) || s2 <= 0) {
    return(list(mu = 0, var = v0, lbf = -Inf,
                prior_variance = v0, intercept = 0))
  }
  s <- sqrt(s2)
  z <- b1 / s
  lrt <- 2 * (-op1$value + op0$value)
  lbf <- 0.5 * log(s2 / (v0 + s2)) +
    0.5 * z^2 * v0 / (v0 + s2)
  lbf <- lbf - 0.5 * z^2 + 0.5 * lrt
  v1 <- 1 / (1 / v0 + 1 / s2)
  m1 <- v1 * b1 / s2
  list(mu = m1, var = v1, lbf = lbf,
       prior_variance = m1^2 + v1, intercept = op1$par[1])
}

zip_fit_z <- function(y, Z) {
  Zdf <- as.data.frame(Z)
  colnames(Zdf) <- paste0("Z", seq_len(ncol(Zdf)))
  rhs <- paste(colnames(Zdf), collapse = " + ")
  fml <- stats::as.formula(paste("y ~", rhs))
  fit <- mgcv::gam(fml, data = data.frame(y = y, Zdf, check.names = FALSE),
                   family = mgcv::ziP(), method = "REML")
  list(eta = as.numeric(fit$linear.predictors),
       theta = as.numeric(fit$family$getTheta(FALSE)))
}

family_name <- "zip"
setting <- "pos30_rm10"
out_structure <- "traditional"
out_file <- file.path("example", "final_zip_results.rds")

ns_env <- Sys.getenv("ZIP_NS", "")
if (nzchar(ns_env)) {
  ns <- as.integer(strsplit(ns_env, ",", fixed = TRUE)[[1]])
} else {
  ns <- c(250L, 500L, 1000L)
}
nrep <- as.integer(Sys.getenv("ZIP_NREP", "20"))
irls_threads <- 1L

p <- 50L
q <- 10L
L <- 5L
coverage <- 0.9
true_idx <- c(5L, 20L, 35L)
varZ <- 0.2
varX_ratio_env <- Sys.getenv("ZIP_VARX_RATIO", "")
if (nzchar(varX_ratio_env)) {
  varX_ratio_vec <- as.numeric(strsplit(varX_ratio_env, ",", fixed = TRUE)[[1]])
} else {
  varX_ratio_vec <- c(1, 0.5)
}
case_all <- c("independent", "mediator", "interaction")
case_vec <- case_all
case_filter <- Sys.getenv("ZIP_CASE", "")
if (nzchar(case_filter)) {
  if (!(case_filter %in% case_all)) {
    stop("ZIP_CASE must be one of independent, mediator, interaction")
  }
  case_vec <- case_filter
  out_file <- file.path("example", paste0("final_zip_", case_filter,
                                          "_results.rds"))
}
rho_x <- 0.5
rho_z <- 0.2
r_m <- as.numeric(Sys.getenv("ZIP_RM", "0.10"))
r_i <- as.numeric(Sys.getenv("ZIP_RI", "0.25"))
theta_zip <- c(-1, 0.2)
zip_b <- 0
target_pos <- 0.30
fit_tol <- 1e-4
ibss_maxit <- as.integer(Sys.getenv("ZIP_IBSS_MAXIT", "8"))
susie_iter <- 100L
irls_max_iter <- 6L

Sx <- stats::toeplitz(rho_x ^ (0:(p - 1L)))
Hx <- chol(Sx)
Sz <- stats::toeplitz(rho_z ^ (0:(q - 1L)))
Hz <- chol(Sz)
ser_zip <- logisticsusie::ser_from_univariate(zip_uni_fun)

rows <- list()
k <- 0L

for (n in ns) {
  for (case in case_vec) {
    for (rr in seq_along(varX_ratio_vec)) {
      varX_ratio <- varX_ratio_vec[rr]
      varX <- varZ * varX_ratio

      for (rep in seq_len(nrep)) {
        seed <- 8300000L + n * 1000L + match(case, case_all) * 100000L +
          rr * 10000L + rep
        set.seed(seed)

        X <- matrix(stats::rnorm(n * p), n, p)
        X <- CppMatrix::matrixMultiply(X, Hx)
        X <- as.matrix(scale(X))
        X <- X[, sample.int(p), drop = FALSE]
        colnames(X) <- paste0("X", seq_len(p))

        E <- matrix(stats::rnorm(n * q), n, q)
        E <- CppMatrix::matrixMultiply(E, Hz)
        E <- as.matrix(scale(E))
        E <- E[, sample.int(q), drop = FALSE]
        Z <- E
        colnames(Z) <- paste0("Z", seq_len(q))

        if (case == "mediator") {
          b1_raw <- c(1, -1, 1)
          b2_raw <- c(-1, 1, -1)
          m1 <- as.numeric(CppMatrix::matrixMultiply(
            X[, true_idx, drop = FALSE], matrix(b1_raw, ncol = 1)
          ))
          m2 <- as.numeric(CppMatrix::matrixMultiply(
            X[, true_idx, drop = FALSE], matrix(b2_raw, ncol = 1)
          ))
          b1 <- b1_raw * sqrt((r_m / (1 - r_m)) * stats::var(E[, 1]) /
                                stats::var(m1))
          b2 <- b2_raw * sqrt((r_m / (1 - r_m)) * stats::var(E[, 2]) /
                                stats::var(m2))
          Z[, 1] <- E[, 1] + as.numeric(CppMatrix::matrixMultiply(
            X[, true_idx, drop = FALSE], matrix(b1, ncol = 1)
          ))
          Z[, 2] <- E[, 2] + as.numeric(CppMatrix::matrixMultiply(
            X[, true_idx, drop = FALSE], matrix(b2, ncol = 1)
          ))
          Z <- as.matrix(scale(Z))
          colnames(Z) <- paste0("Z", seq_len(q))
        }

        beta0 <- rep(0, p)
        beta0[true_idx] <- 1
        etaX0 <- as.numeric(CppMatrix::matrixMultiply(X, matrix(beta0,
                                                               ncol = 1)))
        beta <- beta0 * sqrt(varX / stats::var(etaX0))
        etaX <- as.numeric(CppMatrix::matrixMultiply(X, matrix(beta,
                                                               ncol = 1)))

        alpha0 <- stats::runif(q, -1, 1)
        etaZ0 <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(alpha0,
                                                               ncol = 1)))
        alpha <- alpha0 * sqrt(varZ / stats::var(etaZ0))
        etaZ <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(alpha,
                                                               ncol = 1)))

        etaI <- rep(0, n)
        if (case == "interaction") {
          I0 <- X[, true_idx[1]] * Z[, 1] +
            X[, true_idx[2]] * Z[, 2] +
            X[, true_idx[3]] * Z[, 1]
          gamma_i <- sqrt(r_i * stats::var(etaX) / stats::var(I0))
          etaI <- gamma_i * I0
        }

        eta <- etaX + etaZ + etaI
        pve_main <- stats::var(etaX) / stats::var(eta)
        pve_z <- stats::var(etaZ) / stats::var(eta)
        pve_int <- stats::var(etaI) / stats::var(eta)
        max_abs_xz <- max(abs(stats::cor(X, Z)))
        max_abs_true_xz <- max(abs(stats::cor(X[, true_idx, drop = FALSE], Z)))
        slope <- zip_b + exp(theta_zip[2])
        a <- stats::uniroot(
          function(a0) mean(1 - exp(-exp(theta_zip[1] + slope *
                                           (a0 + eta)))) - target_pos,
          interval = c(-20, 20)
        )$root
        gamma <- a + eta
        y <- zip_rzip(gamma, theta = theta_zip, b = zip_b)
        pos_rate <- mean(y > 0)
        mean_y <- mean(y)

        add_common <- function(evals, method, time_sec) {
          evals$method <- method
          evals$time_sec <- time_sec
          evals$family <- family_name
          evals$setting <- setting
          evals$structure <- out_structure
          evals$case <- case
          evals$varX_ratio <- varX_ratio
          evals$varX <- varX
          evals$varZ <- varZ
          evals$n <- n
          evals$rep <- rep
          evals$seed <- seed
          evals$true_signal <- paste(true_idx, collapse = ",")
          evals$target_pos <- target_pos
          evals$pos_rate <- pos_rate
          evals$mean_y <- mean_y
          evals$pve_main <- pve_main
          evals$pve_z <- pve_z
          evals$pve_int <- pve_int
          evals$max_abs_xz <- max_abs_xz
          evals$max_abs_true_xz <- max_abs_true_xz
          rows[[length(rows) + 1L]] <<- evals
        }

        t1 <- proc.time()[["elapsed"]]
        fit_irls <- SuSiEIRLS::SuSiE_IRLS(X = X, y = y, Z = Z, family = mgcv::ziP(),
                        L = L, max.iter = irls_max_iter, min.iter = 2L, max.eps = 1e-04,
                        n_threads = irls_threads, L.init = 1L, suff_block_size = 5000L,
                        susie_para = list(max_iter = susie_iter, coverage = coverage,
                            estimate_residual_variance = TRUE, residual_variance = 0.5,
                            residual_variance_lowerbound = 0.1, residual_variance_upperbound = 1,
                            verbose = FALSE))
        pip <- fit_irls$fitX$pip[seq_len(p)]
        evals <- eval_cs_pip(fit_irls$discovery_summary, pip, true_idx, p, coverage)
        evals$pip_true_mean <- mean(pip[true_idx])
        evals$pip_true_min <- min(pip[true_idx])
        evals$theta1_hat <- fit_irls$theta_raw[1]
        evals$theta2_hat <- fit_irls$theta_raw[2]
        evals$iter <- NA_integer_
        evals$converged <- NA
        evals$hit_maxit <- NA
        evals$maxit <- irls_max_iter
        evals$tol <- NA_real_
        add_common(evals, "zip_irls", proc.time()[["elapsed"]] - t1)

        y_ibss <- y
        attr(y_ibss, "theta") <- theta_zip
        attr(y_ibss, "zip_b") <- zip_b
        t1 <- proc.time()[["elapsed"]]
        fit_ibss <- logisticsusie::ibss_from_ser(
          X = X, y = y_ibss, L = L, tol = fit_tol, maxit = ibss_maxit,
          num_cores = 1L, ser_function = offset_ser(ser_zip, a + etaZ)
        )
        pip <- fit_ibss$pip[seq_len(p)]
        main_index <- ibss_x_main_index(fit_ibss, X, p, coverage)
        evals <- eval_cs_pip(main_index, pip, true_idx, p, coverage)
        evals$pip_true_mean <- mean(pip[true_idx])
        evals$pip_true_min <- min(pip[true_idx])
        evals$theta1_hat <- theta_zip[1]
        evals$theta2_hat <- theta_zip[2]
        evals$iter <- if (!is.null(fit_ibss$iter)) fit_ibss$iter else NA_integer_
        evals$converged <- if (!is.null(fit_ibss$converged)) {
          isTRUE(fit_ibss$converged)
        } else {
          is.finite(evals$iter) && evals$iter < ibss_maxit
        }
        evals$hit_maxit <- is.finite(evals$iter) && evals$iter >= ibss_maxit &&
          !isTRUE(evals$converged)
        evals$maxit <- ibss_maxit
        evals$tol <- fit_tol
        add_common(evals, "zip_ibss_oracle", proc.time()[["elapsed"]] - t1)

        fit_z <- zip_fit_z(y, Z)
        y_ibss <- y
        attr(y_ibss, "theta") <- fit_z$theta
        attr(y_ibss, "zip_b") <- zip_b
        t1 <- proc.time()[["elapsed"]]
        fit_ibss <- logisticsusie::ibss_from_ser(
          X = X, y = y_ibss, L = L, tol = fit_tol, maxit = ibss_maxit,
          num_cores = 1L, ser_function = offset_ser(ser_zip, fit_z$eta)
        )
        pip <- fit_ibss$pip[seq_len(p)]
        main_index <- ibss_x_main_index(fit_ibss, X, p, coverage)
        evals <- eval_cs_pip(main_index, pip, true_idx, p, coverage)
        evals$pip_true_mean <- mean(pip[true_idx])
        evals$pip_true_min <- min(pip[true_idx])
        evals$theta1_hat <- fit_z$theta[1]
        evals$theta2_hat <- fit_z$theta[2]
        evals$iter <- if (!is.null(fit_ibss$iter)) fit_ibss$iter else NA_integer_
        evals$converged <- if (!is.null(fit_ibss$converged)) {
          isTRUE(fit_ibss$converged)
        } else {
          is.finite(evals$iter) && evals$iter < ibss_maxit
        }
        evals$hit_maxit <- is.finite(evals$iter) && evals$iter >= ibss_maxit &&
          !isTRUE(evals$converged)
        evals$maxit <- ibss_maxit
        evals$tol <- fit_tol
        add_common(evals, "zip_ibss_regular", proc.time()[["elapsed"]] - t1)
      }
    }
  }
}

result_df <- bind_rows_fill(rows)
sim_cols <- c("family", "seed", "n", "setting", "structure", "case",
              "varX_ratio", "varX", "varZ", "true_signal", "target_pos",
              "pos_rate", "mean_y", "pve_main", "pve_z", "pve_int",
              "max_abs_xz", "max_abs_true_xz")
sim_df <- unique(result_df[, sim_cols, drop = FALSE])
result_drop <- c("true_signal", "target_pos", "pos_rate", "mean_y",
                 "pve_main", "pve_z", "pve_int", "max_abs_xz",
                 "max_abs_true_xz")
result_df <- result_df[, setdiff(names(result_df), result_drop), drop = FALSE]
summary_df <- make_summary(sim_df, result_df)

group_cols <- intersect(c("family", "n", "setting", "case", "structure",
                          "varX_ratio", "varX", "varZ", "method"),
                        names(result_df))
full_df <- aggregate(
  I(power_cs == 1) ~ .,
  data = result_df[, c(group_cols, "power_cs"), drop = FALSE],
  FUN = mean
)
names(full_df)[ncol(full_df)] <- "full_power_rate"
summary_df <- merge(summary_df, full_df, by = group_cols,
                    all.x = TRUE, sort = FALSE)
diag_df <- aggregate(
  cbind(iter, converged = as.numeric(converged),
        hit_maxit = as.numeric(hit_maxit), maxit, tol) ~ .,
  data = result_df[, c(group_cols, "iter", "converged", "hit_maxit",
                       "maxit", "tol"), drop = FALSE],
  FUN = mean_finite
)
summary_df <- merge(summary_df, diag_df, by = group_cols,
                    all.x = TRUE, sort = FALSE)
fail_df <- data.frame()
settings <- list(
  family = family_name, setting = setting, structure = out_structure,
  ns = ns, nrep = nrep, irls_threads = irls_threads, p = p, q = q,
  L = L, coverage = coverage, true_idx = true_idx, varZ = varZ,
  varX_ratio_vec = varX_ratio_vec, case_vec = case_vec,
  rho_x = rho_x, rho_z = rho_z, r_m = r_m, r_i = r_i,
  theta_zip = theta_zip, zip_b = zip_b, target_pos = target_pos,
  fit_tol = fit_tol, ibss_maxit = ibss_maxit,
  irls_max_iter = irls_max_iter, susie_iter = susie_iter
)
saveRDS(list(settings = settings, sim = sim_df, result = result_df,
             summary = summary_df, fail = fail_df), out_file)
print(summary_df)
print(fail_df)
