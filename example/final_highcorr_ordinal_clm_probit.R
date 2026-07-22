require(CppMatrix)
require(logisticsusie)
require(ordinal)
require(ordinalNet)
require(susieR)
require(SuSiEIRLS)

source(file.path("example", "evaluation_highcorr.R"))

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

clm_ser_probit <- function(X, y, o = NULL, prior_variance = 1,
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
  fit0 <- ordinal::clm(y ~ offset(e), data = dat0, link = "probit",
                       threshold = "flexible", start = a0)
  ll0 <- as.numeric(stats::logLik(fit0))
  a1 <- fit0$alpha

  mu <- var <- lbf <- pv <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    dat <- make_clm_data(y = y, X = X[, j, drop = FALSE], e = e)
    fit <- ordinal::clm(y ~ V1 + offset(e), data = dat, link = "probit",
                        threshold = "flexible", start = c(a1, 0))
    co <- summary(fit)$coefficients
    b <- co["V1", "Estimate"]
    s <- co["V1", "Std. Error"]
    if (!is.finite(b) || !is.finite(s) || s <= 0) {
      mu[j] <- 0
      var[j] <- v0
      lbf[j] <- -Inf
      pv[j] <- v0
    } else {
      z <- b / s
      lrt <- 2 * (as.numeric(stats::logLik(fit)) - ll0)
      lbf[j] <- 0.5 * log(s^2 / (v0 + s^2)) +
        0.5 * z^2 * v0 / (v0 + s^2) - 0.5 * z^2 + 0.5 * lrt
      v1 <- 1 / (1 / v0 + 1 / s^2)
      mu[j] <- v1 * b / s^2
      var[j] <- v1
      pv[j] <- mu[j]^2 + var[j]
    }
  }

  m <- max(lbf)
  alpha <- exp(lbf - m)
  alpha <- alpha / sum(alpha)
  lbf_model <- sum(alpha * lbf) - logisticsusie:::categorical_kl(
    alpha, rep(1 / p, p)
  )
  if (isTRUE(estimate_prior_variance)) v0 <- sum(alpha * pv)

  res <- list(mu = mu, var = var, alpha = alpha, intercept = rep(0, p),
              lbf = lbf, lbf_model = lbf_model, prior_variance = v0)
  class(res) <- "ser"
  res
}

fit_clm_hat_eta <- function(y, Z) {
  Zdf <- as.data.frame(Z)
  colnames(Zdf) <- make.names(colnames(Zdf), unique = TRUE)
  fit <- ordinal::clm(y ~ ., data = data.frame(y = y, Zdf,
                                               check.names = FALSE),
                      link = "probit", threshold = "flexible")
  list(eta = as.numeric(as.matrix(Zdf) %*% fit$beta),
       alpha = fit$alpha)
}

ordinalnet_ebic_lambda <- function(fit, n, p, x_names, gamma = 1) {
  if (!is.finite(gamma) || gamma < 0) stop("EBIC gamma must be nonnegative")
  b <- as.matrix(fit$coefs[, x_names, drop = FALSE])
  df <- rowSums(is.finite(b) & b != 0)
  ebic <- -2 * fit$loglik + (log(n) + gamma * log(p)) * df
  j <- which.min(ebic)
  list(index = j, lambda = fit$lambdaVals[j], ebic = ebic[j],
       df = df[j], gamma = gamma, loglik = fit$loglik[j],
       ordinalnet_bic = fit$bic[j])
}

summarize_diagnostics <- function(result_df) {
  group_cols <- intersect(c("family", "n", "setting", "case", "structure",
                            "x_structure", "xz_ratio", "varX_ratio", "varX", "varZ",
                            "method"), names(result_df))
  rows <- split(result_df, result_df[group_cols], drop = TRUE)
  out <- lapply(rows, function(R1) {
    data.frame(
      as.data.frame(R1[1, group_cols, drop = FALSE]),
      iter = mean_finite(R1$iter),
      converged = mean_finite(as.numeric(R1$converged)),
      hit_maxit = mean_finite(as.numeric(R1$hit_maxit)),
      maxit = mean_finite(R1$maxit),
      tol = mean_finite(R1$tol),
      sigma2 = mean_finite(R1$sigma2),
      lambda = mean_finite(R1$lambda),
      ordinalnet_bic = mean_finite(R1$ordinalnet_bic)
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

family_name <- "ocat"
link_name <- "probit"
ns <- c(250L, 500L, 1000L)
nrep <- 500L
irls_threads <- 1L
out_file <- file.path("example", "final_highcorr_ordinal_clm_probit_results.rds")

p <- 50L
q <- 10L
dim_cs <- 5L
dim_ar <- 10L
L <- 5L
coverage <- 0.9
true_idx <- c(5L, 20L, 35L)
varZ <- 0.2
varX_ratio_vec <- c(1, 0.5)
case_vec <- c("independent", "mediator", "interaction")
rho_ar <- 0.3
rho_within <- 0.999
rho_z <- 0.2
x_structure <- "ar1_0p3_by10_ar1_0p999_by5"
r_m <- 0.25
r_i <- 0.25
ebic_gamma <- 1
fit_tol <- 1e-5
irls_max_iter <- 100L
ibss_maxit <- 100L
ordinalnet_maxiter <- 100L

cutoff_probs <- list(
  balanced = c(0.25, 0.25, 0.25, 0.25),
  imbalanced = c(0.15, 0.30, 0.45, 0.10)
)
cutoff_filter <- Sys.getenv("ORD_CUTOFF", "")
if (nzchar(cutoff_filter)) {
  if (!(cutoff_filter %in% names(cutoff_probs))) {
    stop("ORD_CUTOFF must be one of balanced, imbalanced")
  }
  cutoff_probs <- cutoff_probs[cutoff_filter]
  out_file <- file.path(
    "example",
    paste0("final_highcorr_ordinal_clm_probit_", cutoff_filter, "_results.rds")
  )
}

if (p != dim_cs * dim_ar) stop("p must equal dim_cs * dim_ar")
S_within <- stats::toeplitz(rho_within ^ (0:(dim_cs - 1L)))
S_ar <- stats::toeplitz(rho_ar ^ (0:(dim_ar - 1L)))
Sx <- kronecker(S_ar, S_within)
Hx <- chol(Sx)
Sz <- stats::toeplitz(rho_z ^ (0:(q - 1L)))
Hz <- chol(Sz)
x_cs_id <- rep(seq_len(dim_ar), each = dim_cs)
true_cs <- sort(unique(x_cs_id[true_idx]))

rows <- list()
fails <- list()
k <- 0L
h <- 0L

for (cutoff_name in names(cutoff_probs)) {
  setting <- cutoff_name
  out_structure <- setting
  target_cat <- cutoff_probs[[cutoff_name]]
  target_cum <- cumsum(target_cat)[1:3]

  for (n in ns) {
    for (case in case_vec) {
      for (rr in seq_along(varX_ratio_vec)) {
        varX_ratio <- varX_ratio_vec[rr]
        varX <- varZ * varX_ratio
        rep_success <- 0L
        attempt <- 0L
        fail_consecutive <- 0L

        while (rep_success < nrep) {
          attempt <- attempt + 1L
          if (attempt > nrep + 200L) stop("Too many failed attempts in ordinal simulation")
          rep_id <- rep_success + 1L
          seed <- 8100000L + match(cutoff_name, names(cutoff_probs)) * 1000000L +
            n * 1000L + match(case, case_vec) * 100000L +
            rr * 10000L + attempt
          stage <- "simulate"

          ok <- tryCatch({
            set.seed(seed)
            rows0 <- list()

            X <- matrix(stats::rnorm(n * p), n, p)
            X <- CppMatrix::matrixMultiply(X, Hx)
            X <- as.matrix(scale(X))
            colnames(X) <- paste0("X", seq_len(p))

            E <- matrix(stats::rnorm(n * q), n, q)
            E <- CppMatrix::matrixMultiply(E, Hz)
            E <- as.matrix(scale(E))
            Z <- E
            colnames(Z) <- paste0("Z", seq_len(q))

            if (case == "mediator") {
              b1_raw <- c(1, -1, 1)
              b2_raw <- c(-1, 1, -1)
              m1 <- as.numeric(CppMatrix::matrixMultiply(X[, true_idx, drop = FALSE], matrix(b1_raw, ncol = 1)))
              m2 <- as.numeric(CppMatrix::matrixMultiply(X[, true_idx, drop = FALSE], matrix(b2_raw, ncol = 1)))
              b1 <- b1_raw * sqrt((r_m / (1 - r_m)) * stats::var(E[, 1]) / stats::var(m1))
              b2 <- b2_raw * sqrt((r_m / (1 - r_m)) * stats::var(E[, 2]) / stats::var(m2))
              Z[, 1] <- E[, 1] + as.numeric(CppMatrix::matrixMultiply(X[, true_idx, drop = FALSE], matrix(b1, ncol = 1)))
              Z[, 2] <- E[, 2] + as.numeric(CppMatrix::matrixMultiply(X[, true_idx, drop = FALSE], matrix(b2, ncol = 1)))
              Z <- as.matrix(scale(Z))
              colnames(Z) <- paste0("Z", seq_len(q))
            }

            beta0 <- rep(0, p)
            beta0[true_idx] <- 1
            etaX0 <- as.numeric(CppMatrix::matrixMultiply(X, matrix(beta0, ncol = 1)))
            beta <- beta0 * sqrt(varX / stats::var(etaX0))
            etaX <- as.numeric(CppMatrix::matrixMultiply(X, matrix(beta, ncol = 1)))

            alpha0 <- stats::runif(q, -1, 1)
            etaZ0 <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(alpha0, ncol = 1)))
            alpha <- alpha0 * sqrt(varZ / stats::var(etaZ0))
            etaZ <- as.numeric(CppMatrix::matrixMultiply(Z, matrix(alpha, ncol = 1)))

            etaI <- rep(0, n)
            if (case == "interaction") {
              I0 <- X[, true_idx[1]] * Z[, 1] +
                X[, true_idx[2]] * Z[, 2] +
                X[, true_idx[3]] * Z[, 1]
              gamma_i <- sqrt(r_i * stats::var(etaX) / stats::var(I0))
              etaI <- gamma_i * I0
            }

            eta <- etaX + etaZ + etaI
            var_eta <- stats::var(eta)
            vare <- 1 - var_eta
            if (!is.finite(vare) || vare <= 0) stop("1 - var(eta) must be positive in ordinal simulation.")
            etaZ_out <- etaZ / sqrt(vare)
            latent_total <- var_eta + vare
            pve_main <- stats::var(etaX) / latent_total
            pve_z <- stats::var(etaZ) / latent_total
            pve_int <- stats::var(etaI) / latent_total
            pve_total <- var_eta / latent_total

            cuts <- stats::qnorm(target_cum)
            eps <- as.numeric(scale(stats::rnorm(n))) * sqrt(vare)
            y_lat <- eta + eps
            y_int <- findInterval(y_lat, c(-Inf, cuts, Inf))
            y <- ordered(y_int, levels = 1:4)
            ct <- tabulate(y_int, nbins = 4L)
            observed_cat <- ct / n

            add_common <- function(evals, method, time_sec) {
              evals$method <- method
              evals$time_sec <- time_sec
              evals$family <- family_name
              evals$setting <- setting
              evals$structure <- out_structure
              evals$x_structure <- x_structure
              evals$case <- case
              evals$varX_ratio <- varX_ratio
              evals$varX <- varX
              evals$varZ <- varZ
              evals$n <- n
              evals$rep <- rep_id
              evals$attempt <- attempt
              evals$seed <- seed
              evals$true_signal <- paste(true_idx, collapse = ",")
              evals$true_cs <- paste(true_cs, collapse = ",")
              evals$n_true_cs <- length(true_cs)
              evals$pve_main <- pve_main
              evals$pve_z <- pve_z
              evals$pve_int <- pve_int
              evals$pve_total <- pve_total
              evals$var_eta <- stats::var(eta)
              evals$vare <- vare
              evals$link_gen <- link_name
              evals$link_fit <- link_name
              evals$target_cat1 <- target_cat[1]
              evals$target_cat2 <- target_cat[2]
              evals$target_cat3 <- target_cat[3]
              evals$target_cat4 <- target_cat[4]
              evals$observed_cat1 <- observed_cat[1]
              evals$observed_cat2 <- observed_cat[2]
              evals$observed_cat3 <- observed_cat[3]
              evals$observed_cat4 <- observed_cat[4]
              evals$cat_min <- min(ct)
              evals$cut1 <- cuts[1]
              evals$cut2 <- cuts[2]
              evals$cut3 <- cuts[3]
              rows0[[length(rows0) + 1L]] <<- evals
            }

            Xnet <- cbind(Z, X)
            colnames(Xnet) <- c(paste0("Z", seq_len(q)), paste0("X", seq_len(p)))
            x_names <- paste0("X", seq_len(p))

            stage <- "ordinalNet_path"
            t1 <- proc.time()[["elapsed"]]
            fit_onet <- ordinalNet::ordinalNet(
              x = Xnet,
              y = y,
              family = "cumulative",
              link = "probit",
              parallelTerms = TRUE,
              nonparallelTerms = FALSE,
              alpha = 1,
              standardize = FALSE,
              penaltyFactors = c(rep(0, q), rep(1, p)),
              stopThresh = fit_tol,
              threshOut = fit_tol,
              threshIn = fit_tol,
              maxiterOut = ordinalnet_maxiter,
              maxiterIn = ordinalnet_maxiter
            )
            path_time <- proc.time()[["elapsed"]] - t1

            bb <- ordinalnet_ebic_lambda(fit_onet, n, p, x_names, gamma = 0)
            co <- fit_onet$coefs[bb$index, , drop = TRUE]
            bx <- co[x_names]
            sel <- which(is.finite(bx) & abs(bx) > 1e-8)
            evals <- eval_selected_xcs(sel, true_idx, x_cs_id, p)
            evals$bic <- bb$ebic
            evals$ebic_df <- bb$df
            evals$lambda <- bb$lambda
            evals$ordinalnet_bic <- bb$ordinalnet_bic
            evals$iter <- fit_onet$iterOut[bb$index]
            evals$converged <- fit_onet$iterOut[bb$index] < ordinalnet_maxiter
            evals$hit_maxit <- fit_onet$iterOut[bb$index] >= ordinalnet_maxiter
            evals$maxit <- ordinalnet_maxiter
            evals$tol <- fit_tol
            evals$sigma2 <- NA_real_
            add_common(evals, "ordinalNet_bic", path_time)

            bb <- ordinalnet_ebic_lambda(fit_onet, n, p, x_names, gamma = ebic_gamma)
            co <- fit_onet$coefs[bb$index, , drop = TRUE]
            bx <- co[x_names]
            sel <- which(is.finite(bx) & abs(bx) > 1e-8)
            evals <- eval_selected_xcs(sel, true_idx, x_cs_id, p)
            evals$ebic <- bb$ebic
            evals$ebic_df <- bb$df
            evals$ebic_gamma <- bb$gamma
            evals$lambda <- bb$lambda
            evals$ordinalnet_bic <- bb$ordinalnet_bic
            evals$iter <- fit_onet$iterOut[bb$index]
            evals$converged <- fit_onet$iterOut[bb$index] < ordinalnet_maxiter
            evals$hit_maxit <- fit_onet$iterOut[bb$index] >= ordinalnet_maxiter
            evals$maxit <- ordinalnet_maxiter
            evals$tol <- fit_tol
            evals$sigma2 <- NA_real_
            add_common(evals, "ordinalNet_ebic", path_time)

            stage <- "irls"
            t1 <- proc.time()[["elapsed"]]
            fit_irls <- SuSiEIRLS::SuSiE_IRLS(X = X, Z = Z, y = y, family = "clm_probit",
                            L = L, L.init = 1L, max.iter = irls_max_iter, min.iter = 2L,
                            max.eps = fit_tol, n_threads = irls_threads, scale_data = FALSE,
                            susie_para = list(max_iter = 300L, coverage = coverage, estimate_residual_variance = TRUE,
                                residual_variance = 0.5, residual_variance_lowerbound = 0.1,
                                residual_variance_upperbound = 1, verbose = FALSE))
            pip <- fit_irls$fitX$pip[seq_len(p)]
            evals <- eval_main_index_xcs(fit_irls$discovery_summary, pip, true_idx, x_cs_id, p, coverage)
            evals$iter <- fit_irls$diagnostics$iterations
            evals$converged <- fit_irls$diagnostics$eps < fit_tol
            evals$hit_maxit <- fit_irls$diagnostics$iterations >= irls_max_iter &&
              fit_irls$diagnostics$eps >= fit_tol
            evals$maxit <- irls_max_iter
            evals$tol <- fit_tol
            evals$sigma2 <- fit_irls$fitX$sigma2
            evals$lambda <- NA_real_
            add_common(evals, "irls", proc.time()[["elapsed"]] - t1)

            stage <- "irls_fixed"
            t1 <- proc.time()[["elapsed"]]
            fit_irls <- SuSiEIRLS::SuSiE_IRLS(X = X, Z = Z, y = y, family = "clm_probit",
                            L = L, L.init = 1L, max.iter = irls_max_iter, min.iter = 2L,
                            max.eps = fit_tol, n_threads = irls_threads, scale_data = FALSE,
                            susie_para = list(max_iter = 300L, coverage = coverage, estimate_residual_variance = FALSE,
                                residual_variance = 1, verbose = FALSE))
            pip <- fit_irls$fitX$pip[seq_len(p)]
            evals <- eval_main_index_xcs(fit_irls$discovery_summary, pip, true_idx, x_cs_id, p, coverage)
            evals$iter <- fit_irls$diagnostics$iterations
            evals$converged <- fit_irls$diagnostics$eps < fit_tol
            evals$hit_maxit <- fit_irls$diagnostics$iterations >= irls_max_iter &&
              fit_irls$diagnostics$eps >= fit_tol
            evals$maxit <- irls_max_iter
            evals$tol <- fit_tol
            evals$sigma2 <- fit_irls$fitX$sigma2
            evals$lambda <- NA_real_
            add_common(evals, "irls_fixed_sigma2_1", proc.time()[["elapsed"]] - t1)

            stage <- "ibss_oracle"
            t1 <- proc.time()[["elapsed"]]
            y_ibss <- y
            attr(y_ibss, "hat_etaZ") <- etaZ_out
            attr(y_ibss, "clm_alpha") <- cuts / sqrt(vare)
            fit_ibss <- logisticsusie::ibss_from_ser(
              X = X, y = y_ibss, L = L, tol = fit_tol, maxit = ibss_maxit,
              num_cores = 1,
              ser_function = clm_ser_probit
            )
            main_index <- ibss_x_main_index(fit_ibss, X, p, coverage)
            pip <- fit_ibss$pip[seq_len(p)]
            evals <- eval_main_index_xcs(main_index, pip, true_idx, x_cs_id, p, coverage)
            evals$iter <- fit_ibss$iter
            evals$converged <- if (!is.null(fit_ibss$converged)) {
              isTRUE(fit_ibss$converged)
            } else {
              fit_ibss$iter < ibss_maxit
            }
            evals$hit_maxit <- fit_ibss$iter >= ibss_maxit && !isTRUE(evals$converged)
            evals$maxit <- ibss_maxit
            evals$tol <- fit_tol
            evals$sigma2 <- NA_real_
            evals$lambda <- NA_real_
            add_common(evals, "ibss_oracle", proc.time()[["elapsed"]] - t1)

            stage <- "ibss_hat_offset"
            t1 <- proc.time()[["elapsed"]]
            clm_hat <- fit_clm_hat_eta(y, Z)
            z_time <- proc.time()[["elapsed"]] - t1
            t1 <- proc.time()[["elapsed"]]
            y_ibss <- y
            attr(y_ibss, "hat_etaZ") <- clm_hat$eta
            attr(y_ibss, "clm_alpha") <- clm_hat$alpha
            fit_ibss <- logisticsusie::ibss_from_ser(
              X = X, y = y_ibss, L = L, tol = fit_tol, maxit = ibss_maxit,
              num_cores = 1,
              ser_function = clm_ser_probit
            )
            main_index <- ibss_x_main_index(fit_ibss, X, p, coverage)
            pip <- fit_ibss$pip[seq_len(p)]
            evals <- eval_main_index_xcs(main_index, pip, true_idx, x_cs_id, p, coverage)
            evals$iter <- fit_ibss$iter
            evals$converged <- if (!is.null(fit_ibss$converged)) {
              isTRUE(fit_ibss$converged)
            } else {
              fit_ibss$iter < ibss_maxit
            }
            evals$hit_maxit <- fit_ibss$iter >= ibss_maxit && !isTRUE(evals$converged)
            evals$maxit <- ibss_maxit
            evals$tol <- fit_tol
            evals$sigma2 <- NA_real_
            evals$lambda <- NA_real_
            add_common(evals, "ibss_hat_offset", z_time + proc.time()[["elapsed"]] - t1)

            for (ii in seq_along(rows0)) {
              k <- k + 1L
              rows[[k]] <- rows0[[ii]]
            }
            rep_success <- rep_success + 1L
            fail_consecutive <- 0L
            TRUE
          }, error = function(e) {
            fail_consecutive <<- fail_consecutive + 1L
            h <<- h + 1L
            fails[[h]] <<- data.frame(
              family = family_name, setting = setting, structure = out_structure,
              x_structure = x_structure, case = case, varX_ratio = varX_ratio,
              varX = varX, varZ = varZ,
              n = n, rep = rep_id, attempt = attempt, seed = seed,
              stage = stage, error = conditionMessage(e), stringsAsFactors = FALSE
            )
            if (fail_consecutive >= 10L) stop("Ten consecutive failed attempts in ordinal simulation")
            FALSE
          })
        }
      }
    }
  }
}

result_df <- bind_rows_fill(rows)
sim_cols <- c("family", "seed", "n", "setting", "structure", "case",
              "x_structure", "varX_ratio", "varX", "varZ", "true_signal",
              "true_cs", "n_true_cs", "pve_main",
              "pve_z", "pve_int", "pve_total", "var_eta", "vare",
              "link_gen", "link_fit", "target_cat1", "target_cat2",
              "target_cat3", "target_cat4", "observed_cat1",
              "observed_cat2", "observed_cat3", "observed_cat4",
              "cat_min", "cut1", "cut2", "cut3")
sim_df <- unique(result_df[, sim_cols, drop = FALSE])
result_drop <- c("pve_main", "pve_z", "pve_int", "pve_total", "var_eta",
                 "vare", "link_gen", "link_fit", "target_cat1",
                 "target_cat2", "target_cat3", "target_cat4",
                 "observed_cat1", "observed_cat2", "observed_cat3",
                 "observed_cat4", "cat_min", "cut1", "cut2", "cut3")
result_df <- result_df[, setdiff(names(result_df), result_drop), drop = FALSE]
summary_df <- make_summary(sim_df, result_df)
diag_df <- summarize_diagnostics(result_df)
diag_cols <- intersect(c("family", "n", "setting", "case", "structure",
                         "x_structure", "xz_ratio", "varX_ratio", "varX", "varZ",
                         "method"), names(summary_df))
summary_df <- merge(summary_df, diag_df, by = diag_cols,
                    all.x = TRUE, sort = FALSE)
fail_df <- if (length(fails)) bind_rows_fill(fails) else data.frame()
settings <- list(
  family = family_name, link = link_name, cutoff_probs = cutoff_probs,
  ns = ns, nrep = nrep, irls_threads = irls_threads, p = p, q = q,
  L = L, coverage = coverage, true_idx = true_idx, varZ = varZ,
  varX_ratio_vec = varX_ratio_vec, case_vec = case_vec,
  x_structure = x_structure, dim_cs = dim_cs, dim_ar = dim_ar,
  rho_within = rho_within, rho_ar = rho_ar, rho_z = rho_z,
  true_cs = true_cs, r_m = r_m, r_i = r_i,
  ebic_gamma = ebic_gamma, fit_tol = fit_tol,
  irls_max_iter = irls_max_iter, ibss_maxit = ibss_maxit,
  ordinalnet_maxiter = ordinalnet_maxiter
)
saveRDS(list(settings = settings, sim = sim_df, result = result_df,
             summary = summary_df, fail = fail_df), out_file)
print(summary_df)
print(fail_df)
