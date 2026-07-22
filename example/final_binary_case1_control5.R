require(CppMatrix)
require(glmnet)
require(logisticsusie)
require(susieR)
require(SuSiEIRLS)

source(file.path("example", "evaluation.R"))

family_name <- "binary"
setting <- "case1_control5"
out_structure <- setting
case_prev <- 1 / 6
ns <- c(250L, 500L, 1000L)
nrep <- 500L
irls_threads <- 1L
out_file <- file.path("example", "final_binary_case1_control5_results.rds")

p <- 50L
q <- 10L
L <- 5L
coverage <- 0.9
true_idx <- c(5L, 20L, 35L)
varZ <- 0.2
varX_ratio_vec <- c(1, 0.5)
case_vec <- c("independent", "mediator", "interaction")
rho_x <- 0.5
rho_z <- 0.2
r_m <- 0.25
r_i <- 0.25
ebic_gamma <- 1
glmnet_nlambda <- 50L
glmnet_lambda_min_ratio <- 0.05
glmnet_maxit <- 1000000L

Sx <- stats::toeplitz(rho_x ^ (0:(p - 1L)))
Hx <- chol(Sx)
Sz <- stats::toeplitz(rho_z ^ (0:(q - 1L)))
Hz <- chol(Sz)
ser_logit <- logisticsusie::ser_from_univariate(logit_uni_fun)

rows <- list()
fails <- list()
k <- 0L
h <- 0L

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
        if (attempt > nrep + 200L) stop("Too many failed attempts in binary simulation")
        rep_id <- rep_success + 1L
        seed <- 1100000L + n * 1000L + match(case, case_vec) * 100000L +
          rr * 10000L + attempt
        stage <- "simulate"

        ok <- tryCatch({
          set.seed(seed)
          rows0 <- list()

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
          if (!is.finite(vare) || vare <= 0) stop("1 - var(eta) must be positive in binary simulation.")
          eta_unit <- eta / sqrt(vare)
          eta_out <- eta_unit * sqrt(pi^2 / 3)
          etaZ_out <- etaZ / sqrt(vare) * sqrt(pi^2 / 3)
          latent_total <- var_eta + vare
          pve_main <- stats::var(etaX) / latent_total
          pve_z <- stats::var(etaZ) / latent_total
          pve_int <- stats::var(etaI) / latent_total
          pve_total <- var_eta / latent_total

          a <- stats::uniroot(
            function(a0) mean(stats::plogis(clip_eta(a0 + eta_out))) - case_prev,
            interval = c(-50, 50)
          )$root
          pr <- stats::plogis(clip_eta(a + eta_out))
          y <- stats::rbinom(n, size = 1L, prob = pr)
          observed_case <- mean(y == 1)

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
            evals$rep <- rep_id
            evals$attempt <- attempt
            evals$seed <- seed
            evals$true_signal <- paste(true_idx, collapse = ",")
            evals$pve_main <- pve_main
            evals$pve_z <- pve_z
            evals$pve_int <- pve_int
            evals$pve_total <- pve_total
            evals$var_eta <- stats::var(eta)
            evals$vare <- vare
            evals$observed_case <- observed_case
            rows0[[length(rows0) + 1L]] <<- evals
          }

          X_aug <- cbind(X, Z)
          penalty_factor <- c(rep(1, p), rep(0, q))
          foldid <- sample(rep(seq_len(5L), length.out = n))

          stage <- "glmnet_cv"
          t1 <- proc.time()[["elapsed"]]
          fit_cv <- glmnet::cv.glmnet(
            x = X_aug, y = y, family = "binomial", alpha = 1,
            penalty.factor = penalty_factor, standardize = FALSE,
            foldid = foldid, nlambda = glmnet_nlambda,
            lambda.min.ratio = glmnet_lambda_min_ratio,
            maxit = glmnet_maxit
          )
          j <- match(fit_cv$lambda.1se, fit_cv$lambda)
          if (!is.finite(j)) j <- which.min(abs(fit_cv$lambda - fit_cv$lambda.1se))
          bhat <- as.numeric(fit_cv$glmnet.fit$beta[, j])
          sel <- which(is.finite(bhat[seq_len(p)]) & bhat[seq_len(p)] != 0)
          add_common(eval_selected(sel, true_idx), "glmnet_cv_1se", proc.time()[["elapsed"]] - t1)

          stage <- "glmnet_path"
          t1 <- proc.time()[["elapsed"]]
          fit_path <- glmnet::glmnet(
            x = X_aug, y = y, family = "binomial", alpha = 1,
            penalty.factor = penalty_factor, standardize = FALSE,
            nlambda = glmnet_nlambda,
            lambda.min.ratio = glmnet_lambda_min_ratio,
            maxit = glmnet_maxit
          )
          path_time <- proc.time()[["elapsed"]] - t1

          bb <- glmnet_ebic_lambda(fit_path, n, p, gamma = 0)
          j <- match(bb$lambda, fit_path$lambda)
          if (!is.finite(j)) j <- which.min(abs(fit_path$lambda - bb$lambda))
          bhat <- as.numeric(fit_path$beta[, j])
          sel <- which(is.finite(bhat[seq_len(p)]) & bhat[seq_len(p)] != 0)
          evals <- eval_selected(sel, true_idx)
          evals$bic <- bb$ebic
          evals$ebic_df <- bb$df
          add_common(evals, "glmnet_bic", path_time)

          bb <- glmnet_ebic_lambda(fit_path, n, p, gamma = ebic_gamma)
          j <- match(bb$lambda, fit_path$lambda)
          if (!is.finite(j)) j <- which.min(abs(fit_path$lambda - bb$lambda))
          bhat <- as.numeric(fit_path$beta[, j])
          sel <- which(is.finite(bhat[seq_len(p)]) & bhat[seq_len(p)] != 0)
          evals <- eval_selected(sel, true_idx)
          evals$ebic <- bb$ebic
          evals$ebic_df <- bb$df
          evals$ebic_gamma <- bb$gamma
          add_common(evals, "glmnet_ebic", path_time)

          stage <- "irls"
          t1 <- proc.time()[["elapsed"]]
          fit_irls <- SuSiEIRLS::SuSiE_IRLS(X = X, Z = Z, y = y, family = stats::binomial(link = "logit"),
                          L = L, L.init = 1L, max.iter = 8L, min.iter = 2L, max.eps = 1e-04,
                          n_threads = irls_threads, susie_para = list(max_iter = 300L,
                              coverage = coverage, estimate_residual_variance = TRUE,
                              residual_variance = 0.5, residual_variance_lowerbound = 0.1,
                              residual_variance_upperbound = 1, verbose = FALSE))
          pip <- fit_irls$fitX$pip[seq_len(p)]
          add_common(eval_cs_pip(fit_irls$discovery_summary, pip, true_idx, p, coverage), "irls", proc.time()[["elapsed"]] - t1)

          stage <- "irls_fixed"
          t1 <- proc.time()[["elapsed"]]
          fit_irls <- SuSiEIRLS::SuSiE_IRLS(X = X, Z = Z, y = y, family = stats::binomial(link = "logit"),
                          L = L, L.init = 1L, max.iter = 8L, min.iter = 2L, max.eps = 1e-04,
                          n_threads = irls_threads, susie_para = list(max_iter = 300L,
                              coverage = coverage, estimate_residual_variance = FALSE,
                              residual_variance = 1, verbose = FALSE))
          pip <- fit_irls$fitX$pip[seq_len(p)]
          add_common(eval_cs_pip(fit_irls$discovery_summary, pip, true_idx, p, coverage), "irls_fixed_sigma2_1", proc.time()[["elapsed"]] - t1)

          stage <- "ibss_oracle"
          t1 <- proc.time()[["elapsed"]]
          attr(y, "hat_etaZ") <- etaZ_out
          fit_ibss <- logisticsusie::ibss_from_ser(
            X = X, y = y, L = L, tol = 1e-4, maxit = 100,
            num_cores = 1,
            ser_function = ser_logit
          )
          main_index <- ibss_x_main_index(fit_ibss, X, p, coverage)
          pip <- fit_ibss$pip[seq_len(p)]
          add_common(eval_cs_pip(main_index, pip, true_idx, p, coverage), "ibss_oracle", proc.time()[["elapsed"]] - t1)

          stage <- "ibss_hat_offset"
          t1 <- proc.time()[["elapsed"]]
          etaZ_hat <- fit_logit_hat_eta(y, Z)
          attr(y, "hat_etaZ") <- etaZ_hat
          fit_ibss <- logisticsusie::ibss_from_ser(
            X = X, y = y, L = L, tol = 1e-4, maxit = 100,
            num_cores = 1,
            ser_function = ser_logit
          )
          main_index <- ibss_x_main_index(fit_ibss, X, p, coverage)
          pip <- fit_ibss$pip[seq_len(p)]
          add_common(eval_cs_pip(main_index, pip, true_idx, p, coverage), "ibss_hat_offset", proc.time()[["elapsed"]] - t1)

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
            case = case, varX_ratio = varX_ratio, varX = varX, varZ = varZ,
            n = n, rep = rep_id, attempt = attempt, seed = seed,
            stage = stage, error = conditionMessage(e), stringsAsFactors = FALSE
          )
          if (fail_consecutive >= 10L) stop("Ten consecutive failed attempts in binary simulation")
          FALSE
        })
      }
    }
  }
}

result_df <- bind_rows_fill(rows)
sim_cols <- c("family", "seed", "n", "setting", "structure", "case",
              "varX_ratio", "varX", "varZ", "true_signal", "pve_main",
              "pve_z", "pve_int", "pve_total", "var_eta", "vare", "observed_case")
sim_df <- unique(result_df[, sim_cols, drop = FALSE])
result_drop <- c("pve_main", "pve_z", "pve_int", "pve_total", "var_eta", "vare",
                 "observed_case")
result_df <- result_df[, setdiff(names(result_df), result_drop), drop = FALSE]
summary_df <- make_summary(sim_df, result_df)
fail_df <- if (length(fails)) bind_rows_fill(fails) else data.frame()
settings <- list(
  family = family_name, setting = setting, structure = out_structure,
  ns = ns, nrep = nrep, irls_threads = irls_threads, p = p, q = q,
  L = L, coverage = coverage, true_idx = true_idx, varZ = varZ,
  varX_ratio_vec = varX_ratio_vec, case_prev = case_prev,
  case_vec = case_vec, rho_x = rho_x, rho_z = rho_z,
  r_m = r_m, r_i = r_i, ebic_gamma = ebic_gamma,
  glmnet_nlambda = glmnet_nlambda,
  glmnet_lambda_min_ratio = glmnet_lambda_min_ratio,
  glmnet_maxit = glmnet_maxit
)
saveRDS(list(settings = settings, sim = sim_df, result = result_df,
             summary = summary_df, fail = fail_df), out_file)
print(summary_df)
print(fail_df)
