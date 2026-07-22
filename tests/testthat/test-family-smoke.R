standardize_matrix <- function(X) {
  X <- scale(X)
  storage.mode(X) <- "double"
  X
}

expect_minimal_output <- function(fit) {
  expect_identical(
    names(fit), c("diagnostics", "fitX", "fitJoint", "discovery_summary")
  )
  expect_s3_class(fit$diagnostics, "data.frame")
  expect_equal(nrow(fit$diagnostics), 1)
  expect_identical(
    names(fit$diagnostics), c("iterations", "eps", "runtime_seconds")
  )
  expect_true(is.numeric(fit$fitJoint$n_eff))
  expect_length(fit$fitJoint$n_eff, 1)
}

test_that("public inputs keep package and external-package controls separate", {
  public <- names(formals(SuSiE_IRLS))

  expect_false(any(c(
    "zip_theta", "zip_b", "zip_info", "ridge", "init_cor_method",
    "refit_noncs"
  ) %in% public))
  expect_true(all(c(
    "L.init", "noncs_var", "noncs_max_abs_cor"
  ) %in% public))
})

smoke_args <- function() {
  list(
    L = 2, max.iter = 1, min.iter = 1,
    susie_para = list(max_iter = 20),
    n_threads = 1, verbose = FALSE, scale_data = FALSE,
    noncs_max_abs_cor = 0.9
  )
}

test_that("GLM non-CS runner completes a real smoke fit", {
  set.seed(101)
  n <- 80L
  X <- standardize_matrix(matrix(rnorm(n * 6L), nrow = n))
  Z <- standardize_matrix(matrix(rnorm(n), ncol = 1L))
  y <- rbinom(n, 1L, plogis(0.5 * X[, 1L] + 0.2 * Z[, 1L]))

  fit <- do.call(SuSiE_IRLS, c(list(X = X, Z = Z, y = y), smoke_args()))
  expect_s3_class(fit$fitX, "susie")
  expect_minimal_output(fit)
})

test_that("the gaussian string uses the GLM runner", {
  set.seed(105)
  n <- 80L
  X <- standardize_matrix(matrix(rnorm(n * 6L), nrow = n))
  Z <- standardize_matrix(matrix(rnorm(n), ncol = 1L))
  y <- 0.5 * X[, 1L] + 0.2 * Z[, 1L] + rnorm(n)

  fit <- do.call(SuSiE_IRLS, c(
    list(X = X, Z = Z, y = y, family = "gaussian"), smoke_args()
  ))
  expect_s3_class(fit$fitJoint, "gam")
  expect_minimal_output(fit)
})

test_that("OCAT and CLM are separate runners", {
  expect_true(exists("Run_OCAT", mode = "function"))
  expect_true(exists("Run_CLM", mode = "function"))
  expect_false(identical(body(Run_OCAT), body(Run_CLM)))
  expect_false("clm_link" %in% names(formals(Run_OCAT)))
  expect_true("clm_link" %in% names(formals(Run_CLM)))
})

test_that("Cox non-CS runner completes a real smoke fit", {
  skip_if_not_installed("survival")
  set.seed(102)
  n <- 80L
  X <- standardize_matrix(matrix(rnorm(n * 6L), nrow = n))
  Z <- standardize_matrix(matrix(rnorm(n), ncol = 1L))
  event_time <- rexp(n, rate = exp(0.4 * X[, 1L]))
  censor_time <- rexp(n, rate = 0.7)
  y <- survival::Surv(pmin(event_time, censor_time), event_time <= censor_time)

  fit <- do.call(SuSiE_IRLS, c(
    list(X = X, Z = Z, y = y, family = "cox"), smoke_args()
  ))
  expect_s3_class(fit$fitX, "susie")
  expect_minimal_output(fit)
})

test_that("ordinal non-CS runner completes a real smoke fit", {
  skip_if_not_installed("ordinal")
  set.seed(103)
  n <- 90L
  X <- standardize_matrix(matrix(rnorm(n * 6L), nrow = n))
  Z <- standardize_matrix(matrix(rnorm(n), ncol = 1L))
  latent <- 0.5 * X[, 1L] + 0.2 * Z[, 1L] + rlogis(n)
  y <- ordered(cut(latent, breaks = c(-Inf, -0.5, 0.5, Inf), labels = FALSE))

  fit <- do.call(SuSiE_IRLS, c(
    list(X = X, Z = Z, y = y, family = "clm_logit"), smoke_args()
  ))
  expect_s3_class(fit$fitX, "susie")
  expect_s3_class(fit$fitJoint, "gam")
  expect_minimal_output(fit)
})

test_that("non-logit ordinal input uses the CLM runner", {
  skip_if_not_installed("ordinal")
  set.seed(106)
  n <- 90L
  X <- standardize_matrix(matrix(rnorm(n * 6L), nrow = n))
  Z <- standardize_matrix(matrix(rnorm(n), ncol = 1L))
  latent <- 0.5 * X[, 1L] + 0.2 * Z[, 1L] + rnorm(n)
  y <- ordered(cut(latent, breaks = c(-Inf, -0.5, 0.5, Inf), labels = FALSE))

  fit <- do.call(SuSiE_IRLS, c(
    list(X = X, Z = Z, y = y, family = "clm_probit"), smoke_args()
  ))
  expect_s3_class(fit$fitX, "susie")
  expect_s3_class(fit$fitJoint, "clm")
  expect_minimal_output(fit)
})

test_that("ZIP non-CS runner completes a real smoke fit", {
  skip_if_not_installed("mgcv")
  set.seed(104)
  n <- 100L
  X <- standardize_matrix(matrix(rnorm(n * 6L), nrow = n))
  Z <- standardize_matrix(matrix(rnorm(n), ncol = 1L))
  mu <- exp(0.2 + 0.35 * X[, 1L] + 0.15 * Z[, 1L])
  y <- rpois(n, mu)
  y[rbinom(n, 1L, 0.3) == 1L] <- 0L

  fit <- do.call(SuSiE_IRLS, c(
    list(X = X, Z = Z, y = y, family = mgcv::ziP()), smoke_args()
  ))
  expect_s3_class(fit$fitX, "susie")
  expect_minimal_output(fit)
})
