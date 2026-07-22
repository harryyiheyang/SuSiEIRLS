correlated_pair <- function(r, n = 101L) {
  x <- seq_len(n) - mean(seq_len(n))
  z <- (seq_len(n) - mean(seq_len(n)))^2
  z <- z - mean(z)
  z <- z - x * sum(x * z) / sum(x^2)
  x <- x / sqrt(sum(x^2))
  z <- z / sqrt(sum(z^2))
  cbind(x, r * x + sqrt(1 - r^2) * z)
}

test_that("non-CS correlation threshold has a hard 0.9 cap", {
  xy <- correlated_pair(0.95)

  expect_false(SuSiEIRLS:::noncs_correlation_ok(
    xy[, 2L], xy[, 1L, drop = FALSE], max_abs_cor = 1
  ))
  expect_equal(SuSiEIRLS:::validate_noncs_max_abs_cor(1), 0.9)
})

test_that("correlations at the threshold are rejected", {
  xy <- correlated_pair(0.9)

  expect_false(SuSiEIRLS:::noncs_correlation_ok(
    xy[, 2L], xy[, 1L, drop = FALSE], max_abs_cor = 0.9
  ))
})

test_that("stricter thresholds remain effective", {
  xy <- correlated_pair(0.6)

  expect_false(SuSiEIRLS:::noncs_correlation_ok(
    xy[, 2L], xy[, 1L, drop = FALSE], max_abs_cor = 0.5
  ))
  expect_true(SuSiEIRLS:::noncs_correlation_ok(
    xy[, 2L], xy[, 1L, drop = FALSE], max_abs_cor = 0.7
  ))
})

test_that("no-CS terms are gated against existing refit covariates", {
  xy <- correlated_pair(0.95)
  X <- xy[, 2L, drop = FALSE]
  fit <- structure(
    list(
      alpha = matrix(1, nrow = 1L, ncol = 1L),
      mu = matrix(1, nrow = 1L, ncol = 1L),
      X_column_scale_factors = 1,
      intercept = 0,
      theta = NULL
    ),
    class = "susie"
  )

  expect_null(SuSiEIRLS:::build_no_cs_noncs_refit_term(
    X, fit, cor_design = xy[, 1L, drop = FALSE],
    noncs_max_abs_cor = 1
  ))
})

test_that("projected non-CS terms are gated against Z as well as XCS", {
  xy <- correlated_pair(0)
  X <- cbind(xy[, 1L], xy[, 2L])
  fit <- structure(
    list(
      alpha = matrix(1, nrow = 1L, ncol = 2L),
      mu = matrix(c(1, 1), nrow = 1L),
      X_column_scale_factors = c(1, 1),
      intercept = 0,
      theta = NULL
    ),
    class = "susie"
  )

  expect_null(SuSiEIRLS:::build_noncs_refit_term(
    X = X, fitX = fit, CSdt = data.frame(variable = 1L, cs = 1L),
    cs_indices = 1L, XCS = X[, 1L, drop = FALSE], noncs_var = 0,
    noncs_max_abs_cor = 0.9, cor_design = X[, 2L, drop = FALSE]
  ))
})

test_that("all non-CS runners expose the correlation gate", {
  runners <- c("Run_GLM", "Run_Cox", "Run_OCAT", "Run_ZIP")
  for (runner in runners) {
    fn <- getFromNamespace(runner, "SuSiEIRLS")
    expect_true("noncs_max_abs_cor" %in% names(formals(fn)), info = runner)
    expect_equal(formals(fn)$noncs_max_abs_cor, 0.9, info = runner)
  }
})

test_that("the public correlation threshold is validated", {
  X <- matrix(seq_len(20), nrow = 10L, ncol = 2L)
  y <- rep(c(0, 1), 5L)

  expect_error(
    SuSiE_IRLS(X, y = y, noncs_max_abs_cor = 0),
    "positive finite numeric scalar"
  )
  expect_error(
    SuSiE_IRLS(X, y = y, noncs_max_abs_cor = NA_real_),
    "positive finite numeric scalar"
  )
})
