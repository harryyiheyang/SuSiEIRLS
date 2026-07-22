test_that("Gaussian fit stops on the third consecutive no-CS iteration", {
  set.seed(20260719)
  n <- 100L
  p <- 6L
  Q <- qr.Q(qr(cbind(Intercept = 1, matrix(rnorm(n * (p + 1L)), nrow = n))))
  X <- Q[, 2:(p + 1L), drop = FALSE] * sqrt(n)
  y <- Q[, p + 2L] * sqrt(n)
  colnames(X) <- paste0("X", seq_len(p))

  fit <- SuSiE_IRLS(
    X = X, y = y, family = "gaussian", scale_data = FALSE,
    L = 2L,
    susie_para = list(coverage = 0.99, min_abs_corr = 1, max_iter = 20L),
    max.iter = 8L, min.iter = 6L, max.eps = 0,
    n_threads = 1L, verbose = FALSE
  )

  expect_equal(fit$diagnostics$iterations, 3L)
  expect_null(fit$discovery_summary)
  expect_s3_class(fit$fitX, "susie")
  expect_false(is.null(fit$fitJoint))
})
