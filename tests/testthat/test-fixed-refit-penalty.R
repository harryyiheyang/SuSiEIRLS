test_that("fixed mgcv ridge uses 1/V and leaves Z unpenalized", {
  set.seed(20260717)
  n <- 160L
  z <- stats::rnorm(n)
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1L, stats::plogis(0.3 * z + 0.8 * x))
  dat <- data.frame(y = y, Z1 = z, Main_CS1 = x)

  b <- SuSiEIRLS:::.mgcv_fit_fixed_ridge(
    "y", c("Z1", "Main_CS1"), dat, stats::binomial(),
    penalty_V = c(Main_CS1 = 0.4), mgcv_model = "gam"
  )
  pen <- attr(b, "refit_penalty")

  expect_s3_class(b, "gam")
  expect_equal(pen$V, c(Main_CS1 = 0.4))
  expect_equal(pen$precision, c(Main_CS1 = 2.5))
  expect_equal(pen$sp, 1)
  expect_equal(unname(b$full.sp), 1)
  expect_true("Z1" %in% names(stats::coef(b)))
  expect_false("Z1" %in% names(pen$V))
  expect_true("Main_CS1" %in% names(stats::coef(b)))
  expect_true(all(is.finite(stats::predict(b, type = "link"))))
})

test_that("Cox ridge uses theta 1/V and leaves Z unpenalized", {
  skip_if_not_installed("survival")
  set.seed(20260718)
  n <- 180L
  z <- stats::rnorm(n)
  x <- stats::rnorm(n)
  te <- stats::rexp(n, exp(0.25 * z + 0.7 * x))
  tc <- stats::rexp(n, 0.5)

  fit <- SuSiEIRLS:::.cox_fit_fixed_ridge(
    pmin(te, tc), te <= tc, matrix(z, ncol = 1L,
                                   dimnames = list(NULL, "Z1")),
    matrix(x, ncol = 1L, dimnames = list(NULL, "Main_CS1")),
    penalty_V = c(Main_CS1 = 0.5)
  )
  fit0 <- SuSiEIRLS:::.cox_fit_fixed_ridge(
    pmin(te, tc), te <= tc, matrix(z, ncol = 1L,
                                   dimnames = list(NULL, "Z1")),
    matrix(x, ncol = 1L, dimnames = list(NULL, "Main_CS1"))
  )
  pen <- attr(fit, "refit_penalty")

  expect_s3_class(fit, "coxph.penal")
  expect_equal(pen$theta, c(Main_CS1 = 2))
  expect_false(pen$scale)
  expect_true(all(c("Z1", "Main_CS1") %in% names(stats::coef(fit))))
  expect_false("Z1" %in% names(pen$V))
  expect_lt(abs(stats::coef(fit)["Main_CS1"]),
            abs(stats::coef(fit0)["Main_CS1"]))
  expect_true(all(is.finite(stats::predict(fit, type = "lp"))))
})

test_that("multi-column fixed ridge preserves coefficient order", {
  set.seed(202607181)
  n <- 140L
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  y <- stats::rbinom(n, 1L, stats::plogis(x1 - 0.5 * x2))
  dat <- data.frame(y = y, Main_CS2 = x1, Main_CS4 = x2)
  fit <- SuSiEIRLS:::.mgcv_fit_fixed_ridge(
    "y", c("Main_CS2", "Main_CS4"), dat, stats::binomial(),
    penalty_V = c(Main_CS2 = 0.25, Main_CS4 = 0.5), mgcv_model = "gam"
  )

  expect_equal(
    attr(fit, "refit_penalty")$precision,
    c(Main_CS2 = 4, Main_CS4 = 2)
  )
  expect_true(all(c("Main_CS2", "Main_CS4") %in% names(stats::coef(fit))))
})

test_that("Gaussian paraPen uses previous dispersion for absolute 1/V", {
  set.seed(202607182)
  n <- 90L
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  phi <- 3.2
  V <- c(Main_CS1 = 0.7, Main_CS2 = 0.4)
  y <- as.numeric(X %*% c(0.8, -0.4) + stats::rnorm(n, sd = sqrt(phi)))
  dat <- data.frame(y = y, Main_CS1 = X[, 1L], Main_CS2 = X[, 2L])

  fit <- SuSiEIRLS:::.mgcv_fit_fixed_ridge(
    "y", names(V), dat, stats::gaussian(), V,
    dispersion = phi, mgcv_model = "gam"
  )
  fit_bam <- SuSiEIRLS:::.mgcv_fit_fixed_ridge(
    "y", names(V), dat, stats::gaussian(), V,
    dispersion = phi, mgcv_model = "bam"
  )
  A <- cbind(Intercept = 1, X)
  P <- diag(c(0, 1 / V), nrow = 3L)
  closed <- solve(crossprod(A) / phi + P, crossprod(A, y) / phi)
  pen <- attr(fit, "refit_penalty")

  expect_equal(unname(stats::coef(fit)[names(V)]), as.numeric(closed[-1L]),
               tolerance = 1e-10)
  expect_equal(
    unname(stats::coef(fit_bam)[names(V)]), as.numeric(closed[-1L]),
    tolerance = 1e-7
  )
  expect_equal(pen$precision, 1 / V)
  expect_equal(pen$working_precision, phi / V)
  expect_equal(pen$dispersion, phi)
})

test_that("refit variance mapping handles one column and invalid V", {
  fit <- list(V = c(0.25, 0))
  cs <- SuSiEIRLS:::.refit_penalty_variance(fit, 1L, "Main_CS1")
  expect_equal(as.numeric(cs), 0.25)
  expect_equal(names(cs), "Main_CS1")
  noncs <- SuSiEIRLS:::.refit_penalty_variance(
    fit, integer(0), "Main_noncs_res"
  )
  expect_equal(as.numeric(noncs), 0.25)
  expect_equal(names(noncs), "Main_noncs_res")
  fallback <- SuSiEIRLS:::.refit_penalty_variance(
    list(V = c(NA_real_, 0)), integer(0), "Main_noncs_res"
  )
  expect_equal(fallback, c(Main_noncs_res = 2))
  component_noncs <- SuSiEIRLS:::.refit_penalty_variance(
    list(V = c(0.2, 0.4)), integer(0), "Main_noncs_res"
  )
  expect_equal(component_noncs, c(Main_noncs_res = 0.4))
})

test_that("ordinal family routing is explicit", {
  route <- function(family, y) {
    X <- scale(matrix(stats::rnorm(length(y) * 3L), nrow = length(y)))
    SuSiE_IRLS(
      X = X, y = y, family = family, L = 1L, max.iter = 1L,
      min.iter = 1L, n_threads = 1L, scale_data = FALSE, verbose = FALSE,
      susie_para = list(max_iter = 10L)
    )
  }
  set.seed(20260719)
  y <- ordered(rep(1:3, each = 20L))

  fit_default <- route("ordinal", y)
  fit_logit <- route("clm_logit", y)
  fit_probit <- route("clm_probit", y)
  X <- scale(matrix(stats::rnorm(length(y) * 3L), nrow = length(y)))
  fit_ordered_default <- SuSiE_IRLS(
    X = X, y = y, L = 1L, max.iter = 1L, min.iter = 1L,
    n_threads = 1L, scale_data = FALSE, verbose = FALSE,
    susie_para = list(max_iter = 10L)
  )

  expect_s3_class(fit_default$fitJoint, "gam")
  expect_s3_class(fit_logit$fitJoint, "gam")
  expect_s3_class(fit_ordered_default$fitJoint, "gam")
  expect_s3_class(fit_probit$fitJoint, "clm")
  expect_error(route("clm", y), "incomplete")
})

test_that("ocat fixed ridge uses paraPen without penalizing Z", {
  set.seed(20260720)
  n <- 150L
  z <- stats::rnorm(n)
  x <- stats::rnorm(n)
  latent <- 0.2 * z + x + stats::rlogis(n)
  y <- findInterval(latent, c(-Inf, -0.5, 0.5, Inf))
  dat <- data.frame(y = y, Z1 = z, Main_CS1 = x)

  fit <- SuSiEIRLS:::.mgcv_fit_fixed_ridge(
    "y", c("Z1", "Main_CS1"), dat, mgcv::ocat(R = 3L),
    penalty_V = c(Main_CS1 = 0.5), mgcv_model = "gam"
  )
  pen <- attr(fit, "refit_penalty")

  expect_s3_class(fit, "gam")
  expect_match(fit$family$family, "Ordered Categorical")
  expect_equal(pen$precision, c(Main_CS1 = 2))
  expect_false("Z1" %in% names(pen$precision))
  expect_equal(unname(fit$full.sp), 1)
})

test_that("penalized coefficient names map back to main-index rows", {
  idx <- data.frame(CS = "Main_CS1")
  G <- matrix(
    c(0.2, 0.1, 0.1, 4, 1, 0.0455), nrow = 1L,
    dimnames = list("Main_CS1", c("coef", "se(coef)", "se2", "Chisq", "DF", "p"))
  )
  out <- SuSiEIRLS:::safe_add_p(idx, G)
  expect_equal(out$Pvalue, 0.0455)
})

test_that("ocat factor levels must match R", {
  y <- ordered(rep(letters[1:3], each = 3L), levels = letters[1:3])
  out <- SuSiEIRLS:::.ocat_prepare_response(y, mgcv::ocat(R = 3L))
  expect_equal(out$y_int, rep(1:3, each = 3L))
  expect_equal(out$ncat, 3L)
  expect_error(
    SuSiEIRLS:::.ocat_prepare_response(y, mgcv::ocat(R = 4L)),
    "does not match"
  )
})
