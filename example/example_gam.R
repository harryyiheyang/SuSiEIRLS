library(SuSiEIRLS)

set.seed(3)
n <- 3000L
p <- 200L
age <- runif(n, 40, 75)
sex <- factor(sample(c("F", "M"), n, replace = TRUE))
X <- scale(matrix(rnorm(n * p), n, p))
colnames(X) <- paste0("rs", seq_len(p))

u <- ecdf(age)(age)
male <- as.numeric(sex == "M")
f_age <- scale(plogis(6 * (u - 0.4)) * (1 + 0.5 * male))[, 1] * sqrt(0.3)
eta <- f_age + 0.3 * male + 0.15 * (X[, 5] + X[, 105] + X[, 155])
dat <- data.frame(age = age, sex = sex)

# Gaussian: a sex-specific age curve as the null model.
# Any bs is replaced by AMatern (with a warning).
dat$y <- eta + rnorm(n, sd = sqrt(0.6))
fit <- SuSiE_IRLS_GAM(y ~ s(age, by = sex) + sex, data = dat, X = X,
                      family = gaussian(), mgcv_model = "bam")
fit$discovery_summary

# Binary outcome
dat$yb <- rbinom(n, 1, plogis(-0.5 + 2 * eta))
fit_b <- SuSiE_IRLS_GAM(yb ~ s(age, by = sex) + sex, data = dat, X = X,
                        family = binomial(), verbose = FALSE)
fit_b$discovery_summary

# Poisson with gam (REML)
dat$yp <- rpois(n, exp(0.2 + 0.5 * eta))
fit_p <- SuSiE_IRLS_GAM(yp ~ s(age) + sex, data = dat, X = X,
                        family = poisson(), mgcv_model = "gam", verbose = FALSE)
fit_p$discovery_summary
