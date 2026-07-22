library(SuSiEIRLS)
library(mgcv)

set.seed(45)
n <- 160L
p <- 10L
X <- scale(matrix(rnorm(n * p), n, p))
Z <- scale(matrix(rnorm(n), n, 1))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- "Z1"
eta <- 0.2 + 1.0 * X[, 2] + 0.6 * Z[, 1]

y_nb <- rnbinom(n, mu = exp(eta), size = 2.75)
fit_nb_gam <- SuSiE_IRLS(X = X, Z = Z, y = y_nb, family = mgcv::nb(theta = NULL),
                  mgcv_model = "gam", scale_data = FALSE, n_threads = 1, L = 2,
                  min.iter = 1, max.iter = 2, susie_para = list(coverage = 0.5,
                      max_iter = 80, verbose = FALSE))
fit_nb_bam <- SuSiE_IRLS(X = X, Z = Z, y = y_nb, family = mgcv::nb(theta = 2.75),
                  mgcv_model = "bam", scale_data = FALSE, n_threads = 1, L = 2,
                  min.iter = 1, max.iter = 2, susie_para = list(coverage = 0.5,
                      max_iter = 80, verbose = FALSE))

y_tw <- mgcv::rTweedie(mu = exp(eta), p = 1.6, phi = 0.8)
fit_tw_gam <- SuSiE_IRLS(X = X, Z = Z, y = y_tw, family = mgcv::tw(theta = NULL,
                  link = "log"), mgcv_model = "gam", scale_data = FALSE, n_threads = 1,
                  L = 2, min.iter = 1, max.iter = 2, susie_para = list(coverage = 0.5,
                      max_iter = 80, verbose = FALSE))
fit_tw_bam <- SuSiE_IRLS(X = X, Z = Z, y = y_tw, family = mgcv::tw(theta = 1.6,
                  link = "log"), mgcv_model = "bam", scale_data = FALSE, n_threads = 1,
                  L = 2, min.iter = 1, max.iter = 2, susie_para = list(coverage = 0.5,
                      max_iter = 80, verbose = FALSE))

fam_zip <- mgcv::ziP(theta = c(-1, 0.2), b = 0.1)
y_zip <- fam_zip$rd(-0.4 + 0.8 * X[, 3] + 0.5 * Z[, 1], rep(1, n), 1)
fit_zip_gam <- SuSiE_IRLS(X = X, Z = Z, y = y_zip, family = mgcv::ziP(), mgcv_model = "gam",
                   scale_data = FALSE, n_threads = 1, L = 2, min.iter = 1, max.iter = 2,
                   susie_para = list(coverage = 0.6, max_iter = 80, verbose = FALSE))
fit_zip_bam <- SuSiE_IRLS(X = X, Z = Z, y = y_zip, family = mgcv::ziP(theta = c(-1,
                   0.2), b = 0.1), mgcv_model = "bam", scale_data = FALSE, n_threads = 1,
                   L = 2, min.iter = 1, max.iter = 2, susie_para = list(coverage = 0.6,
                       max_iter = 80, verbose = FALSE))

stopifnot(
  inherits(fit_nb_gam$fitJoint, "gam"),
  inherits(fit_nb_bam$fitJoint, "bam"),
  inherits(fit_tw_gam$fitJoint, "gam"),
  inherits(fit_tw_bam$fitJoint, "bam"),
  inherits(fit_zip_gam$fitJoint, "gam"),
  inherits(fit_zip_bam$fitJoint, "bam")
)
