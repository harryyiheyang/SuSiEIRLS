require(SuSiEIRLS)
require(mgcv)

set.seed(104)
n <- 300L
p <- 20L
q <- 3L
true_idx <- c(7L, 16L)

X <- matrix(stats::rnorm(n * p), n, p)
Z <- matrix(stats::rnorm(n * q), n, q)
X <- as.matrix(scale(X))
Z <- as.matrix(scale(Z))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- paste0("Z", seq_len(q))

beta <- rep(0, p)
beta[true_idx] <- c(0.25, -0.25)
alpha <- c(0.20, -0.15, 0.10)
eta <- as.numeric(X %*% beta + Z %*% alpha)
a <- stats::uniroot(
  function(a0) mean(exp(a0 + eta)) - 2.5,
  interval = c(-20, 20)
)$root
mu <- exp(a + eta)
y <- mgcv::rTweedie(mu = mu, p = 1.3, phi = 1.2)

fit <- SuSiEIRLS::SuSiE_IRLS(X = X, Z = Z, y = y, family = mgcv::tw(theta = NULL,
           link = "log"), L = 2L, L.init = 1L, max.iter = 3L, min.iter = 1L,
           max.eps = 0.001, n_threads = 1L, scale_data = FALSE, susie_para = list(max_iter = 80L,
               verbose = FALSE))

pip <- fit$fitX$pip[seq_len(p)]
if (length(pip) != p || any(!is.finite(pip))) stop("Tweedie smoke test failed")
print(data.frame(family = "tw", iter = fit$diagnostics$iterations,
                 n_cs = nrow(fit$discovery_summary),
                 top_signal_pip = max(pip[true_idx]),
                 p_hat = fit$theta[1L]))
