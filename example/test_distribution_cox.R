require(SuSiEIRLS)
require(survival)

set.seed(102)
n <- 250L
p <- 20L
q <- 3L
true_idx <- c(5L, 13L)

X <- matrix(stats::rnorm(n * p), n, p)
Z <- matrix(stats::rnorm(n * q), n, q)
X <- as.matrix(scale(X))
Z <- as.matrix(scale(Z))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- paste0("Z", seq_len(q))

beta <- rep(0, p)
beta[true_idx] <- c(0.35, -0.35)
alpha <- c(0.25, -0.20, 0.10)
eta <- as.numeric(X %*% beta + Z %*% alpha)
Ttrue <- -log(stats::runif(n)) / (0.10 * exp(eta))
Ctime <- stats::rexp(n, rate = 0.08)
y <- survival::Surv(pmin(Ttrue, Ctime), as.integer(Ttrue <= Ctime))

fit <- SuSiEIRLS::SuSiE_IRLS(X = X, Z = Z, y = y, L = 2L, L.init = 1L,
           max.iter = 3L, min.iter = 1L, max.eps = 0.001, n_threads = 1L,
           scale_data = FALSE, susie_para = list(max_iter = 80L, verbose = FALSE))

pip <- fit$fitX$pip[seq_len(p)]
if (length(pip) != p || any(!is.finite(pip))) stop("cox smoke test failed")
print(data.frame(family = "cox", iter = fit$diagnostics$iterations,
                 n_cs = nrow(fit$discovery_summary),
                 top_signal_pip = max(pip[true_idx])))
