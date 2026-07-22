require(SuSiEIRLS)
require(ordinal)

set.seed(105)
n <- 300L
p <- 20L
q <- 3L
true_idx <- c(8L, 17L)

X <- matrix(stats::rnorm(n * p), n, p)
Z <- matrix(stats::rnorm(n * q), n, q)
X <- as.matrix(scale(X))
Z <- as.matrix(scale(Z))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- paste0("Z", seq_len(q))

beta <- rep(0, p)
beta[true_idx] <- c(0.45, -0.45)
alpha <- c(0.30, -0.20, 0.10)
eta <- as.numeric(X %*% beta + Z %*% alpha)
liab <- eta + stats::rnorm(n)
y <- cut(liab, breaks = c(-Inf, -0.4, 0.4, Inf), labels = FALSE)
y <- ordered(y, levels = 1:3)

fit <- SuSiEIRLS::SuSiE_IRLS(X = X, Z = Z, y = y, family = "clm_probit",
           L = 2L, L.init = 1L, max.iter = 3L, min.iter = 1L, max.eps = 0.001,
           n_threads = 1L, scale_data = FALSE, susie_para = list(max_iter = 80L,
               verbose = FALSE))

pip <- fit$fitX$pip[seq_len(p)]
if (length(pip) != p || any(!is.finite(pip))) stop("ordinal smoke test failed")
print(data.frame(family = "ordinal", iter = fit$diagnostics$iterations,
                 n_cs = nrow(fit$discovery_summary),
                 top_signal_pip = max(pip[true_idx])))
