# Adaptive Matern ("AMatern") basis for SuSiE_IRLS_GAM, following mgcv.taps
# (ported from the gam branch of SuSiE4I):
# the null space A = (1, z) is kept unpenalized and the Matern part is
# projected to be orthogonal to A (a projection, not a constraint). The basis is
# built once on a quantile grid of z and then frozen: every later gam()/bam()
# call only evaluates it, so repeated refits never redo the QR/eigen steps.

gam_matern_basis <- function(x, kappa, lambda) {
DX <- abs(outer(x, kappa, "-")) / lambda
DK <- abs(outer(kappa, kappa, "-")) / lambda
list(DX = exp(-DX) * (1 + DX), DK = exp(-DK) * (1 + DK))
}

gam_kappa_quantile <- function(x, nk) {
x <- sort(x)
n <- length(x)
stats::quantile(x[2:(n - 1)], seq(0, 1, length.out = nk + 2), names = FALSE)[2:(nk + 1)]
}

# Number of Matern knots: min(300, 0.3 * number of unique values).
gam_amatern_nk <- function(x) {
min(300L, as.integer(floor(0.3 * length(unique(x)))))
}

gam_amatern_base <- function(x, k = 10L, grid_size = 5000L) {
x <- as.numeric(x)
ux <- sort(unique(x))
nk <- gam_amatern_nk(x)
k <- min(as.integer(k), nk)
if (k < 4L) {
stop("A smooth needs at least 14 unique values; enter this variable as a linear term instead.",
     call. = FALSE)
}
grid <- if (length(ux) <= grid_size) ux else
  stats::quantile(x, seq(0, 1, length.out = grid_size), names = FALSE)
kappa <- gam_kappa_quantile(x, nk)
lambda <- 10 * stats::sd(x)
m <- 2L
A <- cbind(1, grid)
mb <- gam_matern_basis(grid, kappa, lambda)
C <- cbind(A, mb$DX)
G <- crossprod(C, A) / nrow(C)
Q <- qr.Q(qr(G), complete = TRUE)[, (m + 1):nrow(G), drop = FALSE]
B <- C %*% Q
v <- eigen(crossprod(B), symmetric = TRUE)$vectors[, seq_len(k - m + 1), drop = FALSE]
Omega <- matrix(0, ncol(C), ncol(C))
Omega[-(1:m), -(1:m)] <- mb$DK
Qv <- Q %*% v
Sr <- crossprod(Qv, Omega %*% Qv)
S <- matrix(0, m + ncol(v), m + ncol(v))
S[-(1:m), -(1:m)] <- (Sr + t(Sr)) / 2
list(kappa = kappa, lambda = lambda, Qv = Qv, S = S, m = m)
}

gam_amatern_predict <- function(base, x) {
x <- as.numeric(x)
A <- cbind(1, x)
cbind(A, cbind(A, gam_matern_basis(x, base$kappa, base$lambda)$DX) %*% base$Qv)
}

#' @importFrom mgcv smooth.construct Predict.matrix
#' @export
smooth.construct.sirlsAM.smooth.spec <- function(object, data, knots) {
base <- object$xt$base
X <- gam_amatern_predict(base, data[[object$term]])
S <- base$S
# A numeric by-variable multiplies the whole basis, so its constant column
# would duplicate the by-variable's own main effect: drop it, leaving z as the
# unpenalized null space of the varying coefficient.
drop_const <- !identical(object$by, "NA")
if (drop_const) {
X <- X[, -1L, drop = FALSE]
S <- S[-1L, -1L, drop = FALSE]
object$C <- matrix(0, 0L, ncol(X))
}
object$X <- X
object$S <- list(S)
object$base <- base
object$drop_const <- drop_const
object$null.space.dim <- base$m - drop_const
object$rank <- ncol(X) - object$null.space.dim
object$df <- ncol(X)
object$bs.dim <- ncol(X)
class(object) <- c("sirlsAM.smooth", "mgcv.smooth")
object
}

#' @export
Predict.matrix.sirlsAM.smooth <- function(object, data) {
X <- gam_amatern_predict(object$base, data[[object$term]])
if (isTRUE(object$drop_const)) X[, -1L, drop = FALSE] else X
}
