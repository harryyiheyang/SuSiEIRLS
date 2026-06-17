#' General GLM IRLS-SuSiE path
#' @export
Run_GLM <- function(X, y, Z = NULL,weight_cutoff=0.005,
                    family = binomial(link = "logit"),
                    L, max.iter, min.iter, max.eps, susie.iter,
                    verbose = TRUE, n_threads = 1, coverage = 0.95,
                    estimate_residual_variance = FALSE,
                    residual_variance = 1,scaled_prior_variance=1,
                    pip.thres=0.005,
                    L.init = 1,
                    init_cor_method = c("pearson", "spearman"), ...) {

n = n_eff = length(y)
p = ncol(X)
init_cor_method <- match.arg(init_cor_method)

# ============================================
# Handle Z edge cases
# ============================================
if (is.null(Z)) {
# Case: Z is NULL (no covariates)
Z = matrix(nrow = n, ncol = 0)  # Create empty matrix with 0 columns
ZI = matrix(1, nrow = n, ncol = 1)
colnames(ZI) = "Intercept"

} else {
# Case: Z is not NULL
if (is.null(dim(Z))) {
  # Z is a vector: convert to single-column matrix
  Z = matrix(Z, ncol = 1)
}

# Ensure Z has column names
if (is.null(colnames(Z))) {
  colnames(Z) = paste0("Z", seq_len(ncol(Z)))
}

# Create ZI = [Intercept | Z]
ZI = cbind(1, Z)
colnames(ZI)[1] = "Intercept"
}

# ============================================
# Greedy low-dimensional GLM warm start
# ============================================
fit_final = greedy_glm_warm_start(
  X = X, y = y, Z = Z, family = family, L.init = L.init,
  cor_method = init_cor_method
)

alpha = coef(fit_final)[seq_len(ncol(ZI))]

# Initialize tracking variables
g = c()
beta = rep(0, p)
beta_prev = beta
alpha_prev = alpha * 0

# ============================================
# Main iteration loop
# ============================================
for (iter in 1:max.iter) {
beta_prev = beta
alpha_prev = alpha

# Extract information from current GLM fit
eta <- fit_final$linear.predictors
mu <- fit_final$fitted.values
mu_eta_func <- fit_final$family$mu.eta
g_prime_mu <- 1 / mu_eta_func(eta)
family_info <- fit_final$family
var_mu <- family_info$variance(mu)

# Construct pseudo-response and weights (IRWLS)
pseudo_response <- eta + (y - mu) * g_prime_mu
W_diag <- 1 / (var_mu * g_prime_mu^2)
W_diag = robust_weight(W_diag,cutoff=weight_cutoff)
n_eff <-  (sum(W_diag))^2 / sum(W_diag^2)
phi0 <- summary(fit_final)$dispersion
W_invsqrt = sqrt(W_diag / phi0)

# Weighted transformation
tilde_y = pseudo_response * W_invsqrt
tilde_X = X * W_invsqrt
tilde_Z = ZI * W_invsqrt

# Projection: remove covariate effects
ZtZ = matrixMultiply(t(tilde_Z), tilde_Z)
Zinv = solve(ZtZ)
Zinv = matrixMultiply(Zinv, t(tilde_Z))

tilde_y = tilde_y - matrixVectorMultiply(tilde_Z, matrixVectorMultiply(Zinv, tilde_y))
tilde_X = tilde_X - ProjectRes(A = tilde_X, B = tilde_Z, n_threads = n_threads)

# Compute sufficient statistics for SuSiE
XtX = blockwise_crossprod(tilde_X,n_threads = n_threads)
Xty = matrixVectorMultiply(t(tilde_X), tilde_y)
yty = sum(tilde_y^2)

# Run SuSiE on projected data
fitX <- susie_ss(
  XtX = XtX, Xty = Xty, yty = yty, n = max(n/2,n_eff), L = L,
  scaled_prior_variance = scaled_prior_variance,
  estimate_residual_variance = estimate_residual_variance,
  residual_variance = residual_variance,
  max_iter = susie.iter,
  estimate_prior_method = "EM",
  coverage=coverage,...
)

beta = coef(fitX)[-1]
beta.cs=group.pip.filter(pip.summary=summary(fitX)$var,xQTL.cred.thres=coverage,xQTL.pip.thres=pip.thres)
pip.alive=beta.cs$ind.keep
beta[-pip.alive]=0

# Extract credible sets using summary information
CSdt <- summary(fitX)$vars
cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
cs_indices=sort(cs_indices)
if(length(cs_indices) == 0) {
stop("No credible set detected at iteration ", iter)
}
Alpha_filtered <- fitX$alpha * 0
for(i in cs_indices) {
vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i]
}
# Align within-CS SNP directions while preserving PIP weights.
Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
XCS <- matrixMultiply(X, t(as.matrix(Alpha_filtered)))
XCS <- XCS[, cs_indices, drop = FALSE]
if(is.null(dim(XCS))) {
XCS <- matrix(XCS, ncol = 1)
}
colnames(XCS) <- paste0("Main_CS", cs_indices)
XCS <- as.matrix(XCS)

# ============================================
# Refit GLM with selected credible sets
# ============================================
if (ncol(Z) == 0) {
  # No covariates: only intercept and XCS
  Data = data.frame(y = y, XCS)
} else {
  # With covariates
  Data = cbind(y, Z, XCS)
  Data = as.data.frame(Data)
}

fit_final = glm(y ~ ., data = Data, family = family)

# Extract covariate coefficients (intercept + Z)
alpha = coef(fit_final)[1:(ncol(ZI))]

# Check convergence
err = max(sqrt(mean((beta - beta_prev)^2)),
          sqrt(mean((alpha - alpha_prev)^2)))
g[iter] = err

if (verbose) {
  cat(sprintf("Iteration %d: err = %.3e\n", iter, err))
}

if (err < max.eps && iter > min.iter) {
  if (verbose) cat("Converged!\n")
  break
}
}

# ============================================
# Post-processing
# ============================================
MainIndex = Identifying_MainEffect(fitX, colnames(X))
G = summary(fit_final)$coefficients
MainIndex <- safe_add_p(MainIndex, G)

if (verbose) {
plot(g, type = "o", col = "black", pch = 16,
     xlab = "Iteration",
     ylab = "Max Parameter Change",
     main = "Convergence Trace (Max |Delta| in alpha and beta)")
for (i in seq_along(g)) {
  text(x = i, y = g[i],
       labels = formatC(g[i], format = "e", digits = 1),
       pos = 3, cex = 0.7, col = "red")
}
}

AA = list(
iter = iter,
error = g,
converged = (iter < max.iter && err < max.eps),
fitX = fitX,
fitJoint = fit_final,
main_index = MainIndex,
JointCoef = G
)

return(AA)
}
