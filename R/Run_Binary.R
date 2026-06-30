#' Binary logistic IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @export
Run_Binary <- function(X, y, Z = NULL,
                    family = binomial(link = "logit"),weight_cutoff=0.005,
                    L, max.iter, min.iter, max.eps, susie.iter,
                    verbose = TRUE, n_threads = 1, coverage = 0.9,
                    estimate_residual_variance = TRUE,scaled_prior_variance=1,
                    residual_variance = 0.5,
                    residual_variance_lowerbound = 0.5,
                    residual_variance_upperbound = 1,
                    L.init = 1,
                    init_cor_method = NULL,
                    refit_noncs = TRUE,
                    noncs_var = 0.2,
                    suff_block_size = 10000L, ...) {

n = n_eff= length(y)
p = ncol(X)
suff_block_size <- validate_suff_block_size(suff_block_size)

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
init_cor_method = init_cor_method
)
alpha = clean_coef(coef(fit_final)[seq_len(ncol(ZI))])

# Initialize tracking variables
g = c()
beta = rep(0, p)
beta_prev = beta
alpha_prev = alpha * 0
early_no_cs <- FALSE
XCS <- NULL

# ============================================
# Main iteration loop
# ============================================
for (iter in 1:max.iter) {
beta_prev = beta
alpha_prev = alpha

## ===== Logistic (Bernoulli, n_i = 1) : PG  =====

eta <- fit_final$linear.predictors
y01 <- as.numeric(fit_final$y)
eta_clip <- pmin(pmax(eta, -20), 20)
omega <- ifelse(abs(eta_clip) < 1e-8, 0.25, 0.5 * tanh(eta_clip / 2) / eta_clip)
z      <- (y01 - 0.5) / omega
W_diag <- omega
W_diag <- pmax(W_diag, 1e-8)
n_eff <- (sum(W_diag))^2 / sum(W_diag^2)

# Compute projected sufficient statistics without materializing tilde_X.
suff = weighted_residual_suffstats(
X = X,
y = z,
ZI = ZI,
weights = W_diag,
n_threads = n_threads,
block_size = suff_block_size
)
XtX = suff$XtX
Xty = suff$Xty
yty = suff$yty
rm(suff)

# Run SuSiE on projected data
fitX <- susie_ss(
XtX = XtX, Xty = Xty, yty = yty, n = max(n/2,n_eff), L = L,
scaled_prior_variance = scaled_prior_variance,
estimate_residual_variance = FALSE,
residual_variance = 1,
max_iter = susie.iter,
estimate_prior_method = "EM",
coverage = coverage,...
)

beta = clean_coef(coef(fitX)[-1])

# Extract credible sets using summary information
CSdt <- summary(fitX)$vars
cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
cs_indices=sort(cs_indices)
if(length(cs_indices) == 0) {
if (iter <= min.iter) {
noncs_res <- build_no_cs_noncs_refit_term(X, fitX)
if (is.null(noncs_res)) {
early_no_cs <- TRUE
if (verbose) {
cat("No credible set detected; returning current no-CS fit.\n")
}
break
}
XCS <- matrix(noncs_res, ncol = 1)
colnames(XCS) <- "Main_CS_noncs"
XCS <- as.matrix(XCS)
XCS_refit <- XCS
} else {
early_no_cs <- TRUE
if (verbose) {
cat("No credible set detected; returning current no-CS fit.\n")
}
break
}
} else {
Alpha_filtered <- fitX$alpha * 0
for(i in cs_indices) {
vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i]
}
# Align within-CS SNP directions while preserving PIP weights.
Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
XCS <- matrixMultiply(X, as.matrix(Alpha_filtered), transB = TRUE)
XCS <- XCS[, cs_indices, drop = FALSE]
if(is.null(dim(XCS))) {
XCS <- matrix(XCS, ncol = 1)
}
colnames(XCS) <- paste0("Main_CS", cs_indices)
XCS <- as.matrix(XCS)
XCS_refit <- XCS
if (isTRUE(refit_noncs)) {
noncs_term <- build_noncs_refit_term(
X = X, fitX = fitX, CSdt = CSdt, cs_indices = cs_indices,
XCS = XCS, noncs_var = noncs_var
)
if (!is.null(noncs_term)) {
XCS_refit <- cbind(XCS_refit, Main_CS_noncs = noncs_term)
}
}
}

# ============================================
# Refit GLM with selected credible sets
# ============================================
if (ncol(Z) == 0) {
# No covariates: only intercept and XCS
Data = data.frame(y = y, XCS_refit)
} else {
# With covariates
Data = cbind(y, Z, XCS_refit)
Data = as.data.frame(Data)
}

fit_final = glm(y ~ ., data = Data, family = family)

# Extract covariate coefficients (intercept + Z)
alpha = clean_coef(coef(fit_final)[1:(ncol(ZI))])

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
if (early_no_cs) {
MainIndex <- Identifying_MainEffect(fitX, colnames(X))
G <- summary(fit_final)$coefficients
MainIndex <- safe_add_p(MainIndex, G)
fit_final <- clean_model_environment(fit_final)
return(list(
iter = iter,
error = g,
converged = FALSE,
fitX = fitX,
fitJoint = fit_final,
main_index = MainIndex,
JointCoef = G
))
}

pg_fitX <- fitX
pg_XCS <- XCS
pg_converged <- (iter < max.iter && err < max.eps)
eta <- fit_final$linear.predictors
mu <- fit_final$fitted.values
g_prime_mu <- 1 / fit_final$family$mu.eta(eta)
var_mu <- fit_final$family$variance(mu)
pseudo_response <- eta + (y - mu) * g_prime_mu
W_diag <- 1 / (var_mu * g_prime_mu^2)
bad <- !is.finite(pseudo_response) | !is.finite(W_diag) | (W_diag <= 0)
if (mean(bad) > 0.9) {
stop("Too many invalid observations in final logit correction")
}
if (any(bad)) {
W_diag[bad] <- 0
pseudo_response[bad] <- 0
}
W_diag <- robust_weight(W_diag, cutoff = weight_cutoff)
weight_denom <- sum(W_diag^2)
if (!is.finite(weight_denom) || weight_denom <= 0) {
stop("All working weights are zero in final logit correction")
}
n_eff <- (sum(W_diag))^2 / weight_denom
phi0 <- summary(fit_final)$dispersion
suff <- weighted_residual_suffstats(
X = X,
y = pseudo_response,
ZI = ZI,
weights = W_diag / phi0,
n_threads = n_threads,
block_size = suff_block_size
)
fitX <- susie_ss(
XtX = suff$XtX, Xty = suff$Xty, yty = suff$yty,
n = max(n/2, n_eff), L = L,
scaled_prior_variance = scaled_prior_variance,
estimate_residual_variance = estimate_residual_variance,
residual_variance = residual_variance,
residual_variance_lowerbound = residual_variance_lowerbound,
residual_variance_upperbound = residual_variance_upperbound,
max_iter = susie.iter,
estimate_prior_method = "EM",
coverage = coverage,...
)
rm(suff)
CSdt <- summary(fitX)$vars
cs_indices <- sort(unique(CSdt$cs[CSdt$cs > 0]))
if (length(cs_indices) == 0) {
if (ncol(Z) == 0) {
Data = data.frame(y = y, pg_XCS)
} else {
Data = cbind(y, Z, pg_XCS)
Data = as.data.frame(Data)
}
fit_final = glm(y ~ ., data = Data, family = family)
MainIndex = Identifying_MainEffect(pg_fitX, colnames(X))
G = summary(fit_final)$coefficients
MainIndex <- safe_add_p(MainIndex, G)
if (!is.null(MainIndex)) MainIndex$status <- "PG IRLS"
fit_final <- clean_model_environment(fit_final)
return(list(
iter = iter,
error = g,
converged = pg_converged,
fitX = pg_fitX,
fitJoint = fit_final,
main_index = MainIndex,
JointCoef = G
))
}
Alpha_filtered <- fitX$alpha * 0
for (i in cs_indices) {
vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i]
}
Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
XCS <- matrixMultiply(X, as.matrix(Alpha_filtered), transB = TRUE)
XCS <- XCS[, cs_indices, drop = FALSE]
if (is.null(dim(XCS))) {
XCS <- matrix(XCS, ncol = 1)
}
colnames(XCS) <- paste0("Main_CS", cs_indices)
XCS <- as.matrix(XCS)

if (ncol(Z) == 0) {
Data = data.frame(y = y, XCS)
} else {
Data = cbind(y, Z, XCS)
Data = as.data.frame(Data)
}
fit_final = glm(y ~ ., data = Data, family = family)
MainIndex = Identifying_MainEffect(fitX, colnames(X))
G = summary(fit_final)$coefficients
MainIndex <- safe_add_p(MainIndex, G)
fit_final <- clean_model_environment(fit_final)

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
