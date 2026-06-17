get_consecutive <- function(max, number) {
start <- sample(1:(max - number + 1), 1)
return(start:(start + number - 1))
}
###############################################################################
demean=function(X){
Xname=colnames(X)
p=ncol(X)
if(is.null(p)){
X=X-mean(X)
}else{
Xmean=colMeans(X)
X=sweep(X,2,Xmean,"-")
}
colnames(X)=Xname
return(X)
}
###############################################################################
get_pairwise_interactions <- function(W, Z=NULL) {
n <- nrow(W)
p <- ncol(W)
q <- ncol(Z)
#
colnames_W <- colnames(W)
colnames_Z <- colnames(Z)
# --- W by W interactions ---
ww_cols <- p * (p + 1) / 2
WW <- matrix(NA, n, ww_cols)
colnames_WW <- character(ww_cols)
#
col_idx <- 1
for (i in 1:p) {
for (j in i:p) {
WW[, col_idx] <- W[, i] * W[, j]
colnames_WW[col_idx] <- paste0(colnames_W[i], "*", colnames_W[j])
col_idx <- col_idx + 1
}
}
# --- Z by W interactions ---
if(is.null(q)==F){
zw_cols <- q * p
ZW <- matrix(NA, n, zw_cols)
colnames_ZW <- character(zw_cols)
#
col_idx <- 1
for (i in 1:q) {
for (j in 1:p) {
ZW[, col_idx] <- Z[, i] * W[, j]
colnames_ZW[col_idx] <- paste0(colnames_Z[i], "*", colnames_W[j])
col_idx <- col_idx + 1
}
}
#
out <- cbind(WW, ZW)
colnames(out) <- c(colnames_WW, colnames_ZW)
}else{
out = WW
colnames(out)=colnames_WW
}
return(out)
}
###############################################################################
get_active_indices <- function(fit) {
cs = tryCatch(summary(fit), error = function(e) NULL)
if (!is.null(cs) && length(cs$cs) > 0) {
active_idx = unique(unlist(cs$cs$cs))
}else{
active_idx=NULL
}
return(active_idx)
}
###############################################################################
Identifying_MainEffect=function(fit,nam){
summ=summary(fit)$vars
g=unique(summ$cs[which(summ$cs>0)])
if(length(g)==0){
stop("No credible set detected")
}
bb=summary(fit)$cs
S=list()
for(i in g){
indi=which(summ$cs==i)
a=summ$variable[indi]
b=data.frame(Index=a,Variable=nam[summ$variable[indi]],CS=paste0("Main_CS",i),log10BF=bb$cs_log10bf[bb$cs==i],PIP=summ$variable_prob[indi])
S[[i]]=b
}
return(do.call(rbind,S))
}
###############################################################################
Identifying_EnvEffect=function(fit,nam){
  summ=summary(fit)$vars
  g=unique(summ$cs[which(summ$cs>0)])
  bb=summary(fit)$cs
  S=list()
  for(i in g){
    indi=which(summ$cs==i)
    a=summ$variable[indi]
    b=data.frame(Index=a,Variable=nam[summ$variable[indi]],CS=paste0("Env_CS",i),log10BF=bb$cs_log10bf[bb$cs==i],PIP=summ$variable_prob[indi])
    S[[i]]=b
  }
  return(do.call(rbind,S))
}
###############################################################################
Identifying_IntEffect=function(fitW,namW){
summ=summary(fitW)$vars
if(length(which(summ$cs>0))>0){
bb=summary(fitW)$cs
g=unique(summ$cs[which(summ$cs>0)])
S=list()
for(i in g){
indi=which(summ$cs==i)
a=summ$variable[indi]
b=data.frame(Index=a,Variable=namW[summ$variable[indi]],CS=paste0("Int_CS",i),log10BF=bb$cs_log10bf[bb$cs==i],PIP=summ$variable_prob[indi])
S[[i]]=b
}
return(do.call(rbind,S))
}else{
return(NULL)
}
}

ProjectRes=function(A,B,inercept=F,n_threads){
if(inercept==T){
B=cbind(1,B)
}
BtB = blockwise_crossprod(X=B,n_threads=n_threads)
BtA = blockwise_crossprod(B,A,n_threads)
ProjPart = matrixMultiply(B,(solve(BtB)%*%(BtA)))
return(ProjPart)
}

safe_add_p <- function(idx, Coefmat) {
  if (is.null(idx)) return(NULL)
  if (is.data.frame(idx) && nrow(idx) == 0) return(idx)
  if (!("CS" %in% names(idx))) return(idx)
  cs <- as.character(idx$CS)
  pos <- match(cs, rownames(Coefmat))
  p   <- rep(NA_real_, length(cs))
  p[!is.na(pos)] <- Coefmat[pos[!is.na(pos)], 4]

  idx$Pvalue <- p
  idx
}

################################################################################
robust_weight <- function(w, cutoff = 0.01) {
n <- length(w)
if (n == 0L || all(is.na(w))) return(w)
w <- as.numeric(w)
w[!is.finite(w)] <- NA
w_na <- is.na(w)

n_eff <- sum(!w_na)
if (n_eff <= 1L) return(w)

if (n_eff < 1 / cutoff) {
lo <- min(w, na.rm = TRUE)
hi <- max(w, na.rm = TRUE)
} else {
lo <- quantile(w, probs = cutoff, na.rm = TRUE, names = FALSE, type = 7)
hi <- quantile(w, probs = 1 - cutoff, na.rm = TRUE, names = FALSE, type = 7)
}
w_trim <- pmin(pmax(w, lo), hi)
w[w_na] <- NA
w_trim
}

init_k_from_L <- function(L.init, p) {
  if (!is.numeric(L.init) || length(L.init) != 1L ||
      !is.finite(L.init) || L.init < 1) {
    stop("L.init must be a positive numeric scalar.")
  }
  min(as.integer(ceiling(L.init)), as.integer(p))
}

make_init_data <- function(y = NULL, Z = NULL, X = NULL, selected = integer(0)) {
  if (!is.null(y)) {
    n <- length(y)
    Data <- data.frame(y = y)
  } else if (!is.null(Z) && ncol(Z) > 0) {
    n <- nrow(Z)
    Data <- data.frame(row.names = seq_len(n))
  } else {
    n <- nrow(X)
    Data <- data.frame(row.names = seq_len(n))
  }

  if (!is.null(Z) && ncol(Z) > 0) {
    Zdf <- as.data.frame(Z)
    colnames(Zdf) <- paste0("Z", seq_len(ncol(Z)))
    Data <- cbind(Data, Zdf)
  }

  if (length(selected) > 0) {
    Xdf <- as.data.frame(X[, selected, drop = FALSE])
    colnames(Xdf) <- paste0("InitX", seq_along(selected))
    Data <- cbind(Data, Xdf)
  }

  Data
}

fit_init_glm <- function(X, y, Z, selected, family) {
  Data <- make_init_data(y = y, Z = Z, X = X, selected = selected)
  rhs_n <- ncol(Data) - 1L
  if (rhs_n == 0L) {
    stats::glm(y ~ 1, data = Data, family = family)
  } else {
    stats::glm(y ~ ., data = Data, family = family)
  }
}

fit_init_nb <- function(X, y, Z, selected, theta_init, estimate_theta) {
  Data <- make_init_data(y = y, Z = Z, X = X, selected = selected)
  rhs_n <- ncol(Data) - 1L
  fml <- if (rhs_n == 0L) y ~ 1 else y ~ .

  if (estimate_theta) {
    tryCatch(
      MASS::glm.nb(fml, data = Data, link = "log"),
      error = function(e) {
        stats::glm(
          fml, data = Data,
          family = MASS::negative.binomial(theta_init, link = "log")
        )
      }
    )
  } else {
    stats::glm(
      fml, data = Data,
      family = MASS::negative.binomial(theta_init, link = "log")
    )
  }
}

fit_init_cox <- function(X, y, status, Z, selected) {
  surv_y <- survival::Surv(y, status)
  Data <- make_init_data(Z = Z, X = X, selected = selected)
  if (ncol(Data) == 0L) {
    survival::coxph(surv_y ~ 1, ties = "breslow", model = TRUE)
  } else {
    survival::coxph(surv_y ~ ., data = Data, ties = "breslow", model = TRUE)
  }
}

select_by_residual_cor <- function(X, residual, available,
                                   cor_method = c("pearson", "spearman")) {
  cor_method <- match.arg(cor_method)
  if (!any(available)) return(NA_integer_)

  r <- as.numeric(residual)
  ok <- is.finite(r)
  if (sum(ok) < 3L) return(NA_integer_)
  if (stats::sd(r[ok]) == 0) return(NA_integer_)

  Xsub <- X[ok, available, drop = FALSE]
  scores <- suppressWarnings(
    as.numeric(stats::cor(
      Xsub, r[ok],
      use = "pairwise.complete.obs",
      method = cor_method
    ))
  )
  scores[!is.finite(scores)] <- NA_real_
  if (all(is.na(scores))) return(NA_integer_)

  which_available <- which(available)
  which_available[which.max(abs(scores))]
}

greedy_glm_warm_start <- function(X, y, Z, family, L.init = 1,
                                  cor_method = c("pearson", "spearman")) {
  cor_method <- match.arg(cor_method)
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- fit_init_glm(X = X, y = y, Z = Z, selected = selected, family = family)

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "response")
    j <- select_by_residual_cor(
      X = X, residual = r, available = available, cor_method = cor_method
    )
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    fit <- fit_init_glm(X = X, y = y, Z = Z, selected = selected, family = family)
  }

  fit
}

greedy_nb_warm_start <- function(X, y, Z, L.init = 1, theta_init, estimate_theta,
                                 cor_method = c("pearson", "spearman")) {
  cor_method <- match.arg(cor_method)
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- fit_init_nb(
    X = X, y = y, Z = Z, selected = selected,
    theta_init = theta_init, estimate_theta = estimate_theta
  )

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "response")
    j <- select_by_residual_cor(
      X = X, residual = r, available = available, cor_method = cor_method
    )
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    fit <- fit_init_nb(
      X = X, y = y, Z = Z, selected = selected,
      theta_init = theta_init, estimate_theta = estimate_theta
    )
  }

  fit
}

greedy_cox_warm_start <- function(X, y, status, Z, L.init = 1,
                                  cor_method = c("pearson", "spearman")) {
  cor_method <- match.arg(cor_method)
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- fit_init_cox(X = X, y = y, status = status, Z = Z, selected = selected)

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "martingale")
    j <- select_by_residual_cor(
      X = X, residual = r, available = available, cor_method = cor_method
    )
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    fit <- fit_init_cox(X = X, y = y, status = status, Z = Z, selected = selected)
  }

  fit
}

is_logit_binomial <- function(fam) {
  inherits(fam, "family") &&
    fam$link == "logit" &&
    fam$family %in% c("binomial", "quasibinomial")
}

group.pip.filter=function(pip.summary,xQTL.cred.thres=0.95,xQTL.pip.thres=0.1){
  ind=which(pip.summary$cs>0)
  if(length(ind)>0){
    J=max(pip.summary$cs[ind])
    pip.summary$cs.pip=pip.summary$variable_prob
    for(i in 1:J){
      indi=which(pip.summary$cs==i)
      summaryi=pip.summary[indi,]
      pip.cred=sum(summaryi$variable_prob)
      pip.summary$cs.pip[indi]=pip.cred
    }
    ind.keep=which(pip.summary$cs.pip>=xQTL.cred.thres&pip.summary$variable_prob>=xQTL.pip.thres)
    cs=pip.summary$cs
    cs.pip=pip.summary$cs.pip
    cs->cs[pip.summary$variable]
    cs.pip->cs.pip[pip.summary$variable]
    cs[which(cs==-1)]=0
  }else{
    ind.keep=NULL
    cs=pip.summary$cs.pip*0
    cs.pip=pip.summary$cs.pip*0
  }
  return(list(ind.keep=pip.summary$variable[ind.keep],cs=cs,cs.pip=cs.pip,result=pip.summary))
}


get_category_probs <- function(fit_clm, link = NULL) {
  n <- fit_clm$n
  J <- length(fit_clm$y.levels)
  K <- J - 1  # Number of thresholds = J - 1

  if(J < 2) {
    stop("Need at least 2 categories for ordinal model")
  }

  # Determine link function
  if(is.null(link)) {
    # Try to extract from fit_clm object
    if(!is.null(fit_clm$link)) {
      link <- fit_clm$link
    } else {
      link <- "logit"  # default
      warning("Link function not specified, using 'logit' as default")
    }
  }

  # Validate link function
  link <- match.arg(link, c("logit", "probit", "cloglog"))

  # Select appropriate CDF function
  F_link <- switch(link,
                   "logit" = plogis,
                   "probit" = pnorm,
                   "cloglog" = function(x) 1 - exp(-exp(x))
  )

  alpha <- as.numeric(fit_clm$alpha)
  if(length(alpha) != K) {
    stop(sprintf("Expected %d thresholds but got %d", K, length(alpha)))
  }

  if(length(fit_clm$beta) > 0) {
    formula_terms <- delete.response(terms(fit_clm))
    X_fit <- model.matrix(formula_terms, data = fit_clm$model)
    if("(Intercept)" %in% colnames(X_fit)) {
      X_fit <- X_fit[, -1, drop = FALSE]
    }
    if(ncol(X_fit) != length(fit_clm$beta)) {
      stop(sprintf("Dimension mismatch: X has %d columns but beta has %d elements",
                   ncol(X_fit), length(fit_clm$beta)))
    }
    h <- matrixVectorMultiply(X_fit, fit_clm$beta)
  } else {
    h <- rep(0, n)
  }

  # Compute cumulative probabilities using the appropriate link function
  alpha_full <- c(-Inf, alpha, Inf)
  eta_left  <- outer(alpha_full[1:J], h, function(a, hi) a - hi)     # J by n
  eta_right <- outer(alpha_full[2:(J+1)], h, function(a, hi) a - hi) # J by n

  # Apply the selected link function
  m_mat_full <- t(F_link(eta_right) - F_link(eta_left))  # n by J
  m_mat <- m_mat_full[, 1:K, drop = FALSE]  # n by (J-1)

  if(any(m_mat < -1e-10 | m_mat > 1 + 1e-10)) {
    warning("Some probabilities outside [0, 1]")
    m_mat <- pmax(0, pmin(1, m_mat))
  }

  row_sums_partial <- rowSums(m_mat)
  if(any(row_sums_partial > 1 + 1e-6)) {
    warning("Sum of first J-1 probabilities exceeds 1")
  }

  colnames(m_mat) <- paste0("mu_", 1:K)

  # Add link function as attribute for reference
  attr(m_mat, "link") <- link

  return(m_mat)
}

estimate_sigma2_null <- function(tilde_y, X, W_diag, strata,
                                 ZI = NULL, m = 200,
                                 method = c("median","mean")) {
  method <- match.arg(method)
  n <- length(tilde_y)
  p <- ncol(X)
  Xnull <- X
  levs <- sort(unique(strata))
  for (s in levs) {
    idx <- which(strata == s)
    Xnull[idx, ] <- Xnull[sample(idx, length(idx), replace = FALSE), , drop = FALSE]
  }
  m <- min(m, p)
  cols <- sample.int(p, m, replace = FALSE)
  G0 <- Xnull[, cols, drop = FALSE]
  W_sqrt <- sqrt(W_diag)
  tilde_G0 <- G0 * W_sqrt
  if (!is.null(ZI) && ncol(ZI) > 0) {
    tilde_Z <- ZI * W_sqrt
    tilde_G0 <- tilde_G0 - ProjectRes(A = tilde_G0, B = tilde_Z, n_threads = 1)
  }
  U <- as.numeric(crossprod(tilde_G0, tilde_y))
  V <- as.numeric(colSums(tilde_G0^2))
  good <- is.finite(U) & is.finite(V) & (V > 0)

  Z2 <- (U[good]^2) / V[good]
  if (!length(Z2)) return(1)

  if (method == "mean") {
    sigma2_hat <- mean(Z2)
  } else {
    sigma2_hat <- median(Z2) / 0.4559364
  }

  max(1, sigma2_hat)
}

