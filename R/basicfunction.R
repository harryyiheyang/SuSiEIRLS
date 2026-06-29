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
b=data.frame(Index=a,Variable=nam[summ$variable[indi]],CS=paste0("Main_CS",i),lbf=bb$cs_log10bf[bb$cs==i]*log(10),PIP=summ$variable_prob[indi])
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
    b=data.frame(Index=a,Variable=nam[summ$variable[indi]],CS=paste0("Env_CS",i),lbf=bb$cs_log10bf[bb$cs==i]*log(10),PIP=summ$variable_prob[indi])
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
b=data.frame(Index=a,Variable=namW[summ$variable[indi]],CS=paste0("Int_CS",i),lbf=bb$cs_log10bf[bb$cs==i]*log(10),PIP=summ$variable_prob[indi])
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
ProjPart = matrixMultiply(B, solve_with_ridge(BtB, BtA))
return(ProjPart)
}

solve_with_ridge <- function(A, B = NULL, ridge = 1e-8) {
  A <- as.matrix(A)
  out <- tryCatch(
    {
      if (is.null(B)) CppMatrix::matrixInverse(A) else CppMatrix::matrixSolve(A, as.matrix(B))
    },
    error = function(e) {
      A2 <- A
      diag(A2) <- diag(A2) + ridge
      if (is.null(B)) CppMatrix::matrixInverse(A2) else CppMatrix::matrixSolve(A2, as.matrix(B))
    }
  )
  out
}

weighted_residual_suffstats <- function(X, y, ZI, weights,
                                        n_threads = 1,
                                        ridge = 1e-8,
                                        block_size = 10000L) {
  n <- nrow(X)
  p <- ncol(X)
  q <- if (is.null(ZI)) 0L else ncol(ZI)
  block_size <- max(1L, as.integer(block_size))

  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  y <- as.numeric(y)

  XtX <- matrix(0, nrow = p, ncol = p)
  Xty <- numeric(p)
  yty <- 0

  if (q > 0L) {
    ZtZ <- matrix(0, nrow = q, ncol = q)
    ZtX <- matrix(0, nrow = q, ncol = p)
    Zty <- numeric(q)
  }

  for (start in seq.int(1L, n, by = block_size)) {
    end <- min(n, start + block_size - 1L)
    idx <- start:end
    wi <- weights[idx]
    if (!any(wi > 0)) next

    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]
    Xw <- Xb * wi

    XtX <- XtX + crossprod(Xb, Xw)
    Xty <- Xty + as.numeric(crossprod(Xb, wi * yb))
    yty <- yty + sum(wi * yb^2)

    if (q > 0L) {
      Zb <- ZI[idx, , drop = FALSE]
      Zw <- Zb * wi
      ZtZ <- ZtZ + crossprod(Zb, Zw)
      ZtX <- ZtX + crossprod(Zb, Xw)
      Zty <- Zty + as.numeric(crossprod(Zb, wi * yb))
      rm(Zb, Zw)
    }

    rm(Xb, Xw, yb, wi)
  }

  yty_raw <- yty
  if (q > 0L) {
    Zinv_ZtX <- solve_with_ridge(ZtZ, ZtX, ridge = ridge)
    Zinv_Zty <- solve_with_ridge(ZtZ, matrix(Zty, ncol = 1), ridge = ridge)

    XtX <- XtX - matrixMultiply(ZtX, Zinv_ZtX, transA = TRUE)
    Xty <- Xty - as.numeric(matrixMultiply(ZtX, Zinv_Zty, transA = TRUE))
    yty <- yty - as.numeric(crossprod(Zty, Zinv_Zty))
  }

  XtX <- (XtX + t(XtX)) / 2
  if (is.finite(yty) && yty < 0 &&
      yty > -sqrt(.Machine$double.eps) * max(1, abs(yty_raw))) {
    yty <- 0
  }

  dimnames(XtX) <- list(colnames(X), colnames(X))
  names(Xty) <- colnames(X)

  list(XtX = XtX, Xty = Xty, yty = as.numeric(yty))
}

clean_model_environment <- function(fit, env = .GlobalEnv) {
  if (is.null(fit)) return(fit)
  if (!is.null(fit$terms)) attr(fit$terms, ".Environment") <- env
  if (!is.null(fit$formula) && inherits(fit$formula, "formula")) {
    environment(fit$formula) <- env
  }
  if (!is.null(fit$call$formula) && inherits(fit$call$formula, "formula")) {
    environment(fit$call$formula) <- env
  }
  fit
}

clean_coef <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x
}

build_noncs_refit_term <- function(X, fitX, CSdt, cs_indices, XCS,
                                   noncs_var = 0.2) {
  if (is.null(fitX) || is.null(CSdt) || !length(cs_indices)) return(NULL)
  if (is.null(XCS) || ncol(as.matrix(XCS)) == 0L) return(NULL)

  beta_total <- clean_coef(stats::coef(fitX)[-1L])
  if (!length(beta_total) || length(beta_total) != ncol(X)) return(NULL)

  eta_total <- as.numeric(CppMatrix::matrixVectorMultiply(X, beta_total))
  var_total <- stats::var(eta_total)
  if (!is.finite(var_total) || var_total <= 1e-12) return(NULL)

  noncs_var <- as.numeric(noncs_var)[1L]
  if (!is.finite(noncs_var)) noncs_var <- 0.2
  noncs_var <- min(max(noncs_var, 0), 1)

  cs_parameter <- numeric(length(cs_indices))
  for (ii in seq_along(cs_indices)) {
    i_cs <- cs_indices[ii]
    vars_in_cs_i <- CSdt$variable[CSdt$cs == i_cs]
    vars_in_cs_i <- vars_in_cs_i[vars_in_cs_i >= 1L & vars_in_cs_i <= ncol(X)]
    if (!length(vars_in_cs_i)) next

    a <- as.numeric(fitX$alpha[i_cs, vars_in_cs_i])
    m <- as.numeric(fitX$mu[i_cs, vars_in_cs_i])
    ok <- is.finite(a) & is.finite(m) & a > 0
    if (any(ok)) {
      cs_parameter[ii] <- sum(a[ok] * abs(m[ok])) / sum(a[ok])
    }
  }

  eta_cs <- as.numeric(as.matrix(XCS) %*% cs_parameter)
  cs_var_ratio <- stats::var(eta_cs) / var_total
  eta_noncs <- eta_total - eta_cs
  if (!is.finite(cs_var_ratio) || cs_var_ratio >= (1 - noncs_var)) return(NULL)
  if (!is.finite(stats::sd(eta_noncs)) || stats::sd(eta_noncs) <= 1e-8) return(NULL)

  eta_noncs
}

safe_add_p <- function(idx, Coefmat) {
  if (is.null(idx)) return(NULL)
  if (is.data.frame(idx) && nrow(idx) == 0) return(idx)
  if (!("CS" %in% names(idx))) return(idx)
  if (is.null(Coefmat) || is.null(dim(Coefmat))) return(idx)
  cs <- as.character(idx$CS)
  pos <- match(cs, rownames(Coefmat))
  p   <- rep(NA_real_, length(cs))
  if (ncol(Coefmat) >= 4) {
    p[!is.na(pos)] <- Coefmat[pos[!is.na(pos)], 4]
  }

  idx$Pvalue <- p
  idx
}

################################################################################
robust_weight <- function(w, cutoff = 0.01) {
n <- length(w)
if (n == 0L || all(is.na(w))) return(rep(0, n))
w <- as.numeric(w)
w[!is.finite(w) | w < 0] <- NA
w_na <- is.na(w)

n_eff <- sum(!w_na)
if (n_eff <= 1L) {
  w[w_na] <- 0
  return(w)
}

if (n_eff < 1 / cutoff) {
lo <- min(w, na.rm = TRUE)
hi <- max(w, na.rm = TRUE)
} else {
lo <- quantile(w, probs = cutoff, na.rm = TRUE, names = FALSE, type = 7)
hi <- quantile(w, probs = 1 - cutoff, na.rm = TRUE, names = FALSE, type = 7)
}
w_trim <- pmin(pmax(w, lo), hi)
w_trim[w_na] <- 0
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

  if (identical(cor_method, "pearson")) {
    row_idx <- which(ok)
    r0 <- r[row_idx]
    r0 <- r0 - mean(r0)
    r_ss <- sum(r0^2)
    if (!is.finite(r_ss) || r_ss <= 0) return(NA_integer_)

    which_available <- which(available)
    scores <- rep(NA_real_, length(which_available))
    block_cols <- 64L

    for (start in seq.int(1L, length(which_available), by = block_cols)) {
      end <- min(length(which_available), start + block_cols - 1L)
      cols <- which_available[start:end]
      Xb <- X[row_idx, cols, drop = FALSE]

      if (any(!is.finite(Xb))) {
        for (j in seq_along(cols)) {
          xj <- Xb[, j]
          good <- is.finite(xj)
          if (sum(good) < 3L) next
          scores[start + j - 1L] <- suppressWarnings(stats::cor(
            xj[good], r0[good], method = "pearson"
          ))
        }
      } else {
        x_sum <- colSums(Xb)
        x_ss <- colSums(Xb^2) - x_sum^2 / length(r0)
        num <- as.numeric(crossprod(Xb, r0))
        good <- is.finite(x_ss) & (x_ss > 0)
        target <- start:end
        scores[target[good]] <- num[good] / sqrt(x_ss[good] * r_ss)
      }

      rm(Xb)
    }

    scores[!is.finite(scores)] <- NA_real_
    if (all(is.na(scores))) return(NA_integer_)
    return(which_available[which.max(abs(scores))])
  }

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
  U <- as.numeric(matrixMultiply(tilde_G0, matrix(tilde_y, ncol = 1), transA = TRUE))
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
