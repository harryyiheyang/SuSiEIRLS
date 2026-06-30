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
return(NULL)
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
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (is.finite(ridge) && ridge > 0) {
    diag(A) <- diag(A) + ridge
  }
  if (is.null(B)) CppMatrix::matrixInverse(A) else CppMatrix::matrixSolve(A, as.matrix(B))
}

weighted_residual_suffstats <- function(X, y, ZI, weights,
                                        n_threads = 1,
                                        ridge = 1e-8,
                                        block_size = 10000L) {
  if (!is.null(ZI)) ZI <- as.matrix(ZI)
  q <- if (is.null(ZI)) 0L else ncol(ZI)
  block_size <- max(1L, as.integer(block_size))

  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  y <- as.numeric(y)

  tilde_X <- X * sqrt(weights)
  XtX <- blockwise_crossprod(tilde_X, n_threads = n_threads,
                             block_size = block_size)
  rm(tilde_X)
  gc(FALSE)

  wy <- weights * y
  Xty <- as.numeric(matrixMultiply(X, matrix(wy, ncol = 1), transA = TRUE))
  yty <- sum(weights * y^2)
  yty_raw <- yty
  if (q > 0L) {
    Zw <- ZI * weights
    ZtZ <- matrixMultiply(ZI, Zw, transA = TRUE)
    ZtX <- matrixMultiply(Zw, X, transA = TRUE)
    Zty <- as.numeric(matrixMultiply(ZI, matrix(wy, ncol = 1), transA = TRUE))
    rm(Zw)

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

  eta_x <- as.numeric(CppMatrix::matrixVectorMultiply(X, beta_total))
  var_eta_x <- stats::var(eta_x)
  if (!is.finite(var_eta_x) || var_eta_x <= 1e-12) return(NULL)

  noncs_var <- as.numeric(noncs_var)[1L]
  if (!is.finite(noncs_var)) noncs_var <- 0.2
  noncs_var <- min(max(noncs_var, 0), 1)

  XCS <- as.matrix(XCS)
  if (nrow(XCS) != length(eta_x)) return(NULL)
  if (any(!is.finite(XCS)) || any(!is.finite(eta_x))) return(NULL)

  keep <- apply(XCS, 2, function(z) {
    stats::sd(as.numeric(z)) > 1e-8
  })
  if (!any(keep)) return(NULL)
  XCS <- XCS[, keep, drop = FALSE]

  X_full <- cbind(Intercept = 1, XCS)
  XtX <- crossprod(X_full)
  if (qr(XtX)$rank < ncol(XtX)) return(NULL)

  proj_coef <- solve(XtX, crossprod(X_full, eta_x))
  if (any(!is.finite(proj_coef))) return(NULL)

  eta_noncs <- eta_x - as.numeric(X_full %*% proj_coef)
  var_noncs <- stats::var(eta_noncs)
  if (!is.finite(var_noncs) || var_noncs <= 1e-12) return(NULL)
  if (var_noncs / var_eta_x < noncs_var) return(NULL)

  eta0 <- eta_noncs - mean(eta_noncs)
  x0 <- sweep(XCS, 2, colMeans(XCS), "-")
  denom <- sqrt(sum(eta0^2) * colSums(x0^2))
  cors <- as.numeric(crossprod(x0, eta0)) / denom
  if (any(is.finite(cors) & abs(cors) >= 0.999)) return(NULL)

  eta_noncs
}

build_no_cs_noncs_refit_term <- function(X, fitX) {
  if (is.null(fitX)) return(NULL)

  beta_total <- clean_coef(stats::coef(fitX)[-1L])
  if (!length(beta_total) || length(beta_total) != ncol(X)) return(NULL)

  noncs_res <- as.numeric(CppMatrix::matrixVectorMultiply(X, beta_total))
  if (length(noncs_res) != nrow(X)) return(NULL)
  if (!is.finite(stats::sd(noncs_res)) || stats::sd(noncs_res) <= 1e-8) {
    return(NULL)
  }

  noncs_res
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

validate_suff_block_size <- function(suff_block_size) {
  if (!is.numeric(suff_block_size) || length(suff_block_size) != 1L ||
      !is.finite(suff_block_size) || suff_block_size < 1) {
    stop("suff_block_size must be a positive numeric scalar.")
  }
  as.integer(suff_block_size)
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

select_by_residual_score <- function(X, residual, available) {
  if (!any(available)) return(NA_integer_)

  r <- as.numeric(residual)
  ok <- is.finite(r)
  if (!any(ok)) return(NA_integer_)
  r[!ok] <- 0
  scores <- as.numeric(CppMatrix::matrixMultiply(
    X, matrix(r, ncol = 1), transA = TRUE
  ))
  scores[!available] <- NA_real_
  scores[!is.finite(scores)] <- NA_real_
  if (all(is.na(scores))) return(NA_integer_)
  which.max(abs(scores))
}

greedy_glm_warm_start <- function(X, y, Z, family, L.init = 1,
                                  init_cor_method = NULL) {
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- fit_init_glm(X = X, y = y, Z = Z, selected = selected, family = family)

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "response")
    j <- select_by_residual_score(X = X, residual = r, available = available)
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    fit <- fit_init_glm(X = X, y = y, Z = Z, selected = selected, family = family)
  }

  fit
}

greedy_nb_warm_start <- function(X, y, Z, L.init = 1, theta_init, estimate_theta,
                                 init_cor_method = NULL) {
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
    j <- select_by_residual_score(X = X, residual = r, available = available)
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
                                  init_cor_method = NULL) {
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- fit_init_cox(X = X, y = y, status = status, Z = Z, selected = selected)

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "martingale")
    j <- select_by_residual_score(X = X, residual = r, available = available)
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
