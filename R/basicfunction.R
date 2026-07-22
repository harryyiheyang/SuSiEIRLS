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
solve_with_ridge <- function(A, B = NULL, ridge = 1e-8) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (is.finite(ridge) && ridge > 0) {
    diag(A) <- diag(A) + ridge
  }
  if (is.null(B)) CppMatrix::matrixInverse(A) else CppMatrix::matrixSolve(A, as.matrix(B))
}

make_diagnostics <- function(iterations, eps, start_time) {
  final_eps <- if (length(eps)) as.numeric(utils::tail(eps, 1L)) else NA_real_
  data.frame(
    iterations = as.integer(iterations),
    eps = final_eps,
    runtime_seconds = unname(proc.time()[["elapsed"]] - start_time)
  )
}

weighted_residual_suffstats <- function(X, y, ZI, weights,
                                        n_threads = 1,
                                        ridge = 1e-8,
                                        block_size = 10000L,
                                        projection_solver = c("ridge", "ginv"),
                                        ginv_tol = sqrt(.Machine$double.eps)) {
  if (!is.null(ZI)) ZI <- as.matrix(ZI)
  q <- if (is.null(ZI)) 0L else ncol(ZI)
  block_size <- max(1L, as.integer(block_size))
  projection_solver <- match.arg(projection_solver)

  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  y <- as.numeric(y)

  tilde_X <- X * sqrt(weights)
  XtX <- SuSiE4I::blockwise_crossprod(
    tilde_X, n_threads = n_threads, block_size = block_size
  )
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

    if (projection_solver == "ginv") {
      Zinv <- MASS::ginv(as.matrix(ZtZ), tol = ginv_tol)
      Zinv_ZtX <- CppMatrix::matrixMultiply(Zinv, ZtX)
      Zinv_Zty <- CppMatrix::matrixMultiply(Zinv, matrix(Zty, ncol = 1))
    } else {
      Zinv_ZtX <- solve_with_ridge(ZtZ, ZtX, ridge = ridge)
      Zinv_Zty <- solve_with_ridge(ZtZ, matrix(Zty, ncol = 1), ridge = ridge)
    }

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

.susie_default_para <- function() {
  list(
    standardize = FALSE,
    scaled_prior_variance = 2,
    estimate_residual_variance = TRUE,
    residual_variance = 0.5,
    residual_variance_lowerbound = 0.1,
    residual_variance_upperbound = 1.01,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    max_iter = 300,
    coverage = 0.9
  )
}

.validate_susie_para <- function(x, arg = "susie_para") {
  if (is.null(x)) return(list())
  if (!is.list(x)) stop(arg, " must be NULL or a named list.")
  if (!length(x)) return(list())
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop(arg, " must be a named list.")
  }
  if (anyDuplicated(nm)) stop(arg, " must not contain duplicate names.")

  protected <- c("XtX", "Xty", "yty", "n", "L")
  blocked <- intersect(nm, protected)
  if (length(blocked)) {
    stop(
      arg, " cannot set structural SuSiE inputs: ",
      paste(blocked, collapse = ", "), "."
    )
  }

  if (all(c("prior_variance", "scaled_prior_variance") %in% nm)) {
    stop(arg, " cannot contain both prior_variance and scaled_prior_variance.")
  }
  if ("prior_variance" %in% nm) {
    V <- x$prior_variance
    if (!is.numeric(V) || length(V) != 1L || !is.finite(V) || V <= 0) {
      stop(arg, "$prior_variance must be a positive finite numeric scalar.")
    }
  }
  if ("scaled_prior_variance" %in% nm) {
    V <- x$scaled_prior_variance
    if (!is.numeric(V) || length(V) != 1L || !is.finite(V) || V <= 0) {
      stop(arg, "$scaled_prior_variance must be a positive finite numeric scalar.")
    }
    warning(
      "scaled_prior_variance is interpreted as an absolute coefficient prior variance and re-scaled at every IRLS iteration; use prior_variance instead.",
      call. = FALSE
    )
    x$prior_variance <- V
    x$scaled_prior_variance <- NULL
    nm <- names(x)
  }

  valid <- c(names(formals(susieR::susie_ss)), "prior_variance")
  unknown <- setdiff(nm, valid)
  if (length(unknown)) {
    stop(
      "Unknown susieR::susie_ss parameter in ", arg, ": ",
      paste(unknown, collapse = ", "), "."
    )
  }
  x
}

.resolve_susie_para <- function(susie_para = NULL) {
  .validate_susie_para(susie_para)
}

.susie_iteration_args <- function(susie_para, structural, iter, min.iter) {
  args <- .susie_default_para()
  overrides <- susie_para[!vapply(susie_para, is.null, logical(1))]
  if (iter <= min.iter && length(overrides)) {
    warm_V_controls <- c(
      "prior_variance", "scaled_prior_variance", "estimate_prior_variance"
    )
    overrides <- overrides[setdiff(names(overrides), warm_V_controls)]
  }
  if (length(overrides)) args[names(overrides)] <- overrides
  if (iter <= min.iter) args$estimate_prior_variance <- FALSE

  epv <- args$estimate_prior_variance
  if (!is.logical(epv) || length(epv) != 1L || is.na(epv)) {
    stop("susie_para$estimate_prior_variance must be TRUE or FALSE.")
  }
  if ("prior_variance" %in% names(args)) {
    y_scale <- structural$yty / (structural$n - 1)
    if (!is.numeric(y_scale) || length(y_scale) != 1L ||
        !is.finite(y_scale) || y_scale <= 0) {
      stop("structural$yty / (structural$n - 1) must be positive and finite to convert prior_variance.")
    }
    args$scaled_prior_variance <- args$prior_variance / y_scale
    args$prior_variance <- NULL
  }
  args[names(structural)] <- structural
  args
}

validate_noncs_max_abs_cor <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("noncs_max_abs_cor must be a positive finite numeric scalar.")
  }
  min(as.numeric(x), 0.9)
}

noncs_correlation_ok <- function(term, design, max_abs_cor = 0.9) {
  if (is.null(design)) return(TRUE)

  term <- as.numeric(term)
  design <- as.matrix(design)
  if (!length(term) || nrow(design) != length(term)) return(FALSE)
  if (any(!is.finite(term)) || any(!is.finite(design))) return(FALSE)
  if (!ncol(design)) return(TRUE)

  term0 <- term - mean(term)
  term_ss <- sum(term0^2)
  if (!is.finite(term_ss) || term_ss <= 1e-12) return(FALSE)

  design0 <- sweep(design, 2L, colMeans(design), "-")
  design_ss <- colSums(design0^2)
  keep <- is.finite(design_ss) & design_ss > 1e-12
  if (!any(keep)) return(TRUE)

  cors <- as.numeric(crossprod(design0[, keep, drop = FALSE], term0)) /
    sqrt(design_ss[keep] * term_ss)
  max_abs_cor <- validate_noncs_max_abs_cor(max_abs_cor)
  !any(is.finite(cors) & abs(cors) >= max_abs_cor)
}

build_noncs_refit_term <- function(X, fitX, CSdt, cs_indices, XCS,
                                   noncs_var = 0.1,
                                   noncs_max_abs_cor = 0.9,
                                   cor_design = NULL) {
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

  gate_design <- if (is.null(cor_design)) XCS else cbind(XCS, cor_design)
  if (!noncs_correlation_ok(
    eta_noncs, gate_design, max_abs_cor = noncs_max_abs_cor
  )) return(NULL)

  eta_noncs
}

build_no_cs_noncs_refit_term <- function(X, fitX, cor_design = NULL,
                                         noncs_max_abs_cor = 0.9) {
  if (is.null(fitX)) return(NULL)
  if (!length(fitX$V)) return(NULL)

  beta_total <- clean_coef(stats::coef(fitX)[-1L])
  if (length(beta_total) != ncol(X)) {
    stop("The SuSiE coefficient vector does not match ncol(X).")
  }

  noncs_res <- as.numeric(CppMatrix::matrixVectorMultiply(X, beta_total))
  if (length(noncs_res) != nrow(X)) {
    stop("The no-CS rescue term does not match nrow(X).")
  }
  noncs_res[!is.finite(noncs_res)] <- 0

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
  p_col <- match(
    c("p", "pr(>|z|)", "pr(>|t|)"), tolower(colnames(Coefmat)),
    nomatch = 0L
  )
  p_col <- p_col[p_col > 0L]
  if (!length(p_col) && ncol(Coefmat) >= 4L) p_col <- 4L
  if (length(p_col)) {
    p[!is.na(pos)] <- Coefmat[pos[!is.na(pos)], p_col[1L]]
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

fit_init_cox <- function(X, y, status, Z, selected) {
  # The warm start precedes SuSiE and contains no Main_CS/non-CS refit terms.
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

greedy_cox_warm_start <- function(X, y, status, Z, L.init = 1) {
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
