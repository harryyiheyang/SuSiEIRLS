tptn_evaulate=function(true_main_index,hat_beta){
  select_index=which(hat_beta!=0)
tp_main=ifelse(length(setdiff(true_main_index,select_index))==0,1,0)
tn_main=ifelse(length(setdiff(select_index,true_main_index))==0,1,0)
g=data.frame(tp=tp_main,tn=tn_main)
return(g)
}

scale_to_target_var <- function(x, target_var) {
  x <- as.numeric(scale(x, center = TRUE, scale = FALSE))
  vx <- stats::var(x)
  if (!is.finite(vx) || vx <= 0) {
    stop("Cannot scale a zero-variance effect.")
  }
  x * sqrt(target_var / vx)
}

make_eta_components <- function(X, Z, beta, alpha, total_h2 = 0.3,
                                z_to_x_ratio = 2) {
  x_var <- total_h2 / (1 + z_to_x_ratio)
  z_var <- total_h2 - x_var
  etaX <- scale_to_target_var(as.numeric(X %*% beta), x_var)
  etaZ <- scale_to_target_var(as.numeric(Z %*% alpha), z_var)
  eta <- etaZ + etaX
  veta <- stats::var(eta)
  if (!is.finite(veta) || veta <= 0) {
    stop("Cannot scale a zero-variance linear predictor.")
  }
  common_scale <- sqrt(total_h2 / veta)
  etaZ <- etaZ * common_scale
  etaX <- etaX * common_scale
  list(etaZ = etaZ, etaX = etaX, eta = etaZ + etaX)
}

scale_logit_liability_h2 <- function(eta, h2 = 0.3) {
  target_var <- h2 / (1 - h2) * (pi^2 / 3)
  scale_to_target_var(eta, target_var)
}

scale_cox_liability_h2 <- function(eta, h2 = 0.3) {
  target_var <- h2 / (1 - h2) * (pi^2 / 6)
  scale_to_target_var(eta, target_var)
}

quiet_eval <- function(expr) {
  capture.output(value <- eval.parent(substitute(expr)))
  value
}

susie_to_main_index_x_only <- function(fit_susie, X_aug, p, coverage = 0.95,
                                       min_abs_cor = 0.1) {
  class(fit_susie) <- c("susie", "list")
  cs_out <- susieR::susie_get_cs(
    fit_susie, X = X_aug, coverage = coverage, min_abs_cor = min_abs_cor
  )
  cs_list <- lapply(cs_out$cs, function(idx) idx[idx <= p])
  cs_list <- cs_list[vapply(cs_list, length, integer(1)) > 0]
  if (!length(cs_list)) return(data.frame())

  idx <- unlist(cs_list)
  out <- data.frame(
    Index = idx,
    Variable = colnames(X_aug)[idx],
    CS = rep(paste0("Main_CS", seq_along(cs_list)),
             times = vapply(cs_list, length, integer(1))),
    stringsAsFactors = FALSE
  )
  out$PIP <- fit_susie$pip[out$Index]
  row.names(out) <- NULL
  out
}

cs_contains_truth <- function(main_index, true_idx) {
  if (is.null(main_index) || !nrow(main_index)) {
    return(list(power = 0, false_cs = 0, n_cs = 0, detected = integer(0)))
  }
  detected <- sort(unique(main_index$Index))
  cs_has_signal <- vapply(
    split(main_index$Index, main_index$CS),
    function(idx) any(idx %in% true_idx),
    logical(1)
  )
  list(
    power = mean(true_idx %in% detected),
    false_cs = mean(!cs_has_signal),
    n_cs = length(cs_has_signal),
    detected = detected
  )
}

benchmark_summary <- function(per_run) {
  summary <- aggregate(
    cbind(power, false_cs, n_cs, time_sec) ~ method,
    data = per_run,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(summary) <- c("method", "mean_power", "mean_false_cs",
                      "mean_n_cs", "mean_time_sec")
  failed <- aggregate(
    is.na(error) ~ method,
    data = per_run,
    FUN = function(x) sum(!x)
  )
  names(failed)[2] <- "n_failed"
  merge(summary, failed, by = "method", all.x = TRUE)
}

