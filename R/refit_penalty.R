.refit_penalty_variance <- function(fitX, cs_indices, term_names,
                                    noncs_V_default = 2) {
  term_names <- as.character(term_names)
  if (!length(term_names)) return(stats::setNames(numeric(0), character(0)))

  V <- as.numeric(fitX$V)
  if (!length(V) && any(term_names != "Main_noncs_res")) {
    stop("The SuSiE fit does not contain prior variances V.")
  }

  out <- rep(NA_real_, length(term_names))
  names(out) <- term_names
  cs_names <- if (length(cs_indices)) paste0("Main_CS", cs_indices) else character(0)
  pos <- match(term_names, cs_names)
  is_cs <- !is.na(pos)
  if (any(is_cs)) out[is_cs] <- V[cs_indices[pos[is_cs]]]

  is_noncs <- term_names == "Main_noncs_res"
  if (any(is_noncs)) {
    positive <- V[is.finite(V) & V > 0]
    if (!length(positive)) positive <- noncs_V_default
    out[is_noncs] <- max(positive)
  }

  unknown <- !is_cs & term_names != "Main_noncs_res"
  if (any(unknown)) {
    stop("Unknown penalized refit term: ", paste(term_names[unknown], collapse = ", "), ".")
  }
  keep <- is.finite(out) & out > 0
  invalid <- !keep
  if (any(invalid)) {
    bad <- names(out)[invalid]
    stop(
      "A positive finite SuSiE V is required for penalized refit term(s): ",
      paste(bad, collapse = ", "), "."
    )
  }
  out[keep]
}

.mgcv_refit_dispersion <- function(fit) {
  phi <- as.numeric(fit$sig2)[1L]
  if (!length(phi) || !is.finite(phi) || phi <= 0) {
    phi <- as.numeric(summary(fit)$dispersion)[1L]
  }
  if (!length(phi) || !is.finite(phi) || phi <= 0) {
    stop("The previous mgcv fit must provide a positive finite dispersion.")
  }
  phi
}

.mgcv_fit_fixed_ridge <- function(response, rhs, data, family, penalty_V,
                                  dispersion = 1, mgcv_model = NULL) {
  unpenalized_terms <- attr(penalty_V, "unpenalized_terms")
  penalty_names <- names(penalty_V)
  penalty_V <- as.numeric(penalty_V)
  names(penalty_V) <- penalty_names
  if (!length(penalty_V)) {
    fit <- .mgcv_fit_explicit(
      response, rhs, data, family, mgcv_model = mgcv_model
    )
    if (length(unpenalized_terms)) {
      attr(fit, "refit_penalty") <- list(
        V = numeric(0), precision = numeric(0),
        working_precision = numeric(0), dispersion = dispersion, sp = 1,
        unpenalized_terms = unpenalized_terms
      )
    }
    return(fit)
  }
  if (is.null(penalty_names) || any(!nzchar(penalty_names)) ||
      anyDuplicated(penalty_names)) {
    stop("penalty_V must be a uniquely named vector.")
  }
  if (any(!is.finite(penalty_V) | penalty_V <= 0)) {
    stop("penalty_V must contain positive finite SuSiE variances.")
  }
  if (!is.numeric(dispersion) || length(dispersion) != 1L ||
      !is.finite(dispersion) || dispersion <= 0) {
    stop("dispersion must be a positive finite scalar.")
  }
  if (!all(penalty_names %in% rhs) || !all(penalty_names %in% names(data))) {
    stop("Every penalized refit term must occur in both rhs and data.")
  }

  X_pen <- as.matrix(data[, penalty_names, drop = FALSE])
  dat <- data[, setdiff(names(data), penalty_names), drop = FALSE]
  dat$X_pen <- I(X_pen)
  rhs <- c(setdiff(rhs, penalty_names), "X_pen")
  PP <- list(X_pen = list(
    diag(dispersion / penalty_V, nrow = length(penalty_V),
         ncol = length(penalty_V)),
    sp = 1
  ))

  engine <- .mgcv_fit_engine(nrow(dat), mgcv_model)
  family <- .mgcv_patch_family_environment(family)
  fit <- engine$fit(
    .mgcv_explicit_formula(response, rhs), data = dat,
    family = family, method = engine$method, paraPen = PP
  )

  coef_names <- names(stats::coef(fit))
  bundled <- grep("^X_pen", coef_names)
  if (length(bundled) != length(penalty_names)) {
    stop("mgcv did not preserve every penalized refit coefficient.")
  }
  names(fit$coefficients)[bundled] <- penalty_names
  attr(fit, "refit_penalty") <- list(
    V = stats::setNames(penalty_V, penalty_names),
    precision = stats::setNames(1 / penalty_V, penalty_names),
    working_precision = stats::setNames(
      dispersion / penalty_V, penalty_names
    ),
    dispersion = dispersion,
    sp = 1,
    unpenalized_terms = unpenalized_terms
  )
  fit
}
