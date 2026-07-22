# SuSiEIRLS

SuSiE fine-mapping for non-Gaussian outcomes via IRLS.

The package extends the Sum of Single Effects (SuSiE) framework to generalized linear models (GLMs), ordered outcomes, and Cox proportional-hazards survival outcomes. At each outer iteration the algorithm constructs approximate Gaussian sufficient statistics from the current working model and passes them to `susieR::susie_ss()`.

## Installation

```r
# install.packages("devtools")
devtools::install_github("harryyiheyang/SuSiEIRLS")
```

**Dependencies:** `susieR`, `mgcv`, `ordinal`, `survival`, `MASS`, `CppMatrix`, `SuSiE4I`.

## Main function

`SuSiE_IRLS()` is the user-facing wrapper. It dispatches to one of three internal paths based on the response type and family:

| Response | Family / method | Sufficient statistics |
|---|---|---|
| Continuous or count | Any `stats` or `mgcv` family object (e.g. `poisson()`, `Gamma()`, `mgcv::tw()`, `mgcv::betar()`, `mgcv::nb()`) | IRLS working response and weights |
| Binary | `binomial(link = "logit")` | Standard IRLS working response |
| Ordered categorical | `family = "ocat"` or `mgcv::ocat(R = )` | Cumulative-logit score and observed information |
| Survival | Pass a `survival::Surv` object as `y` | Cox partial-likelihood score and information |

### Basic usage

```r
library(SuSiEIRLS)

## Binary logistic IRLS
fit <- SuSiE_IRLS(X = X, Z = Z, y = y, L = 10)

## Poisson
fit <- SuSiE_IRLS(X = X, Z = Z, y = y, family = poisson(), L = 10)

## Negative binomial (mgcv parameterisation)
fit <- SuSiE_IRLS(X = X, Z = Z, y = y,
                  family = mgcv::nb(theta = NULL), L = 10)

## Zero-inflated Poisson; ziP parameters belong in the family constructor
fit <- SuSiE_IRLS(X = X, Z = Z, y = y,
                  family = mgcv::ziP(theta = NULL, b = 0), L = 10)

## Ordered categorical cumulative-logit outcome
fit <- SuSiE_IRLS(X = X, Z = Z, y = ordered_y,
                  family = mgcv::ocat(R = 4), L = 10)

## Cox proportional hazards
library(survival)
fit <- SuSiE_IRLS(X = X, Z = Z, y = Surv(time, status), L = 10)
```

### Key arguments

| Argument | Default | Description |
|---|---|---|
| `X` | — | $n \times p$ predictor matrix |
| `y` | — | Response vector, or `Surv` object for Cox |
| `Z` | `NULL` | $n \times q$ covariate matrix (projected out before SuSiE) |
| `family` | `binomial(link = "logit")` | GLM or `mgcv` family; ignored when `y` is `Surv` |
| `L` | 10 | Number of single effects |
| `susie_para` | `NULL` | Named list of optional `susieR::susie_ss()` overrides |
| `max.iter` | 10 | Maximum outer IRLS iterations |

For example, use
`susie_para = list(estimate_residual_variance = FALSE, coverage = 0.95)`.
The sufficient statistics (`XtX`, `Xty`, `yty`, and `n`) and `L` are managed by
SuSiEIRLS and cannot be overridden through this list. Parameters omitted from
`susie_para` retain the current SuSiEIRLS settings when the package already has
one; otherwise they use the native `susieR::susie_ss()` default. Use
`prior_variance` for an absolute coefficient prior variance. The legacy
`scaled_prior_variance` name is accepted with a warning and interpreted the
same way.

### Output

A list containing:

- `fitX` — the SuSiE fit object from the final iteration.
- `fitJoint` — the refitted GLM or Cox model with selected variables.
- `discovery_summary` — discovered main effects with PIPs and p-values.
- `diagnostics` — iterations, final convergence error, and runtime.

## Method overview

**GLM / extended-GLM path.** The algorithm uses IRLS to linearise the GLM likelihood around the current estimate. At each outer iteration it forms a working response $z$ and diagonal weight matrix $W$, projects out covariates $Z$, and constructs weighted sufficient statistics $(X^\top W X,\; X^\top W z,\; z^\top W z)$ for `susie_ss()`. The residual variance is estimated within a bounded interval (default $[0.1, 1.01]$).

**Ordered categorical path.** For ordinal logit outcomes, the package uses `mgcv::ocat()`; explicit non-logit `clm_*` families use `ordinal::clm()`. It builds a local quadratic from the score and observed information for the location parameter, projects out covariates and threshold nuisance parameters, and passes the resulting sufficient statistics to `susie_ss()`. The default SuSiE residual variance is initialized at 0.5 and estimated within $[0.1, 1.01]$.

**Cox path.** For survival outcomes, no working response exists. The package instead constructs score-based sufficient statistics $(X^\top M,\; A - B^\top B)$ from the Cox partial likelihood using a single-pass Breslow accumulator implemented in C++ (`cox_suffstat.cpp`). These replace $X^\top y$ and $X^\top X$ in the SuSiE-SS call.

## License

MIT
