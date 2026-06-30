# SuSiEIRLS

SuSiE fine-mapping for non-Gaussian outcomes via IRLS.

The package extends the Sum of Single Effects (SuSiE) framework to generalized linear models (GLMs), binary outcomes with a Pólya-Gamma augmentation path, and Cox proportional-hazards survival outcomes. At each outer iteration the algorithm constructs approximate Gaussian sufficient statistics from the current working model and passes them to `susieR::susie_ss()`.

## Installation

```r
# install.packages("devtools")
devtools::install_github("harryyiheyang/SuSiEIRLS")
```

**Dependencies:** `susieR`, `mgcv`, `survival`, `MASS`, `CppMatrix`, `SuSiE4I`.

## Main function

`SuSiE_IRLS()` is the user-facing wrapper. It dispatches to one of three internal paths based on the response type and family:

| Response | Family / method | Sufficient statistics |
|---|---|---|---|
| Continuous or count | Any `stats` or `mgcv` family object (e.g. `poisson()`, `Gamma()`, `mgcv::tw()`, `mgcv::betar()`, `mgcv::nb()`) | IRLS working response and weights |
| Binary | `binomial(link = "logit")` with `logit_method = "pg"` (default) | | Pólya-Gamma augmented pseudo-response |
| Binary | `binomial(link = "logit")` with `logit_method = "glm"` | Standard IRLS working response |
| Survival | Pass a `survival::Surv` object as `y` | Cox partial-likelihood score and information |

### Basic usage

```r
library(SuSiEIRLS)

## Binary (Pólya-Gamma path, the default)
fit <- SuSiE_IRLS(X = X, Z = Z, y = y, L = 10)

## Binary (standard IRLS path)
fit <- SuSiE_IRLS(X = X, Z = Z, y = y, L = 10, logit_method = "glm")

## Poisson
fit <- SuSiE_IRLS(X = X, Z = Z, y = y, family = poisson(), L = 10)

## Negative binomial (mgcv parameterisation)
fit <- SuSiE_IRLS(X = X, Z = Z, y = y, family = "negbin", L = 10)

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
| `coverage` | 0.9 | Credible-set coverage |
| `max.iter` | 15 | Maximum outer IRLS iterations |
| `logit_method` | `"pg"` | `"pg"` for Pólya-Gamma, `"glm"` for standard IRLS (binary only) |

### Output

A list containing:

- `fitX` — the SuSiE fit object from the final iteration.
- `fitJoint` — the refitted GLM or Cox model with selected variables.
- `main_index` — summary table of credible sets with PIPs and p-values.
- `JointCoef` — coefficient table from the final joint model.
- `converged` — logical convergence flag.
- `iter` — number of outer iterations completed.

## Method overview

**GLM / extended-GLM path.** The algorithm uses IRLS to linearise the GLM likelihood around the current estimate. At each outer iteration it forms a working response $z$ and diagonal weight matrix $W$, projects out covariates $Z$, and constructs weighted sufficient statistics $(X^\top W X,\; X^\top W z,\; z^\top W z)$ for `susie_ss()`. The residual variance is estimated within a bounded interval (default $[0.1, 1]$).

**Binary (Pólya-Gamma) path.** For binary logistic outcomes, the Pólya-Gamma data-augmentation identity represents the logistic likelihood as a conditionally Gaussian problem with augmentation weights $\omega_i = \tanh(\eta_i/2) / (2\eta_i)$. This path improves signal detection when the standard IRLS approximation is weak. The final log Bayes factors are corrected back to the IRLS scale.

**Cox path.** For survival outcomes, no working response exists. The package instead constructs score-based sufficient statistics $(X^\top M,\; A - B^\top B)$ from the Cox partial likelihood using a single-pass Breslow accumulator implemented in C++ (`cox_suffstat.cpp`). These replace $X^\top y$ and $X^\top X$ in the SuSiE-SS call.

## License

MIT
