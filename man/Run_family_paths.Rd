\name{Run_Binary}
\alias{Run_Binary}
\alias{Run_GLM}
\alias{Run_NB}
\alias{Run_Cox}
\title{Family-specific IRLS-SuSiE paths}
\usage{
Run_Binary(
  X,
  y,
  Z = NULL,
  family = binomial(link = "logit"),
  weight_cutoff = 0.005,
  L,
  max.iter,
  min.iter,
  max.eps,
  susie.iter,
  pip.thres = 5e-3,
  verbose = TRUE,
  n_threads = 1,
  coverage = 0.95,
  estimate_residual_variance = FALSE,
  scaled_prior_variance = 1,
  residual_variance = 1,
  L.init = 1,
  init_cor_method = c("pearson", "spearman"),
  ...
)

Run_GLM(
  X,
  y,
  Z = NULL,
  weight_cutoff = 0.005,
  family = binomial(link = "logit"),
  L,
  max.iter,
  min.iter,
  max.eps,
  susie.iter,
  verbose = TRUE,
  n_threads = 1,
  coverage = 0.95,
  estimate_residual_variance = FALSE,
  residual_variance = 1,
  scaled_prior_variance = 1,
  pip.thres = 0.005,
  L.init = 1,
  init_cor_method = c("pearson", "spearman"),
  ...
)

Run_NB(
  X,
  y,
  Z = NULL,
  weight_cutoff = 0.005,
  L,
  max.iter,
  min.iter,
  max.eps,
  susie.iter,
  verbose = TRUE,
  n_threads = 1,
  coverage = 0.95,
  estimate_residual_variance = FALSE,
  residual_variance = 1,
  scaled_prior_variance = 1,
  pip.thres = 0.005,
  theta_init = 10,
  estimate_theta = TRUE,
  L.init = 1,
  init_cor_method = c("pearson", "spearman"),
  ...
)

Run_Cox(
  X,
  y,
  status,
  Z = NULL,
  L,
  max.iter,
  min.iter,
  max.eps,
  susie.iter,
  pip.thres = 5e-3,
  verbose = TRUE,
  n_threads = 1,
  coverage = 0.95,
  estimate_residual_variance = FALSE,
  scaled_prior_variance = 1,
  residual_variance = 1,
  ridge = 1e-6,
  L.init = 1,
  init_cor_method = c("pearson", "spearman"),
  ...
)
}
\arguments{
\item{X}{An n by p numeric matrix of predictors.}
\item{y}{Response vector. For \code{Run_Cox}, event or follow-up time.}
\item{status}{Event indicator for \code{Run_Cox}.}
\item{Z}{Optional covariate matrix. Covariates are adjusted for and are not treated as SuSiE candidate effects.}
\item{family}{A GLM family object.}
\item{weight_cutoff}{Quantile used to clip extreme IRLS weights.}
\item{L}{Number of single effects in SuSiE.}
\item{L.init}{Number of SNPs used in the initial low-dimensional warm start.}
\item{max.iter}{Maximum number of outer IRLS-SuSiE iterations.}
\item{min.iter}{Minimum number of iterations before convergence checking.}
\item{max.eps}{Convergence threshold.}
\item{susie.iter}{Maximum number of iterations for each SuSiE fit.}
\item{pip.thres}{PIP threshold used when filtering credible-set effects for joint refitting.}
\item{verbose}{Whether to print progress and plot convergence traces.}
\item{n_threads}{Number of threads used by internal matrix operations.}
\item{coverage}{Credible-set coverage level.}
\item{estimate_residual_variance}{Whether SuSiE estimates residual variance.}
\item{scaled_prior_variance}{Initial scaled prior variance used by SuSiE.}
\item{residual_variance}{Fixed residual variance when not estimated.}
\item{init_cor_method}{Correlation method for greedy low-dimensional warm starts.}
\item{theta_init}{Initial or fixed negative-binomial dispersion parameter.}
\item{estimate_theta}{Whether to estimate the negative-binomial dispersion parameter.}
\item{ridge}{Diagonal ridge added to Cox information matrices.}
\item{...}{Additional arguments passed to the SuSiE fitting routine.}
}
\description{
Advanced family-specific entry points used by \code{\link{SuSiE_IRLS}}. These
functions expose the binary-logistic, general-GLM, negative-binomial,
and Cox paths directly for benchmarking and method comparisons.
}
