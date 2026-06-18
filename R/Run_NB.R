#' Negative-binomial IRLS-SuSiE path
#' @export
Run_NB=function(X,y,Z=NULL,weight_cutoff=0.005,
                L,max.iter,min.iter,max.eps,susie.iter,
                verbose=TRUE,n_threads=1,coverage=0.95,
                estimate_residual_variance=FALSE,
                residual_variance=1,scaled_prior_variance=1,
                pip.thres=0.005,
                theta_init=10,
                estimate_theta=TRUE,
                L.init = 1,
                init_cor_method = c("pearson", "spearman"),
                ...) {

n=n_eff=length(y)
p=ncol(X)
init_cor_method <- match.arg(init_cor_method)

if (is.null(Z)) {
Z=matrix(nrow=n,ncol=0)
ZI=matrix(1,nrow=n,ncol=1)
colnames(ZI)="Intercept"
has_covariates=FALSE
} else {
if (is.null(dim(Z))) Z=matrix(Z,ncol=1)
if (is.null(colnames(Z))) colnames(Z)=paste0("Z",seq_len(ncol(Z)))
ZI=cbind(1,Z)
colnames(ZI)[1]="Intercept"
has_covariates=TRUE
}

fit_final=greedy_nb_warm_start(
X=X,y=y,Z=Z,L.init=L.init,theta_init=theta_init,
estimate_theta=estimate_theta,cor_method=init_cor_method
)

alpha=clean_coef(coef(fit_final)[seq_len(ncol(ZI))])
g=numeric(0)
beta=rep(0,p)
beta_prev=beta
alpha_prev=alpha*0
fitX=NULL
early_no_cs=FALSE

for (iter in seq_len(max.iter)) {
beta_prev=beta
alpha_prev=alpha

eta=fit_final$linear.predictors
mu =fit_final$fitted.values
dmu_deta=mu

theta=if (!is.null(fit_final$theta)) fit_final$theta else theta_init

var_mu=mu+(mu^2)/theta

g_prime_mu=1/dmu_deta

pseudo_response=eta+(y-mu)*g_prime_mu

W_diag=1/(var_mu*g_prime_mu^2)

bad=!is.finite(pseudo_response) | !is.finite(W_diag) | (W_diag <= 0)
if (mean(bad) > 0.9) {
warning("Too many invalid observations at iteration ",iter)
break
}
if (any(bad)) {
W_diag[bad]=0
pseudo_response[bad]=0
}

W_diag=robust_weight(W_diag,cutoff=weight_cutoff)

weight_denom=sum(W_diag^2)
if (!is.finite(weight_denom) || weight_denom <= 0) {
stop("All working weights are zero at iteration ",iter)
}
n_eff=(sum(W_diag))^2/weight_denom

phi0=summary(fit_final)$dispersion
W_invsqrt=sqrt(W_diag/phi0)

tilde_y=pseudo_response*W_invsqrt
tilde_X=X*W_invsqrt
tilde_Z=ZI*W_invsqrt

ZtZ =matrixMultiply(tilde_Z,tilde_Z,transA=TRUE)+diag(1e-8,ncol(tilde_Z))
ZinvZt=matrixMultiply(solve_with_ridge(ZtZ),tilde_Z,transB=TRUE)
tilde_y=as.numeric(tilde_y-matrixVectorMultiply(tilde_Z,matrixVectorMultiply(ZinvZt,tilde_y)))
tilde_X=tilde_X-ProjectRes(A=tilde_X,B=tilde_Z,n_threads=n_threads)

XtX=blockwise_crossprod(tilde_X,n_threads=n_threads)
Xty=as.numeric(matrixMultiply(tilde_X,matrix(tilde_y,ncol=1),transA=TRUE))
yty=sum(tilde_y^2)

fitX=susie_ss(
XtX=XtX,Xty=Xty,yty=yty,
n=n,
L=L,
scaled_prior_variance=scaled_prior_variance,
estimate_residual_variance=estimate_residual_variance,
residual_variance=residual_variance,
max_iter = susie.iter,
estimate_prior_method="EM",
coverage=coverage,...
)

bcoef=tryCatch(as.numeric(susieR::coef.susie(fitX)),error=function(e) numeric(0))
if (length(bcoef) == p+1) {
beta=clean_coef(bcoef[-1])
} else if (length(bcoef) == p) {
beta=clean_coef(bcoef)
} else {
beta=beta_prev
}

beta.cs=group.pip.filter(pip.summary=summary(fitX)$var,xQTL.cred.thres=coverage,xQTL.pip.thres=pip.thres)
pip.alive=beta.cs$ind.keep
if (length(pip.alive)) beta[-pip.alive]=0

CSdt=summary(fitX)$vars
cs_indices=sort(unique(CSdt$cs[CSdt$cs > 0]))

if (!length(cs_indices)) {
stop("No credible set detected at iteration ", iter)
}

Alpha_filtered=fitX$alpha*0
for (i_cs in cs_indices) {
vars_in_cs_i=CSdt$variable[CSdt$cs == i_cs]
Alpha_filtered[i_cs,vars_in_cs_i]=fitX$alpha[i_cs,vars_in_cs_i]
}
# Align within-CS SNP directions while preserving PIP weights.
Alpha_filtered=Alpha_filtered*sign(fitX$mu)
XCS=matrixMultiply(X,as.matrix(Alpha_filtered),transB=TRUE)
XCS=XCS[,cs_indices,drop=FALSE]
if (is.null(dim(XCS))) XCS=matrix(XCS,ncol=1)
colnames(XCS)=paste0("Main_CS",cs_indices)

if (!has_covariates) {
Data=data.frame(XCS)
formula_str="y ~ ."
} else {
Data=data.frame(Z,XCS)
formula_str="y ~ ."
}

if (estimate_theta) {
fit_final=tryCatch(
MASS::glm.nb(as.formula(formula_str),data=Data,link="log"),
error=function(e) {
if (verbose) warning("glm.nb failed,using fixed theta=",theta)
glm(as.formula(formula_str),data=Data,family=MASS::negative.binomial(theta,link="log"))
}
)
} else {
fit_final=glm(as.formula(formula_str),data=Data,
             family=MASS::negative.binomial(theta,link="log"))
}

alpha=clean_coef(coef(fit_final)[seq_len(ncol(ZI))])

if (has_covariates) {
err=max(
sqrt(mean((beta-beta_prev)^2)),
if (length(alpha)) sqrt(mean((alpha-alpha_prev)^2)) else 0
)
} else {
err=sqrt(mean((beta-beta_prev)^2))
}

g[iter]=err

if (verbose) {
current_theta=ifelse(is.null(fit_final$theta),theta,fit_final$theta)
cat(sprintf("Iteration %d: err=%.3e,n_eff=%.1f,theta=%.4f,mean(w)=%.2f\n",
        iter,err,n_eff,current_theta,mean(W_diag[W_diag > 0])))
}

if (err < max.eps && iter > min.iter) {
if (verbose) cat("Converged!\n")
break
}
}

if (!is.null(fitX)) {
MainIndex=Identifying_MainEffect(fitX,colnames(X))
} else {
MainIndex=data.frame(
Index=integer(0),Variable=character(0),CS=character(0),
log10BF=numeric(0),PIP=numeric(0),Pvalue=numeric(0),
stringsAsFactors=FALSE
)
}

G=tryCatch(summary(fit_final)$coefficients,error=function(e) NULL)
if (!is.null(G)) MainIndex=safe_add_p(MainIndex,G)

if (verbose && length(g)) {
plot(g,type="o",col="black",pch=16,
xlab="Iteration",
ylab="Max Parameter Change",
main="Convergence Trace (Negative Binomial)")
for (i in seq_along(g)) {
graphics::text(x=i,y=g[i],labels=formatC(g[i],format="e",digits=1),pos=3,cex=0.7,col="red")
}
}

last_err=if (length(g)) utils::tail(g,1) else Inf
list(
iter=if (exists("iter")) iter else 0,
error=g,
converged=if (early_no_cs) TRUE else (exists("iter") && iter < max.iter && last_err < max.eps),
fitX=fitX,
fitJoint=fit_final,
main_index=MainIndex,
JointCoef=G,
n_eff=n_eff,
theta=if (!is.null(fit_final$theta)) fit_final$theta else theta_init,
early_no_cs=early_no_cs
)
}
