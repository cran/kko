#' evaluate performance of KKO selection 
#'
#' @description
#' The function computes \{FDP, FPR, TPR\} of selection by knockoff filtering on importance scores of KKO. 
#'
#' @param W importance scores of variables.
#' @param reg_coef true regression coefficient.
#' @param fdr_range FDR control levels of knockoff filter.
#' @param offset 0/1. If 1, knockoff+ filter. Otherwise, knockoff filter.
#'
#'
#'
#' @return {FDP, FPR, TPR} of knockoff filtering at fdr_range.
#'
#' @author Xiaowu Dai, Xiang Lyu, Lexin Li
#'
#' @examples
#' library(knockoff)
#' p=5 # number of predictors
#' sig_mag=100 # signal strength
#' n= 100 # sample size
#' rkernel="laplacian" # kernel choice
#' s=2  # sparsity, number of nonzero component functions
#' rk_scale=1  # scaling paramtere of kernel
#' rfn_range=c(2,3,4)  # number of random features
#' cv_folds=15  # folds of cross-validation in group lasso
#' n_stb=10 # number of subsampling for importance scores
#' n_stb_tune=5 # number of subsampling for tuning random feature number
#' frac_stb=1/2 # fraction of subsample
#' nCores_para=2 # number of cores for parallelization
#' X=matrix(rnorm(n*p),n,p)%*%chol(toeplitz(0.3^(0:(p-1))))   # generate design
#' X_k = create.second_order(X) # generate knockoff
#' reg_coef=c(rep(1,s),rep(0,p-s))  # regression coefficient
#' reg_coef=reg_coef*(2*(rnorm(p)>0)-1)*sig_mag
#' y=X%*% reg_coef + rnorm(n) # response
#'
#' kko_fit=kko(X,y,X_k,rfn_range,n_stb_tune,n_stb,cv_folds,frac_stb,nCores_para,rkernel,rk_scale)
#' W=kko_fit$importance_score
#' fdr_range=c(0.2,0.3,0.4,0.5)
#' KO_evaluation(W,reg_coef,fdr_range,offset=1)
#'
#'
#' @export
#'
#' @import knockoff
#'

KO_evaluation=function(W,reg_coef,fdr_range=0.2,offset=1){

  fdp = function(selected,reg_coef) {sum(reg_coef[selected] == 0) / max(1, length(selected))}
  fpr = function(selected,reg_coef) {sum(reg_coef[selected] == 0) / max(1, sum(reg_coef==0))}
  tpr = function(selected,reg_coef) {sum(reg_coef[selected] != 0) / max(1, sum(reg_coef!=0))}

  perf=c()
  for (i_fdr in fdr_range){
    thres = knockoff.threshold(W, fdr=i_fdr, offset=offset)
    selected = which(W >= thres)
    perf=c(perf,fdp(selected,reg_coef),fpr(selected,reg_coef),tpr(selected,reg_coef))
  }
  perf=t(perf)
  colnames(perf)=paste(rep(c("FDP","FPR","TPR"),length(fdr_range)),"_",rep(fdr_range,each=3),sep="")
  return(perf)
}
