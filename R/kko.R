#' variable selection for additive model via KKO
#'
#' @description
#' The function applys KKO to compute importance scores of components.
#'
#' @param X design matrix of additive model; rows are observations and columns are variables.
#' @param y response of addtive model.
#' @param X_k knockoffs matrix of design; the same size as X.
#' @param rfn_range  a vector of random feature expansion numbers to be tuned.
#' @param n_stb_tune  number of subsampling for tuning random feature numbers.
#' @param n_stb number of subsampling for computing importance scores.
#' @param cv_folds the folds of cross-validation for tuning group lasso penalty.
#' @param frac_stb fraction of subsample size.
#' @param nCores_para number of cores for parallelizing subsampling.
#' @param rkernel kernel choices. Default is "laplacian". Other choices are "cauchy" and "gaussian".
#' @param rk_scale scale parameter of sampling distribution for random feature expansion. For gaussian kernel, it is standard deviation of gaussian sampling distribution.
#'
#' @return a list of selection results.
#'\tabular{ll}{
#' \code{importance_score}  \tab  importance scores of variables for knockoff filtering.  \cr
#' \code{selection_frequency}  \tab a 0/1 matrix of selection results on subsamples.
#' Rows are subsamples, and columns are variables.
#' The first half columns are variables of design X, and the latter are knockoffs X_k  \cr
#' \code{rfn_tune}  \tab  tuned optimal random feature number.  \cr
#' \code{rfn_range}  \tab range of random feature numbers. \cr
#' \code{tune_result} \tab a list of tuning results.  \cr
#' }
#'
#'
#'
#' @author Xiaowu Dai, Xiang Lyu, Lexin Li
#'
#' @examples
#' library(knockoff)
#' p=4 # number of predictors
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
#' kko(X,y,X_k,rfn_range,n_stb_tune,n_stb,cv_folds,frac_stb,nCores_para,rkernel,rk_scale)
#'
#'
#'
#' @export
#'
#'
#'

kko=function(X,y,X_k,rfn_range=c(2,3,4),n_stb_tune=50,n_stb=100,cv_folds=10,frac_stb=1/2,nCores_para=4,rkernel=c("laplacian","gaussian","cauchy"),rk_scale=1){

  # n=dim(X)[1]
  p=dim(X)[2]

  ## tune random feature number
  tune_result=rk_tune(X,y,X_k,rfn_range,n_stb_tune,cv_folds,
                      frac_stb,nCores_para,rkernel,rk_scale)
  rfn_tune=tune_result[["rfn_tune"]]


  ## selection frequencies given the tuned feature number
  Pi= rk_subsample(X,y,X_k,rfn_tune,n_stb,cv_folds,frac_stb,nCores_para,rkernel,rk_scale)

  ##  combine selection frequencies from tuning
  Pi_all=rbind(Pi,tune_result$Pi_list[[as.character(rfn_tune)]])

  ##  importance scores
  W=apply(Pi_all,2,mean)
  W=W[1:p]-W[(1+p):(2*p)]


  return(list(importance_score=W, selection_frequency=Pi_all,
              rfn_tune=rfn_tune, rfn_range=rfn_range,tune_result=tune_result))
}
