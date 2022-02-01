#' tune random feature number for KKO.
#'
#' @description
#' The function applys KKO with different random feature numbers to tune the optimal number.
#'
#'
#' @param X design matrix of additive model; rows are observations and columns are variables.
#' @param y response of addtive model.
#' @param X_k knockoffs matrix of design; the same size as X.
#' @param rfn_range a vector of random feature expansion numbers to be tuned.
#' @param n_stb number of subsampling in KKO.
#' @param cv_folds the folds of cross-validation for tuning group lasso.
#' @param frac_stb fraction of subsample.
#' @param nCores_para number of cores for parallelizing subsampling.
#' @param rkernel kernel choices. Default is "laplacian". Other choices are "cauchy" and "gaussian".
#' @param rk_scale scaling parameter of sampling distribution for random feature expansion. For gaussian kernel, it is standard deviation of gaussian sampling distribution.
#'
#' @return a list of tuning results.
#'\tabular{ll}{
#' \code{rfn_tune}  \tab  tuned optimal random feature number.  \cr
#' \code{rfn_range}  \tab a vector of random feature expansion numbers to be tuned. \cr
#' \code{scores} \tab scores of random feature numbers. rfn_tune has the maximal score. \cr
#' \code{Pi_list} \tab a list of subsample selection results for each random feature number. \cr
#' }
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
#' rfn_range= c(2,3,4)  # number of random features
#' cv_folds=15  # folds of cross-validation in group lasso
#' n_stb=10 # number of subsampling
#' frac_stb=1/2 # fraction of subsample
#' nCores_para=2 # number of cores for parallelization
#' X=matrix(rnorm(n*p),n,p)%*%chol(toeplitz(0.3^(0:(p-1))))   # generate design
#' X_k = create.second_order(X) # generate knockoff
#' reg_coef=c(rep(1,s),rep(0,p-s))  # regression coefficient
#' reg_coef=reg_coef*(2*(rnorm(p)>0)-1)*sig_mag
#' y=X%*% reg_coef + rnorm(n) # response
#'
#' rk_tune(X,y,X_k,rfn_range,n_stb,cv_folds,frac_stb,nCores_para,rkernel,rk_scale)
#'
#'
#' @export
#'
#'
#' @import doParallel parallel foreach
#'
#'



rk_tune=function(X,y,X_k,rfn_range,n_stb,cv_folds,frac_stb=1/2,nCores_para=1,rkernel="laplacian",rk_scale=1){

  score_std=NULL   # standard deviation of selection
  Pi_list=list()   # selection frequency

  n=dim(X)[1]
  p=dim(X)[2]

  # cl=makeCluster(nCores_para,type='SOCK') # register cluster for parallelization to subsample
  # registerDoSNOW(cl)
  cl = parallel::makeCluster(nCores_para)
  doParallel::registerDoParallel(cl)
  for (rfn in rfn_range){
    ptm=proc.time()
    ## stability selection
    i_stb=NA  # visible binding for global variable i_stb in foreach to avoid package check note
    Pi = foreach(i_stb = 1:n_stb,.combine="rbind",#.packages=c("grpreg","ExtDist"),
                 .export=c("rk_fit")) %dopar% {
                   rseed=as.numeric(Sys.time())+(i_stb*10+11)^2
                   subsample_idx=sample(1:n,floor(n*frac_stb))
                   fit_nonzero=rk_fit(X[subsample_idx,],y[subsample_idx],X_k[subsample_idx,],
                                      rfn,cv_folds,rkernel,rk_scale,rseed)
                   ## nonzero coefficients
                   fit_nonzero
                 }
    Pi_list[[as.character(rfn)]]=Pi
    score_std=c(score_std, sd(apply(Pi,2,mean)))
  }
  parallel::stopCluster(cl)

  score_std=score_std*2*p-log(rfn_range)
  rfn_tune=rfn_range[which.max(score_std)]

  return(list(rfn_tune=rfn_tune,rfn_range=rfn_range,scores=score_std,Pi_list=Pi_list))
}
