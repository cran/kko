#' compute selection frequency of rk_fit on subsamples
#'
#' @description
#' The function applys rk_fit on subsamples and record selection results.
#'
#' @param X design matrix of additive model; rows are observations and columns are variables.
#' @param y response of addtive model.
#' @param X_k knockoffs matrix of design; the same size as X.
#' @param rfn  random feature expansion number.
#' @param n_stb number of subsampling.
#' @param cv_folds the folds of cross-validation for tuning group lasso.
#' @param frac_stb fraction of subsample size.
#' @param nCores_para number of cores for parallelizing subsampling.
#' @param rkernel kernel choices. Default is "laplacian". Other choices are "cauchy" and "gaussian".
#' @param rk_scale scaling parameter of sampling distribution for random feature expansion. For gaussian kernel, it is standard deviation of gaussian sampling distribution.
#'
#'
#' @return a 0/1 matrix indicating selection results. Rows are subsamples, and columns are variables.
#' The first half columns are variables of design X, and the latter are knockoffs X_k.
#'
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
#' rfn= 3  # number of random features
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
#' rk_subsample(X,y,X_k,rfn,n_stb,cv_folds,frac_stb,nCores_para,rkernel,rk_scale)
#'
#'
#' @export
#'
#'
#' @import doParallel parallel foreach
#'
#'


rk_subsample=function(X,y,X_k,rfn,n_stb,cv_folds,frac_stb=1/2,nCores_para,rkernel="laplacian",rk_scale=1){


  n=dim(X)[1]
  p=dim(X)[2]

  # cl=makeCluster(nCores_para,type='SOCK') # register cluster for parallelization to subsample
  # registerDoSNOW(cl)
  cl = parallel::makeCluster(nCores_para)
  doParallel::registerDoParallel(cl)
  i_stb=NA  # visible binding for global variable i_stb in foreach to avoid package check note
  Pi =  foreach(i_stb = 1:n_stb,.combine="rbind",#.packages=c("grpreg","ExtDist"),
                .export=c("rk_fit")) %dopar% {
    rseed=as.numeric(Sys.time())+(i_stb*10+11)^2
    subsample_idx=sample(1:n,floor(n*frac_stb))
    fit_nonzero=rk_fit(X[subsample_idx,],y[subsample_idx],X_k[subsample_idx,],
                       rfn,cv_folds,rkernel,rk_scale,rseed)
    fit_nonzero
  }
  parallel::stopCluster(cl)

  return(Pi)
}
