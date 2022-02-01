#' nonparametric additive model seleciton via random kernel
#'
#' @description
#' The function selects additive components via applying group lasso on random feature expansion of data and knockoffs.
#'
#'
#' @param X design matrix of additive model; rows are observations and columns are variables.
#' @param y response of addtive model.
#' @param X_k knockoffs matrix of design; the same size as X.
#' @param rfn  random feature expansion number.
#' @param cv_folds the folds of cross-validation for tuning group lasso penalty.
#' @param rkernel kernel choices. Default is "laplacian". Other choices are "cauchy" and "gaussian".
#' @param rk_scale scaling parameter of sampling distribution for random feature expansion. For gaussian kernel, it is standard deviation of gaussian sampling distribution.
#' @param rseed  seed for random feature expansion.
#'
#' @return a 0/1 vector indicating selected components.
#'
#' @author Xiaowu Dai, Xiang Lyu, Lexin Li
#'
#' @examples
#' library(knockoff)
#' p=5 # number of predictors
#' sig_mag=100 # signal strength
#' n= 200 # sample size
#' rkernel="laplacian" # kernel choice
#' s=2  # sparsity, number of nonzero component functions
#' rk_scale=1  # scaling paramtere of kernel
#' rfn= 3  # number of random features
#' cv_folds=15  # folds of cross-validation in group lasso
#' X=matrix(rnorm(n*p),n,p)%*%chol(toeplitz(0.3^(0:(p-1))))   # generate design
#' X_k = create.second_order(X) # generate knockoff
#' reg_coef=c(rep(1,s),rep(0,p-s))  # regression coefficient
#' reg_coef=reg_coef*(2*(rnorm(p)>0)-1)*sig_mag
#' y=X%*% reg_coef + rnorm(n) # response
#'
#' # the first half is variables of design X, and the latter is knockoffs X_k
#' rk_fit(X,y,X_k,rfn,cv_folds,rkernel,rk_scale)
#'
#'
#' @export
#'
#' @import grpreg stats
#' @importFrom  ExtDist rLaplace
#'

rk_fit=function(X,y,X_k,rfn,cv_folds,rkernel="laplacian",rk_scale=1,rseed=NULL){

  p=dim(X)[2];n=dim(X)[1]

  if (!is.null(rseed)){ # set random seed for random feature expansion
    set.seed(rseed)
  }

  ## kernel choice and random features
  ## expansion by cos(wx+b)
  if (rkernel=="gaussian"){
    w=rnorm(2*p*rfn,sd=rk_scale)
  } else if (rkernel=="laplacian") {
    w=rcauchy(2*p*rfn,scale=rk_scale)
  } else if (rkernel=="cauchy"){
    w=rLaplace(2*p*rfn,b=rk_scale)
  } else {
    stop("Kernel type is unknown.")
  }
  b=runif(2*p*rfn,0,2*pi)


  X_kernel=NULL  # random feature expansion of data
  for ( i_var in 1:(2*p)){
    if (i_var<=p){ # original data
      X_kernel=cbind(X_kernel,cos(t(t(as.matrix(X[,i_var]) %*% w[((i_var-1)*rfn+1):(i_var*rfn)]) +b[((i_var-1)*rfn+1):(i_var*rfn)])))
    } else { # knockoff
      X_kernel=cbind(X_kernel,cos(t(t(as.matrix(X_k[,i_var-p]) %*% w[((i_var-1)*rfn+1):(i_var*rfn)]) +b[((i_var-1)*rfn+1):(i_var*rfn)])))
    }
  }
  rm(X,X_k);

  ### variable selection by group lasso
  cvfit = cv.grpreg(scale(X_kernel,scale=FALSE), scale(y,scale=FALSE),  # scale X and y to remove intercept term
                    factor(rep(1:(2*p),each=rfn)),penalty="grLasso",nfolds=cv_folds )

  beta_mat=cvfit$fit$beta[-1,] # coefficients
  selected_path=0 # solution path
  for ( i in 1:(2*p)){
    selected_path=selected_path + ((apply(abs(beta_mat[((i-1)*rfn+1):(i*rfn),]),2,mean)!=0) +0)
  }

  # tune by BIC
  fit_coef=cvfit$fit$beta[-1,which.min(log(cvfit$cve) +rfn*log(n)/n*selected_path)]


  ## nonzero coefficient
  fit_nonzero=c()
  for (i_var in 1:(2*p)){
    if ( all(fit_coef[((i_var-1)*rfn+1):(i_var*rfn)]==0)){
      fit_nonzero=c(fit_nonzero,0)
    } else {
      fit_nonzero=c(fit_nonzero,1)
    }
  }

  return((fit_nonzero!=0)+0)
}
