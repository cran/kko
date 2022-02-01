#' generate response from nonparametric additive model
#'
#' @description
#' The function generate response from additive models of various components.
#'
#'
#' @param X design matrix of additive model; rows are observations and columns are variables.
#' @param reg_coef regression coefficient vector.
#' @param model  types of components. Default is "linear". Other choices are
#' \tabular{ll}{
#' \code{linear}  \tab  linear regression.  \cr
#' \code{poly}  \tab   polynomial of degree sampled from 2 to 4. \cr
#' \code{sinpoly} \tab sum of polynomial of sin and cos. \cr
#' \code{sinratio}  \tab  ratio of sin. \cr
#' \code{sinmix} \tab   sampled from poly and sinratio. \cr
#' }
#' @param err_sd standard deviation of regression error.
#'
#' @return reponse vector
#'
#' @author Xiaowu Dai, Xiang Lyu, Lexin Li
#'
#' @examples
#' p=5 # number of predictors
#' s=2  # sparsity, number of nonzero component functions
#' sig_mag=100 # signal strength
#' n= 200 # sample size
#' model="poly" # component function type
#' X=matrix(rnorm(n*p),n,p) %*%chol(toeplitz(0.3^(0:(p-1))))   # generate design
#' reg_coef=c(rep(1,s),rep(0,p-s))  # regression coefficient
#' reg_coef=reg_coef*(2*(rnorm(p)>0)-1)*sig_mag
#' y=generate_data(X,reg_coef,model) # reponse vector
#'
#'
#' @export
#'
#'


generate_data=function(X,reg_coef,model="linear",err_sd=1){

  n=dim(X)[1]
  if (model=="linear"){

  } else if (model=="poly"){
    # transfer X into polynomial of 2-4 degree
    for ( i_var  in which(reg_coef!=0)) {
      X[,i_var]=X[,i_var]^sample(2:4,1)
    }
  } else  if (model=="sinmix"){
    ### mix of sinpoly and sinratio
    for ( i_var  in which(reg_coef!=0)) {
      if (rnorm(1)>0){ # 1/2 probability
        X[,i_var]=sin(sample(1:10,1)*X[,i_var]) / (2-sin(sample(1:10,1)*X[,i_var]))
      } else {
        X[,i_var]=runif(1,min=1,max=2)*sin(sample(1:10,1)*X[,i_var])
        + runif(1,min=1,max=2)*cos(sample(1:10,1)*X[,i_var])
        + runif(1,min=1,max=2)*sin(sample(1:10,1)*X[,i_var])^2
        + runif(1,min=1,max=2)*cos(sample(1:10,1)*X[,i_var])^2
      }
    }
  } else if (model=="sinratio") {
    for ( i_var  in which(reg_coef!=0)) {
      X[,i_var]=sin(sample(1:10,1)*X[,i_var]) / (2-sin(sample(1:10,1)*X[,i_var]))
    }
  } else if (model=='sinpoly') {
    for ( i_var  in which(reg_coef!=0)) {
      X[,i_var]=runif(1,min=1,max=2)*sin(sample(1:10,1)*X[,i_var])
      + runif(1,min=1,max=2)*cos(sample(1:10,1)*X[,i_var])
      + runif(1,min=1,max=2)*sin(sample(1:10,1)*X[,i_var])^2
      + runif(1,min=1,max=2)*cos(sample(1:10,1)*X[,i_var])^2
    }
  } else {
    stop("Model type is unknown for regression.")
  }

  y= X %*%reg_coef + rnorm(n,sd=err_sd)

  return(c(y))
}
