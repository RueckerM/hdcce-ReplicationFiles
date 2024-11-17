#' Estimate HD panels with IFE
#'
#' \code{hd_panel_estimator} fits the high-dimensional CCE estimation procedure
#' proposed in Vogt, M., Walsh, C., and Linton, 0. (2022) "Estimation of
#' High-Dimensional Panel Data Models with Interactive Fixed Effects".
#'
#' @param data List containing the balanced panel data with data$y containing
#'   the dependent variables and data$x the regressors. Both are sorted such
#'   that first T observations are those for unit 1, followed by the T
#'   observations for unit 2, etc.
#' @param obs_N,obs_T The number of cross-section units (obs_N) and (obs_T) in
#'   the data. Must both be provided.
#' @param TRUNC The truncation parameter tau used to estimate the number of
#'   factors. Default is 0.05.
#' @param NFACTORS Allows to set the number of factors use in estimation.
#'   Default is NULL so that the data driven choice with TRUNC=0.05 is used.
#' @param variant Choose whether to use the lasso based estimator
#'   (variant="Lasso") or the least squares variant (variant="LS"). Must be
#'   specified.
#' @param penalty_choice Setting penalty_choice = "None" will ensure that no
#'   data driven lasso penalty selection will be done. Default is to run
#'   cross-validation to select the lasso penalty parameter.
#' @param NFOLDS Specify the number of folds to use when using cross-validation
#'   to select lasso penalty based on cv.glmnet. Default is the default in
#'   cv.glmnet, which is 10.
#' @param NLAMBDA Specify the number of lasso penalties to use in glmnet when
#'   penalty_choice= "None". Default is the default in glmnet, which is 100.
#' @param scree_plot Logical variable to indicate whether a scree plot of the
#'   eigendecomposition of \eqn{\Sigma} should be shown. The default is TRUE.
#' @return In general the function returns the estimated coefficients, the
#' number of estimated factors and the lasso penalty parameter. More
#' specifically:
#' \itemize{
#' \item Unless penalty_choice ="None", the function returns a list with
#' \describe{
#'   \item{$coefs}{The coefficient estimates.}
#'   \item{$K_hat}{The estimated number of factors. If NFACTORS was
#'   set by user then this will be returned.}
#'   \item{$Lambda}{The data driven lasso penalty choice. If variant ="LS",
#'   it is set to 0.}
#'   }
#' \item If the option penalty_choice ="None" was set, then the function
#' returns a list with
#'  \describe{
#'   \item{$coefs}{The coefficient estimates.}
#'   \item{$K_hat}{The estimated number of factors. If NFACTORS was
#'   set by user then this will be returned.}
#'   \item{$Lambda}{The values of the lasso penalty parameters used. If
#'   variant ="LS", it is set to 0.}
#'   }
#' }
#' @examples
#' # Load the data set
#' data("data_example")
#'
#' # Set the dimensions of the data
#' obs_N <- 50
#' obs_T <- 50
#'
#' # Do the estimation
#' estimate_model <- hd_panel_estimator(data=data_example,obs_N=obs_N,
#'                       obs_T=obs_T,TRUNC=0.1,penalty_choice = "None",variant="Lasso")
#'
#' @references  Vogt, M., Walsh, C. and Linton, O.  (2022) "Estimation of
#' High-Dimensional Panel Data Models with Interactive Fixed Effects"




#' @export
hd_panel_estimator_fixed_pen <- function(data, obs_N, obs_T,
                               TRUNC = 0.05, NFACTORS = NULL, variant = "LS", 
                               penalty_choice = NULL, NFOLDS=NULL, 
                               fold_vec = NULL, NLAMBDA=NULL, scree_plot = FALSE,
                               pen = NULL){

  # Function to determine whether an input is a natural number
  is.naturalnumber <- function(x, tol = .Machine$double.eps^0.5)
    x > tol & abs(x - round(x)) < tol

  # Initial Checks
  #-------------------------------------------------------------------------- 

  # Check to see whether a correct variant has been supplied
  if ((variant!="LS") * (variant!="Lasso")==1){
    stop('Variant must be set to "LS" or "Lasso"')
  }


  # Pull out the data
  X_data <- data$x
  Y_data <- data$y
  p <- dim(X_data)[2]

  #============================================================================#
  # Step 1: Eigendecomposition of empirical covariance matrix
  #============================================================================#

  # Calculate the cross-sectional averages of the regressors
  X_bar <- matrix(NA, ncol = p, nrow = obs_T)
  for(t in 1:obs_T){
     indices <- seq(t,obs_N*obs_T,by=obs_T)
     X_bar[t,] <- colMeans(X_data[indices,])
  }

  # Calculate empirical covariance matrix and eigenstructure
  Cov_X_bar <- 1/obs_T * t(X_bar) %*% X_bar
  Cov_X_bar_eigen <- eigen(Cov_X_bar, symmetric=TRUE)



  #============================================================================#
  # Step 2: Estimation of number of factors
  #============================================================================#

  # Normalize the eigenvalues by the largest one
  eigen_values <- Cov_X_bar_eigen$values / (Cov_X_bar_eigen$values[1])

  # Check for user-specified fixed number of factors
  if(class(NFACTORS) == "integer"){
    if(is.naturalnumber(NFACTORS) == TRUE){
      K_hat <- NFACTORS
      message(paste("User-supplied number of factors given by 'NFACTORS' = ", NFACTORS," is used in estimation."))

      if( scree_plot == TRUE){
        graphics::plot(eigen_values, ylim=c(0,1), ylab="Normalized Eigenvalues",main=paste("Number of factors set to ", K_hat))
        graphics::points(eigen_values[1:K_hat], col="red")
      }
    }
  } else {
    # Count the number of normalized eigenvalues larger than 5%
    K_hat <- sum((TRUNC<eigen_values))
    message(paste("Number of factors estimated by truncation at level given by 'TRUNC' = ", TRUNC))

    if( scree_plot == TRUE){
      graphics::plot(eigen_values,ylim=c(0,1),ylab="Normalized Eigenvalues",main=paste("Estimated number of factors = ",K_hat))
      graphics::abline(h=TRUNC,col="red")
      graphics::legend("topright",legend=c("Truncation"),lty=c(1),col=c("red"))
    }
  }



  #============================================================================#
  # Step 3: Computation of projection matrix
  #============================================================================#

  W_hat <- X_bar %*% Cov_X_bar_eigen$vectors[,1:K_hat]
  Pi_hat <- diag(obs_T) -  W_hat %*% solve(t(W_hat) %*% W_hat)  %*% t(W_hat)

  #============================================================================#
  # Step 4: Transform the data
  #============================================================================#

  # Get transformed data
  Y_hat <- rep(NA, obs_T * obs_N)
  X_hat <- matrix(NA, nrow = obs_N * obs_T, ncol = p)
  for(i in 1:obs_N){
    indices <- ((i-1)*obs_T + 1):(i*obs_T)
    Y_hat[indices] <- Pi_hat %*% t(t(Y_data[indices]))
    X_hat[indices,] <- Pi_hat %*% t(t(X_data[indices,]))
  }

  #============================================================================#
  # Step 5: Estimate Lasso or OLS on the transformed data
  #============================================================================#

  # Run the LS variant if it was supplied
  if (variant == "LS"){
    # Get least squares estimates
    coef_est <- lm(Y_hat ~ X_hat - 1)$coefficients
    results <- list(coefs = coef_est, K_hat = K_hat)
  # Close the variant == "LS" block
  }

  # Run the Lasso variant if it was supplied
  if(variant == "Lasso"){
    # Check to see whether to NOT do data driven tuning parameter selection
    if(penalty_choice == "None"){

      # Check for user-specified fixed number of lambdas
      if(class(NLAMBDA) == "numeric"){
        if(is.naturalnumber(NLAMBDA)==TRUE){
          NLAMBDA <- NLAMBDA
          message(paste("User-supplied number of penalties given by 'NLAMBDA' = ",
                        NLAMBDA," is used in call to glmnet."))

          fit_Lasso <- glmnet::glmnet(x=X_hat,y=Y_hat,family="gaussian",
                                      alpha=1,nlambda = NLAMBDA,intercept=FALSE)

          results <- list(coefs=fit_Lasso$beta,K_hat = K_hat,Lambda = fit_Lasso$lambda)

        } else {
          stop(paste("Supplied value of NLAMBDA = ",NLAMBDA, "was not an integer."))
        }
      } else {
          # Run glmnet with default number of Lambda values
          fit_Lasso <- glmnet::glmnet(x=X_hat, y=Y_hat, family="gaussian", 
                                      alpha=1, intercept=FALSE)

          results <- list(coefs=fit_Lasso$beta,K_hat = K_hat, 
                          Lambda = fit_Lasso$lambda)
      }
    # Closes penalty == "None" block
    } 
    
    
    else {
      # Run the cross-validated code as a default

      # Check for user-specified folds
      if(class(NFOLDS) == "numeric" & class(fold_vec) != "numeric"){
        if(is.naturalnumber(NFOLDS) == TRUE){
          NFOLDS <- NFOLDS
          message(paste("User-supplied number of folds given by 'NFOLDS' = ", 
                        NFOLDS," is used in call to cv.glmnet."))

          fit_Lasso.cv <- glmnet::cv.glmnet(x=X_hat, y=Y_hat, nfolds=NFOLDS, 
                                            family="gaussian", alpha=1, 
                                            intercept=FALSE)

          coefs <- stats::coef(fit_Lasso.cv,s="lambda.min")[-1]
          Lambda_CV  <- fit_Lasso.cv$lambda.min
          
          results <- list(coefs=coefs, K_hat = K_hat, Lambda = Lambda_CV)
        } 
        else {
          stop(paste("Supplied value of NFOLDS = ", NFOLDS, "was not an integer."))
        }
      # Check if there is user-specified set of folds
      } 
      else if(class(fold_vec) == "numeric"){
          fit_Lasso.cv <- glmnet::cv.glmnet(x=X_hat, y=Y_hat, nfolds = NFOLDS, 
                                            foldid = fold_vec, family="gaussian", 
                                            alpha=1, intercept=FALSE)
          coefs <- stats::coef(fit_Lasso.cv, s = "lambda.min")[-1]
          Lambda_CV  <- fit_Lasso.cv$lambda.min
          print("User specific Fold Vec selected")
          
          # Also calculate the coefficients for a fixed penalty
          coefs_pen <- stats::coef(fit_Lasso.cv, s = pen)[-1]
          
          results <- list(coefs = coefs, K_hat = K_hat, Lambda = Lambda_CV, 
                          coefs_pen = coefs_pen,
                          eigen_values = eigen_values)
     }
      else{
          # Run cv.glmnet with default number of folds
          fit_Lasso.cv <- glmnet::cv.glmnet(x=X_hat, y=Y_hat, family="gaussian", 
                                            alpha=1, intercept=FALSE)
          coefs <- stats::coef(fit_Lasso.cv,s="lambda.min")[-1]
          Lambda_CV  <- fit_Lasso.cv$lambda.min
          results <- list(coefs = coefs, K_hat = K_hat,Lambda = Lambda_CV)
     }
  }

# Close the variant == "Lasso" block
  }
  # Return the estimation results
  return(results)
# Close the estimation function
}

# Run the function
#tt <- hd_panel_estimator(data=data,obs_N=obs_N,obs_T=obs_T,TRUNC=0.1,penalty_choice = "None",variant="Lasso")


