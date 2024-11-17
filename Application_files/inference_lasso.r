inference_lasso <- function(data, obs_N, nFactors = NULL, obs_T, nFolds, foldid, coef_index){
  
  # Pull out the entities needed 
  X_data <- data$x
  Y_data <- data$y
  p <- dim(X_data)[2]
  
  H0_coef_index <- coef_index
  
  #============================================================================#
  # Step 1: Calculate the projection matrix needed for inference
  #============================================================================#
  
  # Calculate the cross-sectional averages of the regressors
  X_bar <- matrix(NA,ncol=p,nrow=obs_T)
  for(t in 1:obs_T){
    indices <- seq(t, obs_N*obs_T, by=obs_T) 
    X_bar[t,] <- colMeans(X_data[indices,])
  }
  
  # Calculate empirical covariance matrix and eigenstructure without the "j-th" regressor 
  Cov_X_bar <- 1/obs_T * t(X_bar[,-H0_coef_index]) %*% X_bar[,-H0_coef_index] 
  Cov_X_bar_eigen <- eigen(Cov_X_bar,symmetric=TRUE)
  
  
  #============================================================================#
  # Step 1(b): Estimation of number of factors
  #============================================================================#
  
  # Normalize the eigenvalues by the largest one
  eigen_values <- Cov_X_bar_eigen$values/Cov_X_bar_eigen$values[1]
  
  if(is.numeric(nFactors)){
    K_hat <- nFactors
  } else{
  ## Count the number of normalized eigenvalues larger than 1% 
    K_hat <- sum((0.01 < eigen_values))
  }
  
  #============================================================================#
  # Step 1(c): Computation of projection matrix
  #============================================================================#
  
  W_tilde <- X_bar[,-H0_coef_index] %*% Cov_X_bar_eigen$vectors[,1:K_hat]
  Pi_tilde <- diag(obs_T) -  W_tilde %*% solve(t(W_tilde) %*% W_tilde)  %*% t(W_tilde)
  
  
  #============================================================================#
  # Step 2: Lasso and nodewise lasso estimation on transformed data
  #============================================================================#
  
  # Get transformed data
  Y_tilde <- rep(NA, obs_T*obs_N)
  X_tilde <- matrix(NA, nrow=obs_N*obs_T, ncol=p)
  for(i in 1:obs_N){
    indices <- ((i-1) * obs_T + 1) : (i * obs_T) 
    Y_tilde[indices] <- Pi_tilde %*% t(t(Y_data[indices]))
    X_tilde[indices,] <- Pi_tilde %*% t(t(X_data[indices,]))
  }
  
  
  #============================================================================#
  # Step 2(a): Run the Lasso estimation on transformed data and get variance estimate
  #============================================================================#
  
  # Get the main equation lasso estimates using Cross-validation 
  library(glmnet)
  
  fit_Lasso.cv <- cv.glmnet(x = X_tilde, y = Y_tilde, nfolds = nFolds, foldid = foldid, standardize = FALSE, family = "gaussian", alpha=1, intercept=FALSE)
  coef_LASSO <- coef(fit_Lasso.cv, s = "lambda.min")[-1]
  Lambda_CV  <- fit_Lasso.cv$lambda.min
  
  # Get the variance estimate using cross-validated penalty
  yhat_Lasso <- predict(fit_Lasso.cv, newx = X_tilde, type = "response", s = "lambda.min")
  resid_Lasso <- Y_tilde - yhat_Lasso
  var_est_Lasso <- obs_T/(obs_T-K_hat) * mean(resid_Lasso^2)
  
  #============================================================================#
  # Step 2(b): Run the nodewise Lasso estimation on transformed data 
  #============================================================================#
  
  # Run the nodewise lasso estimates using CV to get intial estimates
  fit_node_Lasso.cv <- cv.glmnet(x = X_tilde[,-H0_coef_index], y = X_tilde[, H0_coef_index], nfolds = nFolds, foldid = foldid, standardize = FALSE, family="gaussian", alpha = 1, intercept = FALSE)
  
  # Get the kappa grid that cv.glmnet used to perform cross-validation over
  kappa_grid <- fit_node_Lasso.cv$lambda
  
  # Get the index for the CV chosen kappa
  kappa_node_CV_idx <- fit_node_Lasso.cv$index[1]
  
  #============================================================================#
  # Step 3: Calculate the de-sparsified estimator
  #============================================================================#
  
  #============================================================================#
  # Step 3(a): Choose the nodewise penalty parameter
  #============================================================================#
  
  
  # Calculate the "asymptotic variance" of /hat{b}_j//sigma_eps for every value of kappa in the grid
  #-----------------------------------------------------------------------------------------------
  
  # Determine lambda grid length
  kappa_grid_len <- length(kappa_grid)
  
  # Initialize variance for each kappa
  var_scaled <- rep(0, len = kappa_grid_len)
  num <- var_scaled
  denom <- var_scaled
  
  # Loop over the grid and calculate the "scaled asymptotic variance"
  for(i in 1:kappa_grid_len){
    # (i) Get the residuals from the nodewise lasso
    yhat_node_Lasso <- predict(fit_node_Lasso.cv, newx = X_tilde[,-H0_coef_index], type = "response", s = kappa_grid[i])
    resid_node_Lasso <- X_tilde[,H0_coef_index] - yhat_node_Lasso
    
    # (ii) Compute the "scaled asymptotic variance
    num[i] <- t(resid_node_Lasso) %*% resid_node_Lasso
    denom[i] <- (t(X_tilde[,H0_coef_index]) %*% resid_node_Lasso)^2
    var_scaled[i] <-  t(resid_node_Lasso) %*% resid_node_Lasso / (t(X_tilde[,H0_coef_index]) %*% resid_node_Lasso)^2 
    }

  # Extract the lambda grid index of the CV choice
  kappa_node_CV_idx <- fit_node_Lasso.cv$index[1]
  
  # Choice (b): Use the Buehlmann 25% increase in variance above that for CV rule. 
  V_BvG_trunc <- 1.25 * var_scaled[kappa_node_CV_idx]
  # Get the indicator whether the scaled variance is less than the threshold
  BvG_ind <- (var_scaled<V_BvG_trunc)
  # Find the lambda index of the kappa choice and its value
  kappa_node_idx <- sum(BvG_ind)
  kappa_node <- fit_node_Lasso.cv$lambda[kappa_node_idx]
  
  #============================================================================#
  # Step 3(b): Construct the de-sparsified estimator
  #============================================================================#
  
  yhat_node_Lasso <- predict(fit_node_Lasso.cv, newx = X_tilde[,-H0_coef_index], type = "response", s = kappa_grid[kappa_node_idx])
  resid_node_Lasso <- X_tilde[,H0_coef_index] - yhat_node_Lasso

  coef_despar <- coef_LASSO[H0_coef_index] + t(resid_node_Lasso) %*% resid_Lasso / t(resid_node_Lasso) %*% X_tilde[,H0_coef_index]
  
  Avar <- sqrt(var_est_Lasso * var_scaled[kappa_node_idx])

  return(list(coef_despar = coef_despar, Avar = Avar))
  
  # Larger return for diagnostics during coding.
  #return(list(coef_LASSO = coef_LASSO, Lambda_CV = Lambda_CV, K_hat = K_hat, 
  #            resid_Lasso = resid_Lasso, var_scaled = var_scaled, kappa_grid = kappa_grid,
  #            kappa_node_CV_idx=kappa_node_CV_idx, kappa_node_idx = kappa_node_idx,
  #            coef_despar = coef_despar, Avar = Avar)  )
}  


