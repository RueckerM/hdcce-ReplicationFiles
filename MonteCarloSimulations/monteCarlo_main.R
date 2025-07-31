# Monte Carlo Simulations for "Estimation and Inference in High-Dimensional Panel
# Data Models with Interactive Fixed Effects" (2025) by Maximilian Rücker, 
# Michael Vogt, Oliver Linton and Christopher Walsh.

########################### Install dependencies ###############################
# install.packages("mvtnorm")
# install.packages("parallel")
# install.packages("doParallel")
# install.packages("doRNG")
# install.packages("glmnet")
# install.packages("latex2exp")
# install.packages("this.path")

setwd(this.path::here())

library(mvtnorm)
library(parallel)
library(doParallel)
library(doRNG) #For appropriate seeding
library(glmnet)
library(latex2exp)



generate_results <- function(N_Sim, T_dim, HAC, gamma){
  ###################### Set fix parameters ########################
  {
    n_dim <- 50
    p_res_K <- vector(mode = 'list', length = 3)
    
    
    if(T_dim == 15){
      pA_vec <- c(7, 10, 13) # Low Dimensional Setting A 
      pB_vec <- c(31, 151, 301) # Setting B
      pC <- 901 # High-Dimensional Setting C
    }else if (T_dim == 50){
      pA_vec <- c(16, 31, 46) # Low Dimensional Setting A 
      pB_vec <- c(91, 451, 901)# Setting B
      pC <- 3001 # High-Dimensional Setting C
    }
    
    K <- 3
    rho <- 0.25
    c_val <- 0.4
    
    
    # Folds for cross-validation
    n_Folds <- 10
    fold_size <- n_dim/n_Folds
    fold_vec <- rep(1, fold_size * T_dim)
    
    for (folds in 2:n_Folds) {
      tmp <- rep(folds, fold_size * T_dim) 
      fold_vec <- append(fold_vec, tmp)
    }
    
    # Initialize parallelization
    c = detectCores()
    cl = makeCluster(c)
    registerDoParallel(cl)
  }
  
  
  results_Factors = foreach(case = 1:N_Sim, .options.RNG = 123) %dorng% {
    
    source("./hdcce_estimator.R", local = TRUE)
    source("./hdcce_inference.R", local = TRUE)
    
    ###################### Data generation ########################
    for (p in 1:3) {
      # Append data für each i to zero row vector
      XA_tmp <- numeric(pA_vec[p])
      XB_tmp <- numeric(pB_vec[p])
      
      if(p == 3){
        XC_tmp <- numeric(pC)
      }
      
      YA_tmp <-  numeric(1)
      YB_tmp <-  numeric(1)
      
      if(p == 3){
        YC_tmp <-  numeric(1)
      }
      
      YA_INFERENCE_tmp <- numeric(3)
      YB_INFERENCE_tmp <- numeric(3)
      if(p == 3){
        YC_INFERENCE_tmp <- numeric(3)
      }
      
      #Group size
      dA <- (pA_vec[p] - 1)/3 
      dB <- (pB_vec[p] - 1)/3
      if(p == 3){
        dC <- (pC - 1)/3
      }
      # Coefficient vector beta
      # Corresponds to small T case
      if(T_dim == 15){
        if(p == 1){
          beta_A <- c(c_val, rep(rep(c_val,2),3))
        }
        if(p == 2){
          beta_A <- c(c_val, rep(rep(c_val,3),3))
        }
        if(p == 3){
          beta_A <- c(c_val, rep(c(rep(c_val,3),0),3))
        }
      }else{
        betaA_group <- c(rep(c_val, 3), rep(0, dA-3))
        beta_A <- c(c_val, betaA_group, betaA_group, betaA_group)
      }
      
      
      betaB_group <- c(rep(c_val, 3), rep(0, dB-3))
      beta_B <- c(c_val, betaB_group, betaB_group, betaB_group)
      
      if(p == 3){
        betaC_group <- c(rep(c_val, 3), rep(0, dC-3))
        beta_C <- c(c_val, betaC_group, betaC_group, betaC_group)
      }
      
      
      
      
      
      # Generate factor structure
      F_1 <- as.vector(arima.sim(model = list(ar = 0.5), n = T_dim, 
                                 innov = rnorm(n = T_dim, mean = 0, sd = sqrt(0.75))))
      F_2 <- as.vector(arima.sim(model = list(ar = 0.5), n = T_dim, 
                                 innov = rnorm(n = T_dim, mean = 0, sd = sqrt(0.75))))
      F_3 <- as.vector(arima.sim(model = list(ar = 0.5), n = T_dim, 
                                 innov = rnorm(n = T_dim, mean = 0, sd = sqrt(0.75))))
      
      F_matrix <- cbind(F_1, F_2, F_3)
      
      # Oracle projection matrix
      Pi <- diag(1, T_dim, T_dim) - F_matrix %*% solve((t(F_matrix) %*% F_matrix)) %*% t(F_matrix)
      
      
      # Expectation of component vector G_i ----> We spike the first eigenvalue with value gamma
      muG_A <- c(rep(1, K), rep(gamma,dA), rep(1, 2*dA))
      muG_B <- c(rep(1, K), rep(gamma,dB), rep(1, 2*dB))
      
      if(p == 3){
        muG_C <- c(rep(1, K), rep(gamma,dC), rep(1, 2*dC))
      }
      
      
      
      SigmaG_A <- matrix(rho, nrow = (K + pA_vec[p] - 1 ), 
                         ncol = (K + pA_vec[p]-1)) + diag(1 - rho, nrow = (K + pA_vec[p] - 1 ), 
                                                          ncol = (K + pA_vec[p] - 1) )
      SigmaG_B <- matrix(rho, nrow = (K + pB_vec[p] - 1 ), 
                         ncol = (K + pB_vec[p]-1)) +  diag(1 - rho, nrow = ( K + pB_vec[p] - 1 ), 
                                                           ncol = (K + pB_vec[p] - 1))
      if(p == 3){
        SigmaG_C <- matrix(rho, nrow = (K + pC - 1 ),
                           ncol = (K + pC-1)) + diag(1 - rho, nrow = ( K + pC - 1 ),
                                                     ncol = (K + pC - 1))
      }
      
      
      
      #Select for HOMOSCEDASTIC DESIGN
      if(HAC == 1){
        eps <- rnorm(n_dim*T_dim, 0, 1)
        u <- rnorm(n_dim*T_dim, 0, 1)
      }
      
      # HETEROSCEDASTIC AND UNCORRELATED DESIGN
      if(HAC == 2){
        sigma_eps <- runif(n_dim, 0.5, 1.5)
        
        eps <- numeric(n_dim * T_dim)
        u <- rnorm(n_dim * T_dim, 0, 1)
        
        for (i in 1:n_dim) {
          eps[((i-1)*T_dim+1):(i*T_dim)] <- rnorm(T_dim, mean = 0, sd = sqrt(sigma_eps[i]))
        }
      }
      
      
      #Select for HETEROSCEDASTIC AND CORRELATED DESIGN
      if(HAC == 3){
        sigma_eps <- runif(n_dim, 0.5, 1.5)
        sigma_u <- 1
        eps <- numeric(n_dim * T_dim)
        u <- numeric(n_dim * T_dim)
        
        for (i in 1:n_dim) {
          eps[((i-1)*T_dim+1):(i*T_dim)] <- as.vector(arima.sim(model = list(ar = 0.5), n = T_dim,
                                                                innov = rnorm(n = T_dim, mean = 0,
                                                                              sd = sqrt(sigma_eps[i] * (3/4)))))
          u[((i-1)*T_dim+1):(i*T_dim)] <- as.vector(arima.sim(model = list(ar = 0.5), n = T_dim,
                                                              innov = rnorm(n = T_dim, mean = 0,
                                                                            sd = sqrt(sigma_u * (3/4)))))
        }
      }
      
      
      
      # idiosyncratic regressor part across i and t
      Z_A <- mvtnorm::rmvnorm(n = n_dim * T_dim, mean = rep(0, (pA_vec[p] - 1), 
                                                            sd = diag(1, nrow = (pA_vec[p] - 1), 
                                                                      ncol = (pA_vec[p] - 1) ))) 
      Z_B <- mvtnorm::rmvnorm(n = n_dim * T_dim, mean = rep(0, (pB_vec[p] - 1), 
                                                            sd = diag(1, nrow = (pB_vec[p] - 1), 
                                                                      ncol = (pB_vec[p] - 1) )))
      if(p ==3){
        Z_C <- mvtnorm::rmvnorm(n = n_dim * T_dim, mean = rep(0, (pC - 1), 
                                                              sd = diag(1, nrow = (pC - 1), 
                                                                        ncol = (pC - 1) )))
      }
      
      G_Ai <- mvtnorm::rmvnorm(n_dim , mean = muG_A, sigma = SigmaG_A) # Model components
      G_Bi <- mvtnorm::rmvnorm(n_dim, mean = muG_B, sigma = SigmaG_B)
      if(p == 3){
        G_Ci <- mvtnorm::rmvnorm(n_dim, mean = muG_C, sigma = SigmaG_C)
      }
      # Empty lists to append the data for each i
      
      for (i in 1:n_dim) {################# Start iteration i #################
        
        gamma_Ai <- G_Ai[i,1:K]
        gamma_Bi <- G_Bi[i,1:K]
        
        if(p == 3){
          gamma_Ci <- G_Ci[i,1:K]
        }
        
        GammaA_i <- cbind(c(G_Ai[i, (K + 1):(K + dA)], rep(0, 2 * dA)),
                          c(rep(0, dA), G_Ai[i, (K + dA + 1):(K + 2 * dA)], rep(0, dA)),
                          c(rep(0, 2 * dA), G_Ai[i, (K + 2 * dA + 1):(K + 3 * dA)]))
        
        GammaB_i <- cbind(c(G_Bi[i, (K + 1):(K + dB)], rep(0, 2 * dB)),
                          c(rep(0, dB), G_Bi[i, (K + dB + 1):(K + 2 * dB)], rep(0, dB)),
                          c(rep(0, 2 * dB), G_Bi[i, (K + 2 * dB + 1):(K + 3 * dB)]))
        
        if(p == 3){
          GammaC_i <- cbind(c(G_Ci[i, (K + 1):(K + dC)], rep(0, 2 * dC)),
                            c(rep(0, dC), G_Ci[i, (K + dC + 1):(K + 2*dC)], rep(0, dC)),
                            c(rep(0, 2 * dC), G_Ci[i, (K + 2*dC + 1):(K + 3*dC)]))
        }
        
        
        # Regressors generated with factor structure
        XA_i <- F_matrix %*% t(GammaA_i) + Z_A[((i-1) * T_dim + 1) : (i * T_dim),]
        XB_i <- F_matrix %*% t(GammaB_i) + Z_B[((i-1) * T_dim + 1) : (i * T_dim),]
        if(p == 3){
          XC_i <- F_matrix %*% t(GammaC_i) + Z_C[((i-1) * T_dim + 1) : (i * T_dim),]
        }
        
        
        # First regressor generated with nodewise structure
        XA_i1 <- XA_i[, 1] * sqrt(2/3) + u[((i-1) * T_dim + 1) : (i * T_dim)]
        XB_i1 <- XB_i[, 1] * sqrt(2/3) + u[((i-1) * T_dim + 1) : (i * T_dim)]
        if(p == 3){
          XC_i1 <- XC_i[, 1] * sqrt(2/3) + u[((i-1) * T_dim + 1) : (i * T_dim)]
        }
        
        # Collect the design X_i in a first step
        XA_i <- cbind(XA_i1, XA_i)
        XB_i <- cbind(XB_i1, XB_i)
        if(p == 3){
          XC_i <- cbind(XC_i1, XC_i)
        }
        
        # Generate response variable Y  
        YA_i <- XA_i %*% beta_A + F_matrix %*% gamma_Ai + eps[((i-1) * T_dim + 1) : (i * T_dim)]
        
        YB_i <- XB_i %*% beta_B + F_matrix %*% gamma_Bi + eps[((i-1) * T_dim + 1) : (i * T_dim)]
        if(p == 3){
          YC_i <- XC_i %*% beta_C + F_matrix %*% gamma_Ci + eps[((i-1) * T_dim + 1) : (i * T_dim)]
        }
        # Generate response Y for INFERENCE
        YA_INFERENCE_i <- matrix(nrow = T_dim, ncol = 3)
        YB_INFERENCE_i <- matrix(nrow = T_dim, ncol = 3)
        if(p == 3){
          YC_INFERENCE_i <- matrix(nrow = T_dim, ncol = 3)
        }
        k <- 1
        for (signal in c(0, 0.1, 0.2)) {
          beta_A_INFERENCE <- beta_A
          beta_A_INFERENCE[1] <- signal
          YA_INFERENCE_i[,k] <- XA_i %*% beta_A_INFERENCE + F_matrix %*% gamma_Ai + eps[((i-1) * T_dim + 1) : (i * T_dim)]
          
          beta_B_INFERENCE  <- beta_B
          beta_B_INFERENCE[1] <- signal
          YB_INFERENCE_i[,k] <- XB_i %*% beta_B_INFERENCE + F_matrix %*% gamma_Bi + eps[((i-1) * T_dim + 1) : (i * T_dim)]
          if(p == 3){
            beta_C_INFERENCE <- beta_C
            beta_C_INFERENCE[1] <- signal
            YC_INFERENCE_i[,k] <- XC_i %*% beta_C_INFERENCE + F_matrix %*% gamma_Ci + eps[((i-1) * T_dim + 1) : (i * T_dim)]
          }
          k <- k + 1
        }
        
        
        # Collect the data 
        XA_tmp <- rbind(XA_tmp, XA_i)
        XB_tmp <- rbind(XB_tmp, XB_i)
        if(p == 3){
          XC_tmp <- rbind(XC_tmp, XC_i)
        }
        
        YA_tmp <- rbind(YA_tmp, YA_i)
        YB_tmp <- rbind(YB_tmp, YB_i)
        if(p == 3){
          YC_tmp <- rbind(YC_tmp, YC_i)
        }
        
        YA_INFERENCE_tmp <- rbind(YA_INFERENCE_tmp, YA_INFERENCE_i)
        YB_INFERENCE_tmp <- rbind(YB_INFERENCE_tmp, YB_INFERENCE_i)
        if(p == 3){
          YC_INFERENCE_tmp <- rbind(YC_INFERENCE_tmp, YC_INFERENCE_i)
        }
        
        
      }
      ################# End iteration i #################
      XA_tmp <- XA_tmp[-1,] # Remove first row
      XB_tmp <- XB_tmp[-1,]
      if(p == 3){
        XC_tmp <- XC_tmp[-1,]
      }
      
      YA_tmp <- YA_tmp[-1]
      YB_tmp <- YB_tmp[-1]
      if(p == 3){
        YC_tmp <- YC_tmp[-1]
      }
      
      YA_INFERENCE_tmp <- YA_INFERENCE_tmp[-1,]
      YB_INFERENCE_tmp <- YB_INFERENCE_tmp[-1,]
      if(p == 3){
        YC_INFERENCE_tmp <- YC_INFERENCE_tmp[-1,]
      }
      
      
      
      # Oracle Projected Data
      YA_oracle <- rep(NA, T_dim * n_dim)
      XA_oracle <- matrix(NA, nrow = T_dim * n_dim, ncol = pA_vec[p])
      
      YB_oracle <- rep(NA, T_dim * n_dim)
      XB_oracle <- matrix(NA, nrow = T_dim * n_dim, ncol = pB_vec[p])
      
      if(p == 3){
        YC_oracle <- rep(NA, T_dim * n_dim)
        XC_oracle <- matrix(NA, nrow = T_dim * n_dim, ncol = pC)
      }
      
      PIeps <- rep(NA, T_dim * n_dim)
      
      
      for(l in 1:n_dim){
        indices <- ((l-1) * T_dim + 1):(l * T_dim)
        YA_oracle[indices] <- Pi %*% t(t(YA_tmp[indices]))
        XA_oracle[indices,] <- Pi %*% t(t(XA_tmp[indices,]))
        
        YB_oracle[indices] <- Pi %*% t(t(YB_tmp[indices]))
        XB_oracle[indices,] <- Pi %*% t(t(XB_tmp[indices,]))
        if(p == 3){
          YC_oracle[indices] <- Pi %*% t(t(YC_tmp[indices]))
          XC_oracle[indices,] <- Pi %*% t(t(XC_tmp[indices,]))
        }
      }
      
      # ############### Collect estimation results ##############
      # ## Lasso results
      # 1.  Estimated Pi - Lasso results
      resALasso <- hdcce_estimator(data = list(x = XA_tmp, y = YA_tmp),
                                          obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01,
                                          variant = "Lasso", NFOLDS = n_Folds,
                                          foldid = fold_vec, scree_plot = FALSE)
      resBLasso <- hdcce_estimator(data = list(x = XB_tmp, y = YB_tmp),
                                          obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01,
                                          variant = "Lasso", NFOLDS = n_Folds,
                                          foldid = fold_vec, scree_plot = FALSE)
      if(p == 3){
        resCLasso <- hdcce_estimator(data = list(x = XC_tmp, y = YC_tmp),
                                            obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01,
                                            variant = "Lasso", NFOLDS = n_Folds,
                                            foldid = fold_vec, scree_plot = FALSE)
      }
      
      # 2. Oracle Lasso results
      tmp_fit_Lasso.cv <- glmnet::cv.glmnet(x = XA_oracle, y = YA_oracle, nfolds = n_Folds,
                                            foldid = fold_vec, family = "gaussian",
                                            alpha = 1, intercept = FALSE)
      resAOracleLasso <- stats::coef(tmp_fit_Lasso.cv, s = "lambda.min")[-1]
      
      tmp_fit_Lasso.cv <- glmnet::cv.glmnet(x = XB_oracle, y = YB_oracle, nfolds = n_Folds,
                                            foldid = fold_vec, family = "gaussian",
                                            alpha = 1, intercept = FALSE)
      resBOracleLasso <- stats::coef(tmp_fit_Lasso.cv, s = "lambda.min")[-1]
      
      if(p == 3){
        tmp_fit_Lasso.cv <- glmnet::cv.glmnet(x = XC_oracle, y = YC_oracle, nfolds = n_Folds,
                                              foldid = fold_vec, family = "gaussian",
                                              alpha = 1, intercept = FALSE)
        resCOracleLasso <- stats::coef(tmp_fit_Lasso.cv, s = "lambda.min")[-1]
      }
      
      
      # ## Least squares based result
      # 3. Oracle least squares
      resAOracleLS <- lm(YA_oracle ~ XA_oracle - 1)$coefficients
      resBOracleLS <- lm(YB_oracle ~ XB_oracle - 1)$coefficients
      
      # 4. Pi_hat projected least squares
      resALS <- hdcce_estimator(data = list(x = XA_tmp, y = YA_tmp),
                                       obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01,
                                       variant = "LS", NFOLDS = n_Folds,
                                       foldid = fold_vec, scree_plot = FALSE)$coefs
      
      resBLS <- hdcce_estimator(data = list(x = XB_tmp, y = YB_tmp),
                                       obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01,
                                       variant = "LS", NFOLDS = n_Folds,
                                       foldid = fold_vec, scree_plot = FALSE)$coefs
      # 5. Oracle Oracle LS
      S_A <- as.logical(beta_A)
      resAOracleOracleLS <- lm(YA_oracle ~ XA_oracle[ ,S_A] - 1)$coefficients
      
      S_B <- as.logical(beta_B)
      resBOracleOracleLS <- lm(YB_oracle ~ XB_oracle[ ,S_B] - 1)$coefficients
      
      if(p == 3){
        S_C <- as.logical(beta_C)
        resCOracleOracleLS <- lm(YC_oracle ~ XC_oracle[ ,S_C] - 1)$coefficients
      }
      
      # 6. Pesarans CCE Solution (Only applicable in scenario A)
      XA_bar <- matrix(NA, ncol = pA_vec[p], nrow = T_dim)
      YA_bar <- numeric(T_dim)
      
      for(t in 1:T_dim){
        indices <- seq(t, n_dim*T_dim, by = T_dim)
        XA_bar[t,] <- colMeans(XA_tmp[indices,])
        YA_bar[t] <- mean(YA_tmp[indices])
      }
      W_A <- cbind(XA_bar, YA_bar)
      PiA_pesaran <- diag(1, nrow = T_dim, ncol = T_dim) - W_A %*% solve(t(W_A) %*% W_A) %*% t(W_A)
      
      # Pesaran Projected Data
      YA_pesaran <- rep(NA, T_dim * n_dim)
      XA_pesaran<- matrix(NA, nrow = T_dim * n_dim, ncol = pA_vec[p])
      
      for(l in 1:n_dim){
        indices <- ((l-1) * T_dim + 1):(l * T_dim)
        YA_pesaran[indices] <- PiA_pesaran %*% t(t(YA_tmp[indices]))
        XA_pesaran[indices,] <- PiA_pesaran %*% t(t(XA_tmp[indices,]))
      }
      resAPesaran <- lm(YA_pesaran ~ XA_pesaran - 1)$coefficients
      
      
      
      
      # 7. Collect results for High K, i.e. overestimate its value
      if(p   >  1){
        HIGHresALS <- hdcce_estimator(data = list(x = XA_tmp, y = YA_tmp),
                                             obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01, NFACTORS = 6,
                                             variant = "LS", NFOLDS = n_Folds,
                                             foldid = fold_vec, scree_plot = FALSE)
        HIGHresBLS <- hdcce_estimator(data = list(x = XB_tmp, y = YB_tmp),
                                             obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01, NFACTORS = 6,
                                             variant = "LS", NFOLDS = n_Folds,
                                             foldid = fold_vec, scree_plot = FALSE)
        
        
        HIGHresALasso <- hdcce_estimator(data = list(x = XA_tmp, y = YA_tmp),
                                                obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01, NFACTORS = 6,
                                                variant = "Lasso", NFOLDS = n_Folds,
                                                foldid = fold_vec, scree_plot = FALSE)
        
        HIGHresBLasso <- hdcce_estimator(data = list(x = XB_tmp, y = YB_tmp),
                                                obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01, NFACTORS = 6,
                                                variant = "Lasso", NFOLDS = n_Folds,
                                                foldid = fold_vec, scree_plot = FALSE)
        if(p == 3){
          HIGHresCLasso <- hdcce_estimator(data = list(x = XC_tmp, y = YC_tmp),
                                                  obs_N = n_dim, obs_T = T_dim, TRUNC = 0.01, NFACTORS = 6,
                                                  variant = "Lasso", NFOLDS = n_Folds,
                                                  foldid = fold_vec, scree_plot = FALSE)
        }else{
          HIGHresCLasso <- NULL
          HIGHresCLS <- NULL
        }
      }else{
        HIGHresALasso <- NULL
        HIGHresALS <- NULL
        
        HIGHresBLasso <- NULL
        HIGHresBLS <- NULL
        
        HIGHresCLasso <- NULL
        HIGHresCLS <- NULL
      }
      
      
      # 8. Inference results
      statisticA <- numeric(3)
      statisticB <- numeric(3)
      if(p == 3){
        statisticC <- numeric(3)
      }
      
      
      for (k in 1:3) {
        resAInference <- hdcce_inference(data = list(x = XA_tmp, y = YA_INFERENCE_tmp[,k]),
                                                obs_N = n_dim, obs_T = T_dim,  foldid = fold_vec,
                                                COEF_INDEX_VEC = 1, HAC = HAC)
        
        resBInference <- hdcce_inference(data = list(x = XB_tmp, y = YB_INFERENCE_tmp[,k]),
                                                obs_N = n_dim, obs_T = T_dim,  foldid = fold_vec,
                                                COEF_INDEX_VEC = 1, HAC = HAC)
        
        if(p == 3){
          resCInference <- hdcce_inference(data = list(x = XC_tmp, y = YC_INFERENCE_tmp[,k]),
                                                  obs_N = n_dim, obs_T = T_dim, foldid = fold_vec, 
                                                  COEF_INDEX_VEC  = 1, HAC = HAC)
        }
        
        statisticA[k] <- (1/resAInference$Avar)*resAInference$coef_despar
        statisticB[k] <- (1/resBInference$Avar)*resBInference$coef_despar
        if(p == 3){
          statisticC[k] <- (1/resCInference$Avar)*resCInference$coef_despar
        }
      }
      
      
      
      
      
      ### Collect all results for A, B and C ###
      results_A <- list(TestStatistic = statisticA, Lasso = resALasso, OracleLasso = resAOracleLasso, LS =
                          resALS, OracleLS = resAOracleLS, OracleOracleLS = resAOracleOracleLS,
                        PesaranCCE = resAPesaran, HighKLasso = HIGHresALasso, HighKLS = HIGHresALS)
      
      results_B <- list(TestStatistic = statisticB, Lasso = resBLasso, OracleLasso = resBOracleLasso, LS = resBLS,
                        OracleLS = resBOracleLS, OracleOracleLS = resBOracleOracleLS,
                        HighKLasso = HIGHresBLasso, HighKLS = HIGHresBLS)
      
      if(p == 3){
        results_C <- list(TestStatistic = statisticC, Lasso = resCLasso, OracleLasso = resCOracleLasso,
                          OracleOracleLS = resCOracleOracleLS, HighKLasso = HIGHresCLasso)
      }else{
        results_C <- NULL
      }
      
      
      p_res_K[[p]] <- list(results_A, results_B, results_C)
      
      
    }
    
    p_res_K
  }
  
  stopCluster(cl)
  ############################### End data generation ############################
  
  
  
  
  
  ################# 1. Accuracy Plots: Figures S.1,S.2 and S.3 ######################
  {if(HAC == 1 & gamma == 1){
    LassoA <-  vector(mode = 'list', length = 3)
    OracleLassoA <-  vector(mode = 'list', length = 3)
    LSA <- vector(mode = 'list', length = 3)
    OracleLSA <- vector(mode = 'list', length = 3)
    OracleOracleLSA <-  vector(mode = 'list', length = 3)
    PesaranCCE <-  vector(mode = 'list', length = 3)
    
    LassoB <-  vector(mode = 'list', length = 3)
    OracleLassoB <- vector(mode = 'list', length = 3)
    LSB <-  vector(mode = 'list', length = 3)
    OracleLSB <-  vector(mode = 'list', length = 3)
    OracleOracleLSB <-  vector(mode = 'list', length = 3)
    
    
    LassoC <-  vector(mode = 'list', length = 1)
    OracleLassoC <- vector(mode = 'list', length = 1)
    OracleOracleLSC <-  vector(mode = 'list', length = 1)
    
    
    
    for (p in 1:3) {
      tmpLassoA <- numeric(pA_vec[p])
      tmpOracleLassoA <- numeric(pA_vec[p])
      tmpLSA <- numeric(pA_vec[p])
      tmpOracleLSA <- numeric(pA_vec[p])
      if(T_dim == 15){
        if(p == 1){
          tmpOracleOracleLSA<-numeric(7)
        }else{
          tmpOracleOracleLSA<-numeric(10)
        }
      }else{
        tmpOracleOracleLSA<-numeric(10)
      }
      
      
      tmpPesaranCCE <- numeric(pA_vec[p])
      
      tmpLassoB <- numeric(pB_vec[p])
      tmpOracleLassoB <- numeric(pB_vec[p])
      tmpLSB <- numeric(pB_vec[p])
      tmpOracleLSB<- numeric(pB_vec[p])
      tmpOracleOracleLSB <- numeric(10)
      
      if(p == 3){
        tmpLassoC <- numeric(pC)
        tmpOracleLassoC <- numeric(pC)
        tmpOracleOracleLSC <- numeric(10)
      }
      
      for (sim in 1:N_Sim) {
        tmpLassoA <- rbind(tmpLassoA, results_Factors[[sim]][[p]][[1]][["Lasso"]]$coefs)
        tmpOracleLassoA <- rbind(tmpOracleLassoA, results_Factors[[sim]][[p]][[1]][["OracleLasso"]])
        tmpLSA<-rbind(tmpLSA, results_Factors[[sim]][[p]][[1]]$LS)
        tmpOracleLSA <- rbind(tmpOracleLSA, results_Factors[[sim]][[p]][[1]][["OracleLS"]])
        tmpOracleOracleLSA<- rbind(tmpOracleOracleLSA, results_Factors[[sim]][[p]][[1]][["OracleOracleLS"]])
        tmpPesaranCCE <- rbind(tmpPesaranCCE, results_Factors[[sim]][[p]][[1]][["PesaranCCE"]])
        
        tmpLassoB <- rbind(tmpLassoB, results_Factors[[sim]][[p]][[2]][["Lasso"]]$coefs)
        tmpOracleLassoB <- rbind(tmpOracleLassoB, results_Factors[[sim]][[p]][[2]][["OracleLasso"]])
        tmpLSB <- rbind(tmpLSB, results_Factors[[sim]][[p]][[2]]$LS)
        tmpOracleLSB<- rbind(tmpOracleLSB,  results_Factors[[sim]][[p]][[2]][["OracleLS"]])
        tmpOracleOracleLSB <- rbind(tmpOracleOracleLSB, results_Factors[[sim]][[p]][[2]][["OracleOracleLS"]])
        
        if(p == 3){
          tmpLassoC <- rbind(tmpLassoC, results_Factors[[sim]][[p]][[3]][["Lasso"]]$coefs)
          tmpOracleLassoC <- rbind(tmpOracleLassoC,results_Factors[[sim]][[p]][[3]][["OracleLasso"]])
          tmpOracleOracleLSC <- rbind(tmpOracleOracleLSC, results_Factors[[sim]][[p]][[3]][["OracleOracleLS"]])
        }
      }
      LassoA[[p]] <-  tmpLassoA[-1,]
      OracleLassoA[[p]] <- tmpOracleLassoA[-1,]
      LSA[[p]] <- tmpLSA[-1,]
      OracleLSA[[p]] <- tmpOracleLSA[-1,]
      OracleOracleLSA[[p]] <- tmpOracleOracleLSA[-1,]
      PesaranCCE[[p]] <- tmpPesaranCCE[-1,]
      
      LassoB[[p]] <- tmpLassoB[-1,]
      OracleLassoB[[p]] <- tmpOracleLassoB[-1,]
      LSB[[p]] <- tmpLSB[-1,]
      OracleLSB[[p]] <- tmpOracleLSB[-1,]
      OracleOracleLSB[[p]] <- tmpOracleOracleLSB[-1,]
      if(p == 3){
        LassoC[[1]] <- tmpLassoC[-1,]
        OracleLassoC[[1]] <- tmpOracleLassoC[-1,]
        OracleOracleLSC[[1]] <- tmpOracleOracleLSC[-1,]
      }
    }
    
    
    ########################### 1.1 Boxplot for T_dim = 15 ######################
    
    library(latex2exp)
    if(T_dim == 15){
      ## Scenario A
      # Boxplot for group representatives for each estimator
      col_vec2 <- c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 4))
      col_vec_A <-  c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 5))
      
      col_vec_ASMALLT <- c("green", "blue", "red", "magenta", "bisque")
      
      
      pdf(file = paste0("./plots/FigS1/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pA_vec[1],".pdf", sep = ""), width = 6, height = 6)
      
      boxplot(data.frame(OracleOracleLSA[[1]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[1]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[1]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[1]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSA[[1]][,6]-rep(c_val, N_Sim), PesaranCCE[[1]][,1]-rep(c_val, N_Sim), 
                         PesaranCCE[[1]][,2]-rep(c_val, N_Sim),  PesaranCCE[[1]][,3]-rep(c_val, N_Sim), PesaranCCE[[1]][,4]-rep(c_val, N_Sim), PesaranCCE[[1]][,6]-rep(c_val, N_Sim), 
                         OracleLSA[[1]][,1]-rep(c_val, N_Sim), OracleLSA[[1]][,2]-rep(c_val, N_Sim),  OracleLSA[[1]][,3]-rep(c_val, N_Sim), OracleLSA[[1]][,4]-rep(c_val, N_Sim), OracleLSA[[1]][,6]-rep(c_val, N_Sim), 
                         LSA[[1]][,1]-rep(c_val, N_Sim), LSA[[1]][,2]-rep(c_val, N_Sim), LSA[[1]][,3]-rep(c_val, N_Sim), LSA[[1]][,4]-rep(c_val, N_Sim), LSA[[1]][,6]-rep(c_val, N_Sim),
                         OracleLassoA[[1]][,1]-rep(c_val, N_Sim), OracleLassoA[[1]][,2]-rep(c_val, N_Sim), OracleLassoA[[1]][,3]-rep(c_val, N_Sim), OracleLassoA[[1]][,4]-rep(c_val, N_Sim), OracleLassoA[[1]][,6]-rep(c_val, N_Sim), LassoA[[1]][,1]-rep(c_val, N_Sim), LassoA[[1]][,2]-rep(c_val, N_Sim), LassoA[[1]][,3]-rep(c_val, N_Sim), LassoA[[1]][,4]-rep(c_val, N_Sim), LassoA[[1]][,6]-rep(c_val, N_Sim)), yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,10), seq(12,16), seq(18,22), seq(24,28), seq(30, 34)),
              col = col_vec_ASMALLT, cex = 0.2, boxwex = 0.5, lwd = 0.8, whisklty = 1)
      title("p =  7",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 8, 15, 21, 26, 32), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 11, lwd = 1)
      abline(h = 17, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 29, lwd = 1)
      abline(v = 0)
      
      dev.off()
      
      
      
      ## Scenario B
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS2/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pB_vec[1],".pdf", sep = ""), width = 6, height = 6)
      
      boxplot(data.frame(OracleOracleLSB[[1]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[1]][,2]-rep(c_val, N_Sim), OracleOracleLSB[[1]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[1]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSB[[1]][,6]-rep(c_val, N_Sim), OracleLSB[[1]][,1]-rep(c_val, N_Sim), OracleLSB[[1]][,2]-rep(c_val, N_Sim),  OracleLSB[[1]][,3]-rep(c_val, N_Sim), OracleLSB[[1]][,5], OracleLSB[[1]][,(2+10)]-rep(c_val, N_Sim), 
                         OracleLSB[[1]][,(5+10)], OracleLSB[[1]][,(2+2*10)]-rep(c_val, N_Sim),OracleLSB[[1]][,(5+2*10)], LSB[[1]][,1]-rep(c_val, N_Sim), LSB[[1]][,2]-rep(c_val, N_Sim), LSB[[1]][,3]-rep(c_val, N_Sim), LSB[[1]][,5], LSB[[1]][,(2+10)]-rep(c_val, N_Sim), LSB[[1]][,(5+10)], LSB[[1]][,(2+2*10)]-rep(c_val, N_Sim), 
                         LSB[[1]][,(5+2*10)],
                         OracleLassoB[[1]][,1]-rep(c_val, N_Sim), OracleLassoB[[1]][,2]-rep(c_val, N_Sim), OracleLassoB[[1]][,3]-rep(c_val, N_Sim), OracleLassoB[[1]][,5], OracleLassoB[[1]][,(2+10)]-rep(c_val, N_Sim), 
                         OracleLassoB[[1]][,(5+10)], OracleLassoB[[1]][,(2+2*10)]-rep(c_val, N_Sim), OracleLassoB[[1]][,(5+2*10)], 
                         LassoB[[1]][,1]-rep(c_val, N_Sim), LassoB[[1]][,2]-rep(c_val, N_Sim), LassoB[[1]][,3]-rep(c_val, N_Sim), LassoB[[1]][,5], LassoB[[1]][,(2+10)]-rep(c_val, N_Sim), LassoB[[1]][,(5+10)], 
                         LassoB[[1]][,(2+2*10)]-rep(c_val, N_Sim), LassoB[[1]][,(5+2*10)]),
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title("p =  31",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
      
      ## Scenario C
      pdf(file = paste0("./plots/FigS3/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pC,".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame(OracleOracleLSC[[1]][,1]-rep(c_val, N_Sim), OracleOracleLSC[[1]][,2]-rep(c_val, N_Sim), OracleOracleLSC[[1]][,3]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[1]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[1]][,6]-rep(c_val, N_Sim), OracleLassoC[[1]][,1]-rep(c_val, N_Sim), OracleLassoC[[1]][,2]-rep(c_val, N_Sim), OracleLassoC[[1]][,3]-rep(c_val, N_Sim), OracleLassoC[[1]][,5], 
                         OracleLassoC[[1]][,(2+300)]-rep(c_val, N_Sim), 
                         OracleLassoC[[1]][,(5+300)], OracleLassoC[[1]][,(2+2*300)]-rep(c_val, N_Sim), OracleLassoC[[1]][,(5+2*300)], 
                         LassoC[[1]][,1]-rep(c_val, N_Sim), LassoC[[1]][,2]-rep(c_val, N_Sim), LassoC[[1]][,3]-rep(c_val, N_Sim), LassoC[[1]][,5], LassoC[[1]][,(2+300)]-rep(c_val, N_Sim), 
                         LassoC[[1]][,(5+300)], LassoC[[1]][,(2+2*300)]-rep(c_val, N_Sim), LassoC[[1]][,(5+2*300)]),
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22)), col = col_vec2, cex = 0.2, boxwex = 0.4, lwd = 0.8, whisklty = 1)
      title("p =  901",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19), labels = c(TeX(r"(LS-O$^2$)"), "L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      # p = 2
      
      ## Scenario A
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS1/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pA_vec[2],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame(OracleOracleLSA[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[2]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[2]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSA[[2]][,6]-rep(c_val, N_Sim),  
                         PesaranCCE[[2]][,1]-rep(c_val, N_Sim), PesaranCCE[[2]][,2]-rep(c_val, N_Sim),  PesaranCCE[[2]][,3]-rep(c_val, N_Sim), PesaranCCE[[2]][,5]-rep(c_val, N_Sim), PesaranCCE[[2]][,8]-rep(c_val, N_Sim),
                         OracleLSA[[2]][,1]-rep(c_val, N_Sim), OracleLSA[[2]][,2]-rep(c_val, N_Sim),  OracleLSA[[2]][,3]-rep(c_val, N_Sim), OracleLSA[[2]][,5]-rep(c_val, N_Sim), OracleLSA[[2]][,8]-rep(c_val, N_Sim),
                         LSA[[2]][,1]-rep(c_val, N_Sim), LSA[[2]][,2]-rep(c_val, N_Sim), LSA[[2]][,3]-rep(c_val, N_Sim), LSA[[2]][,5]-rep(c_val, N_Sim), LSA[[2]][,8]-rep(c_val, N_Sim),
                         OracleLassoA[[2]][,1]-rep(c_val, N_Sim), OracleLassoA[[2]][,2]-rep(c_val, N_Sim), OracleLassoA[[2]][,3]-rep(c_val, N_Sim), OracleLassoA[[2]][,5]-rep(c_val, N_Sim), OracleLassoA[[2]][,8]-rep(c_val, N_Sim), 
                         LassoA[[2]][,1]-rep(c_val, N_Sim), LassoA[[2]][,2]-rep(c_val, N_Sim), LassoA[[2]][,3]-rep(c_val, N_Sim), LassoA[[2]][,5]-rep(c_val, N_Sim), LassoA[[2]][,8]-rep(c_val, N_Sim)),
              
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,10), seq(12,16), seq(18,22), seq(24,28), seq(30, 34)),
              col = col_vec_ASMALLT, cex = 0.2, boxwex = 0.5, lwd = 0.8, whisklty = 1)
      title("p =  10",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 8, 15, 21, 26, 32), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 11, lwd = 1)
      abline(h = 17, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 29, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      
      ## Scenario B
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS2/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pB_vec[2],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame(OracleOracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,2]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,5]-rep(c_val, N_Sim), 
                         OracleOracleLSB[[2]][,8]-rep(c_val, N_Sim), OracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleLSB[[2]][,2]-rep(c_val, N_Sim),  OracleLSB[[2]][,3]-rep(c_val, N_Sim),  OracleLSB[[2]][,5], OracleLSB[[2]][,(2+50)]-rep(c_val, N_Sim),  
                         OracleLSB[[2]][,(5+50)],  OracleLSB[[2]][,(2+2*50)]-rep(c_val, N_Sim), OracleLSB[[2]][,(5+2*50)], LSB[[2]][,1]-rep(c_val, N_Sim), LSB[[2]][,2]-rep(c_val, N_Sim), LSB[[2]][,3]-rep(c_val, N_Sim), LSB[[2]][,5], LSB[[2]][,(2+50)]-rep(c_val, N_Sim), LSB[[2]][,(5+50)], LSB[[2]][,(2+2*50)]-rep(c_val, N_Sim), 
                         LSB[[2]][,(5+2*50)], OracleLassoB[[2]][,1]-rep(c_val, N_Sim), OracleLassoB[[2]][,2]-rep(c_val, N_Sim), OracleLassoB[[2]][,3]-rep(c_val, N_Sim), OracleLassoB[[2]][,5], OracleLassoB[[2]][,(2+50)]-rep(c_val, N_Sim), 
                         OracleLassoB[[2]][,(5+50)], OracleLassoB[[2]][,(2+2*50)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+2*50)],  LassoB[[2]][,1]-rep(c_val, N_Sim), LassoB[[2]][,2]-rep(c_val, N_Sim), LassoB[[2]][,3]-rep(c_val, N_Sim), LassoB[[2]][,5], LassoB[[2]][,(2+50)]-rep(c_val, N_Sim), LassoB[[2]][,(5+50)], 
                         LassoB[[2]][,(2+2*50)]-rep(c_val, N_Sim), LassoB[[2]][,(5+2*50)]), 
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title("p =  151",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      
      
      # p = 3
      
      ## Scenario A
      # Boxplot for group representatives for each estimator
      
      pdf(file = paste0("./plots/FigS1/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pA_vec[3],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame( 
        OracleOracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,4]-rep(c_val, N_Sim), 
        OracleOracleLSA[[3]][,6]-rep(c_val, N_Sim), PesaranCCE[[3]][,1]-rep(c_val, N_Sim),  PesaranCCE[[3]][,2]-rep(c_val, N_Sim),  PesaranCCE[[3]][,3]-rep(c_val, N_Sim), PesaranCCE[[3]][,5], PesaranCCE[[3]][,6]-rep(c_val, N_Sim), 
        PesaranCCE[[3]][,9], PesaranCCE[[3]][,10]-rep(c_val, N_Sim), PesaranCCE[[3]][,13], OracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleLSA[[3]][,2]-rep(c_val, N_Sim),  OracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleLSA[[3]][,5], OracleLSA[[3]][,6]-rep(c_val, N_Sim), 
        OracleLSA[[3]][,9], OracleLSA[[3]][,10]-rep(c_val, N_Sim),OracleLSA[[3]][,13], LSA[[3]][,1]-rep(c_val, N_Sim), LSA[[3]][,2]-rep(c_val, N_Sim), LSA[[3]][,3]-rep(c_val, N_Sim), LSA[[3]][,5], LSA[[3]][,6]-rep(c_val, N_Sim), LSA[[3]][,9], LSA[[3]][,10]-rep(c_val, N_Sim), 
        LSA[[3]][,13], OracleLassoA[[3]][,1]-rep(c_val, N_Sim), OracleLassoA[[3]][,2]-rep(c_val, N_Sim), OracleLassoA[[3]][,3]-rep(c_val, N_Sim), OracleLassoA[[3]][,5], OracleLassoA[[3]][,6]-rep(c_val, N_Sim), OracleLassoA[[3]][,9], OracleLassoA[[3]][,10]-rep(c_val, N_Sim), OracleLassoA[[3]][,13], LassoA[[3]][,1]-rep(c_val, N_Sim), LassoA[[3]][,2]-rep(c_val, N_Sim), LassoA[[3]][,3]-rep(c_val, N_Sim), LassoA[[3]][,5], LassoA[[3]][,6]-rep(c_val, N_Sim), LassoA[[3]][,9], 
        LassoA[[3]][,10]-rep(c_val, N_Sim), LassoA[[3]][,13]),
        
        yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
        col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
      title("p =  13",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(h = 41, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      ## Scenario B
      pdf(file = paste0("./plots/FigS2/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pB_vec[3],".pdf", sep = ""), width = 6, height = 6)
      # Boxplot for group representatives for each estimator
      boxplot(data.frame( OracleOracleLSB[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[3]][,2]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[3]][,5]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[3]][,10]-rep(c_val, N_Sim), OracleLSB[[3]][,1]-rep(c_val, N_Sim), OracleLSB[[3]][,2]-rep(c_val, N_Sim),  
                          OracleLSB[[3]][,3]-rep(c_val, N_Sim), OracleLSB[[3]][,5], 
                          OracleLSB[[3]][,(2+100)]-rep(c_val,N_Sim),  
                          OracleLSB[[3]][,(5+100)], OracleLSB[[3]][,(2+2*100)]-rep(c_val, N_Sim), OracleLSB[[3]][,(5+2*100)], LSB[[3]][,1]-rep(c_val, N_Sim), LSB[[3]][,2]-rep(c_val, N_Sim), LSB[[3]][,3]-rep(c_val, N_Sim), 
                          LSB[[3]][,5], LSB[[3]][,(2+100)]-rep(c_val, N_Sim), LSB[[3]][,(5+100)], 
                          LSB[[3]][,(2+2*100)]-rep(c_val, N_Sim), 
                          LSB[[3]][,(5+2*100)], 
                          OracleLassoB[[3]][,1]-rep(c_val, N_Sim), OracleLassoB[[3]][,2]-rep(c_val, N_Sim), 
                          OracleLassoB[[3]][,3]-rep(c_val, N_Sim), OracleLassoB[[3]][,5], 
                          OracleLassoB[[3]][,(2+100)]-rep(c_val, N_Sim), OracleLassoB[[3]][,(5+100)], 
                          OracleLassoB[[3]][,(2+2*100)]-rep(c_val, N_Sim), OracleLassoB[[3]][,(5+2*100)],  LassoB[[3]][,1]-rep(c_val, N_Sim), LassoB[[3]][,2]-rep(c_val, N_Sim), LassoB[[3]][,3]-rep(c_val, N_Sim), 
                          LassoB[[3]][,5], LassoB[[3]][,(2+100)]-rep(c_val, N_Sim), LassoB[[3]][,(5+100)], 
                          LassoB[[3]][,(2+2*100)]-rep(c_val, N_Sim), LassoB[[3]][,(5+2*100)]), 
              
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title("p =  301",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
    }
    
    
    ########################### 1.2 Boxplots for T_dim = 50 ######################  
    if(T_dim == 50){
      # p = 1
      
      ## Scenario A
      # Boxplot for group representatives for each estimator
      col_vec2 <- c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 4))
      col_vec_A <-  c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 5))
      
      pdf(file = paste0("./plots/FigS1/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pA_vec[1],".pdf", sep = ""), width = 6, height = 6)
      
      boxplot(data.frame(OracleOracleLSA[[1]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[1]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[1]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSA[[1]][,6]-rep(c_val, N_Sim), PesaranCCE[[1]][,1]-rep(c_val, N_Sim), 
                         PesaranCCE[[1]][,2]-rep(c_val, N_Sim),  PesaranCCE[[1]][,3]-rep(c_val, N_Sim), PesaranCCE[[1]][,5], PesaranCCE[[1]][,(2+5)]-rep(c_val, N_Sim), 
                         PesaranCCE[[1]][,(5+5)], PesaranCCE[[1]][,(2+2*5)]-rep(c_val, N_Sim),PesaranCCE[[1]][,(5+2*5)], 
                         OracleLSA[[1]][,1]-rep(c_val, N_Sim), OracleLSA[[1]][,2]-rep(c_val, N_Sim),  OracleLSA[[1]][,3]-rep(c_val, N_Sim), OracleLSA[[1]][,5], OracleLSA[[1]][,(2+5)]-rep(c_val, N_Sim), 
                         OracleLSA[[1]][,(5+5)], OracleLSA[[1]][,(2+2*5)]-rep(c_val, N_Sim),OracleLSA[[1]][,(5+2*5)],
                         OracleOracleLSA[[1]][,1]-rep(c_val, N_Sim), 
                         LSA[[1]][,1]-rep(c_val, N_Sim), LSA[[1]][,2]-rep(c_val, N_Sim), LSA[[1]][,3]-rep(c_val, N_Sim), LSA[[1]][,5], LSA[[1]][,(2+5)]-rep(c_val, N_Sim), LSA[[1]][,(5+5)], LSA[[1]][,(2+2*5)]-rep(c_val, N_Sim), 
                         LSA[[1]][,(5+2*5)],
                         OracleLassoA[[1]][,1]-rep(c_val, N_Sim), OracleLassoA[[1]][,2]-rep(c_val, N_Sim), OracleLassoA[[1]][,3]-rep(c_val, N_Sim), OracleLassoA[[1]][,5], OracleLassoA[[1]][,(2+5)]-rep(c_val, N_Sim), 
                         OracleLassoA[[1]][,(5+5)], OracleLassoA[[1]][,(2+2*5)]-rep(c_val, N_Sim), OracleLassoA[[1]][,(5+2*5)], 
                         LassoA[[1]][,1]-rep(c_val, N_Sim), LassoA[[1]][,2]-rep(c_val, N_Sim), LassoA[[1]][,3]-rep(c_val, N_Sim), LassoA[[1]][,5], LassoA[[1]][,(2+5)]-rep(c_val, N_Sim), LassoA[[1]][,(5+5)], 
                         LassoA[[1]][,(2+2*5)]-rep(c_val, N_Sim), LassoA[[1]][,(5+2*5)]),
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
              col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
      title("p =  16",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(h = 41, lwd = 1)
      abline(v = 0)
      
      dev.off()
      
      
      
      
      ## Scenario B
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS2/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pB_vec[1],".pdf", sep = ""), width = 6, height = 6)
      
      boxplot(data.frame(OracleOracleLSB[[1]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[1]][,2]-rep(c_val, N_Sim), OracleOracleLSB[[1]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[1]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSB[[1]][,6]-rep(c_val, N_Sim), OracleLSB[[1]][,1]-rep(c_val, N_Sim), OracleLSB[[1]][,2]-rep(c_val, N_Sim),  OracleLSB[[1]][,3]-rep(c_val, N_Sim), OracleLSB[[1]][,5], OracleLSB[[1]][,(2+30)]-rep(c_val, N_Sim), 
                         OracleLSB[[1]][,(5+30)], OracleLSB[[1]][,(2+2*30)]-rep(c_val, N_Sim),OracleLSB[[1]][,(5+2*30)], LSB[[1]][,1]-rep(c_val, N_Sim), LSB[[1]][,2]-rep(c_val, N_Sim), LSB[[1]][,3]-rep(c_val, N_Sim), LSB[[1]][,5], LSB[[1]][,(2+30)]-rep(c_val, N_Sim), LSB[[1]][,(5+30)], LSB[[1]][,(2+2*30)]-rep(c_val, N_Sim), 
                         LSB[[1]][,(5+2*30)],
                         OracleLassoB[[1]][,1]-rep(c_val, N_Sim), OracleLassoB[[1]][,2]-rep(c_val, N_Sim), OracleLassoB[[1]][,3]-rep(c_val, N_Sim), OracleLassoB[[1]][,5], OracleLassoB[[1]][,(2+30)]-rep(c_val, N_Sim), 
                         OracleLassoB[[1]][,(5+30)], OracleLassoB[[1]][,(2+2*30)]-rep(c_val, N_Sim), OracleLassoB[[1]][,(5+2*30)], 
                         LassoB[[1]][,1]-rep(c_val, N_Sim), LassoB[[1]][,2]-rep(c_val, N_Sim), LassoB[[1]][,3]-rep(c_val, N_Sim), LassoB[[1]][,5], LassoB[[1]][,(2+30)]-rep(c_val, N_Sim), LassoB[[1]][,(5+30)], 
                         LassoB[[1]][,(2+2*30)]-rep(c_val, N_Sim), LassoB[[1]][,(5+2*30)]),
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title("p =  91",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      ## Scenario C
      pdf(file = paste0("./plots/FigS3/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pC,".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame(OracleOracleLSC[[1]][,1]-rep(c_val, N_Sim), OracleOracleLSC[[1]][,2]-rep(c_val, N_Sim), OracleOracleLSC[[1]][,3]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[1]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[1]][,6]-rep(c_val, N_Sim), OracleLassoC[[1]][,1]-rep(c_val, N_Sim), OracleLassoC[[1]][,2]-rep(c_val, N_Sim), OracleLassoC[[1]][,3]-rep(c_val, N_Sim), OracleLassoC[[1]][,5], 
                         OracleLassoC[[1]][,(2+1000)]-rep(c_val, N_Sim), 
                         OracleLassoC[[1]][,(5+1000)], OracleLassoC[[1]][,(2+2*1000)]-rep(c_val, N_Sim), OracleLassoC[[1]][,(5+2*1000)], 
                         LassoC[[1]][,1]-rep(c_val, N_Sim), LassoC[[1]][,2]-rep(c_val, N_Sim), LassoC[[1]][,3]-rep(c_val, N_Sim), LassoC[[1]][,5], LassoC[[1]][,(2+1000)]-rep(c_val, N_Sim), 
                         LassoC[[1]][,(5+1000)], LassoC[[1]][,(2+2*1000)]-rep(c_val, N_Sim), LassoC[[1]][,(5+2*1000)]),
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22)), col = col_vec2, cex = 0.2, boxwex = 0.4, lwd = 0.8, whisklty = 1)
      title("p =  3001",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19), labels = c(TeX(r"(LS-O$^2$)"), "L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      # p = 2
      
      ## Scenario A
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS1/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pA_vec[2],".pdf", sep = ""), width = 6, height = 6)
      
      boxplot(data.frame(OracleOracleLSA[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[2]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[2]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSA[[2]][,6]-rep(c_val, N_Sim),  
                         PesaranCCE[[2]][,1]-rep(c_val, N_Sim), PesaranCCE[[2]][,2]-rep(c_val, N_Sim),  PesaranCCE[[2]][,3]-rep(c_val, N_Sim), PesaranCCE[[2]][,5], PesaranCCE[[2]][,(2+10)]-rep(c_val, N_Sim), 
                         PesaranCCE[[2]][,(5+10)], PesaranCCE[[2]][,(2+2*10)]-rep(c_val, N_Sim),PesaranCCE[[2]][,(5+2*10)],  OracleLSA[[2]][,1]-rep(c_val, N_Sim), OracleLSA[[2]][,2]-rep(c_val, N_Sim),  OracleLSA[[2]][,3]-rep(c_val, N_Sim), OracleLSA[[2]][,5], OracleLSA[[2]][,(2+10)]-rep(c_val, N_Sim), 
                         OracleLSA[[2]][,(5+10)], OracleLSA[[2]][,(2+2*10)]-rep(c_val, N_Sim),OracleLSA[[2]][,(5+2*10)], LSA[[2]][,1]-rep(c_val, N_Sim), LSA[[2]][,2]-rep(c_val, N_Sim), LSA[[2]][,3]-rep(c_val, N_Sim), LSA[[2]][,5], LSA[[2]][,(2+10)]-rep(c_val, N_Sim), LSA[[2]][,(5+10)], LSA[[2]][,(2+2*10)]-rep(c_val, N_Sim), 
                         LSA[[2]][,(5+2*10)], OracleLassoA[[2]][,1]-rep(c_val, N_Sim), OracleLassoA[[2]][,2]-rep(c_val, N_Sim), OracleLassoA[[2]][,3]-rep(c_val, N_Sim), OracleLassoA[[2]][,10], OracleLassoA[[2]][,(2+10)]-rep(c_val, N_Sim), 
                         OracleLassoA[[2]][,(5+10)], OracleLassoA[[2]][,(2+2*10)]-rep(c_val, N_Sim), OracleLassoA[[2]][,(5+2*10)], LassoA[[2]][,1]-rep(c_val, N_Sim), LassoA[[2]][,2]-rep(c_val, N_Sim), LassoA[[2]][,3]-rep(c_val, N_Sim), LassoA[[2]][,5], LassoA[[2]][,(2+10)]-rep(c_val, N_Sim), LassoA[[2]][,(5+10)], 
                         LassoA[[2]][,(2+2*10)]-rep(c_val, N_Sim), LassoA[[2]][,(5+2*10)]),
              
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
              col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
      title("p =  31",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(h = 41, lwd = 1)
      abline(v = 0)
      
      dev.off()
      
      
      
      ## Scenario B
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS2/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pB_vec[2],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame(OracleOracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,2]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,5]-rep(c_val, N_Sim), 
                         OracleOracleLSB[[2]][,8]-rep(c_val, N_Sim), OracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleLSB[[2]][,2]-rep(c_val, N_Sim),  OracleLSB[[2]][,3]-rep(c_val, N_Sim),  OracleLSB[[2]][,5], OracleLSB[[2]][,(2+150)]-rep(c_val, N_Sim),  
                         OracleLSB[[2]][,(5+150)],  OracleLSB[[2]][,(2+2*150)]-rep(c_val, N_Sim), OracleLSB[[2]][,(5+2*150)], LSB[[2]][,1]-rep(c_val, N_Sim), LSB[[2]][,2]-rep(c_val, N_Sim), LSB[[2]][,3]-rep(c_val, N_Sim), LSB[[2]][,5], LSB[[2]][,(2+150)]-rep(c_val, N_Sim), LSB[[2]][,(5+150)], LSB[[2]][,(2+2*150)]-rep(c_val, N_Sim), 
                         LSB[[2]][,(5+2*150)], OracleLassoB[[2]][,1]-rep(c_val, N_Sim), OracleLassoB[[2]][,2]-rep(c_val, N_Sim), OracleLassoB[[2]][,3]-rep(c_val, N_Sim), OracleLassoB[[2]][,5], OracleLassoB[[2]][,(2+150)]-rep(c_val, N_Sim), 
                         OracleLassoB[[2]][,(5+150)], OracleLassoB[[2]][,(2+2*150)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+2*150)],  LassoB[[2]][,1]-rep(c_val, N_Sim), LassoB[[2]][,2]-rep(c_val, N_Sim), LassoB[[2]][,3]-rep(c_val, N_Sim), LassoB[[2]][,5], LassoB[[2]][,(2+150)]-rep(c_val, N_Sim), LassoB[[2]][,(5+150)], 
                         LassoB[[2]][,(2+2*150)]-rep(c_val, N_Sim), LassoB[[2]][,(5+2*150)]), 
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title("p =  451",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      
      # p = 3
      
      ## Scenario A
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS1/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pA_vec[3],".pdf", sep = ""), width = 6, height = 6)
      
      boxplot(data.frame( 
        OracleOracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,4]-rep(c_val, N_Sim), 
        OracleOracleLSA[[3]][,6]-rep(c_val, N_Sim), PesaranCCE[[3]][,1]-rep(c_val, N_Sim),  PesaranCCE[[3]][,2]-rep(c_val, N_Sim),  PesaranCCE[[3]][,3]-rep(c_val, N_Sim), PesaranCCE[[3]][,5], PesaranCCE[[3]][,(2+15)]-rep(c_val, N_Sim), 
        PesaranCCE[[3]][,(5+15)], PesaranCCE[[3]][,(2+2*15)]-rep(c_val, N_Sim),PesaranCCE[[3]][,(5+2*15)], OracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleLSA[[3]][,2]-rep(c_val, N_Sim),  OracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleLSA[[3]][,5], OracleLSA[[3]][,(2+15)]-rep(c_val, N_Sim), 
        OracleLSA[[3]][,(5+15)], OracleLSA[[3]][,(2+2*15)]-rep(c_val, N_Sim),OracleLSA[[3]][,(5+2*15)], LSA[[3]][,1]-rep(c_val, N_Sim), LSA[[3]][,2]-rep(c_val, N_Sim), LSA[[3]][,3]-rep(c_val, N_Sim), LSA[[3]][,5], LSA[[3]][,(2+15)]-rep(c_val, N_Sim), LSA[[3]][,(5+15)], LSA[[3]][,(2+2*15)]-rep(c_val, N_Sim), 
        LSA[[3]][,(5+2*15)], OracleLassoA[[3]][,1]-rep(c_val, N_Sim), OracleLassoA[[3]][,2]-rep(c_val, N_Sim), OracleLassoA[[3]][,3]-rep(c_val, N_Sim), OracleLassoA[[3]][,15], OracleLassoA[[3]][,(2+15)]-rep(c_val, N_Sim), OracleLassoA[[3]][,(5+15)], OracleLassoA[[3]][,(2+2*15)]-rep(c_val, N_Sim), OracleLassoA[[3]][,(5+2*15)], LassoA[[3]][,1]-rep(c_val, N_Sim), LassoA[[3]][,2]-rep(c_val, N_Sim), LassoA[[3]][,3]-rep(c_val, N_Sim), LassoA[[3]][,5], LassoA[[3]][,(2+15)]-rep(c_val, N_Sim), LassoA[[3]][,(5+15)], 
        LassoA[[3]][,(2+2*15)]-rep(c_val, N_Sim), LassoA[[3]][,(5+2*15)]),
        
        yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
        col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
      title("p =  46",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(h = 41, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      ## Scenario B
      pdf(file = paste0("./plots/FigS2/", "T_dim", T_dim, "gamma", gamma,
                        "HAC", HAC, "p", pB_vec[3],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame( OracleOracleLSB[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[3]][,2]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[3]][,5]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[3]][,10]-rep(c_val, N_Sim), OracleLSB[[3]][,1]-rep(c_val, N_Sim), OracleLSB[[3]][,2]-rep(c_val, N_Sim),  
                          OracleLSB[[3]][,3]-rep(c_val, N_Sim), OracleLSB[[3]][,5], 
                          OracleLSB[[3]][,(2+300)]-rep(c_val,N_Sim),  
                          OracleLSB[[3]][,(5+300)], OracleLSB[[3]][,(2+2*300)]-rep(c_val, N_Sim), OracleLSB[[3]][,(5+2*300)], LSB[[3]][,1]-rep(c_val, N_Sim), LSB[[3]][,2]-rep(c_val, N_Sim), LSB[[3]][,3]-rep(c_val, N_Sim), 
                          LSB[[3]][,5], LSB[[3]][,(2+300)]-rep(c_val, N_Sim), LSB[[3]][,(5+300)], 
                          LSB[[3]][,(2+2*300)]-rep(c_val, N_Sim), 
                          LSB[[3]][,(5+2*300)], 
                          OracleLassoB[[3]][,1]-rep(c_val, N_Sim), OracleLassoB[[3]][,2]-rep(c_val, N_Sim), 
                          OracleLassoB[[3]][,3]-rep(c_val, N_Sim), OracleLassoB[[3]][,5], 
                          OracleLassoB[[3]][,(2+300)]-rep(c_val, N_Sim), OracleLassoB[[3]][,(5+300)], 
                          OracleLassoB[[3]][,(2+2*300)]-rep(c_val, N_Sim), OracleLassoB[[3]][,(5+2*300)],  LassoB[[3]][,1]-rep(c_val, N_Sim), LassoB[[3]][,2]-rep(c_val, N_Sim), LassoB[[3]][,3]-rep(c_val, N_Sim), 
                          LassoB[[3]][,5], LassoB[[3]][,(2+300)]-rep(c_val, N_Sim), LassoB[[3]][,(5+300)], 
                          LassoB[[3]][,(2+2*300)]-rep(c_val, N_Sim), LassoB[[3]][,(5+2*300)]), 
              
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title("p =  901",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
    }
  }
    
    
    
  }
  
  
  
  ################# 2. Inference: Tables S.3,S.4,S.5 and S.8-S.13 #################
  { if(gamma == 1){
    rejA <- numeric(3)
    rejB <- numeric(3)
    rejC <- 0
    
    rej <- matrix(NA, nrow = 7, ncol = 9)
    
    
    count <- 1
    for (signal in 1:3) {
      
      for (alpha in c(0.01, 0.05, 0.1)) {
        
        z_alpha <- qnorm(1-alpha/2)
        
        for (p in 1:3) {
          
          tmpA <- 0
          tmpB <- 0
          tmpC <- 0
          
          for (sim in 1:N_Sim) {
            
            tmpA <- tmpA + as.numeric(abs(results_Factors[[sim]][[p]][[1]][["TestStatistic"]][signal]) > 
                                        z_alpha)
            tmpB <- tmpB + as.numeric(abs(results_Factors[[sim]][[p]][[2]][["TestStatistic"]][signal]) > 
                                        z_alpha)
            if(p == 3){
              tmpC <- tmpC + as.numeric(abs(results_Factors[[sim]][[p]][[3]][["TestStatistic"]][signal]) > 
                                          z_alpha)
            }
            
          }
          
          rejA[p] <- tmpA
          rejB[p] <- tmpB
          if(p == 3){
            rejC <- tmpC
          }
        }
        tmp <- c(rejA, rejB, rejC)/N_Sim
        
        rej[, count] <- tmp
        
        count <- count + 1
      }
    }
    write.csv(rej, paste0("./tables/Inference",
                          "T_dim", T_dim, "gamma", gamma, "HAC", HAC,".csv", sep = ""))
  }
  }
  
  
  
  
  ################# 3. Estimated Factors Khat: Tables  S.6 and S.7 ###############
  {if(HAC == 1){
    
    estimated_KhatA <- matrix(NA, nrow = 3, ncol = 9)
    estimated_KhatB <- matrix(NA, nrow = 3, ncol = 9)
    estimated_KhatC <- numeric(9)
    
    
    for (p in 1:3) {
      tmpestimated_KhatA <- numeric(9)
      tmpestimated_KhatB <- numeric(9)
      tmpestimated_KhatC <- numeric(9)
      for (sim in 1:N_Sim){
        for (j in 1:9) {
          if(results_Factors[[sim]][[p]][[1]][["Lasso"]][["K_hat"]] == j){
            tmpestimated_KhatA[j] <- tmpestimated_KhatA[j] + results_Factors[[sim]][[p]][[1]][["Lasso"]][["K_hat"]]
          }
          
          if(results_Factors[[sim]][[p]][[2]][["Lasso"]][["K_hat"]] == j){
            tmpestimated_KhatB[j] <- tmpestimated_KhatB[j] + results_Factors[[sim]][[p]][[2]][["Lasso"]][["K_hat"]]
          }
          if(p == 3){
            if(results_Factors[[sim]][[p]][[3]][["Lasso"]][["K_hat"]]  == j){
              tmpestimated_KhatC[j] <- tmpestimated_KhatC[j] + results_Factors[[sim]][[p]][[3]][["Lasso"]][["K_hat"]]
            }
          }
        }
      }
      estimated_KhatA[p, ] <- tmpestimated_KhatA / c(1,2,3,4,5,6,7,8,9)
      estimated_KhatB[p, ]<- tmpestimated_KhatB / c(1,2,3,4,5,6,7,8,9)
      
      if(p == 3){
        estimated_KhatC <- tmpestimated_KhatC / c(1,2,3,4,5,6,7,8,9)
      }
    }
    df_estimated_factors <- rbind(estimated_KhatA, estimated_KhatB, estimated_KhatC)
    write.csv(df_estimated_factors, paste0("./tables/Factors",
                                           "T_dim", T_dim, "gamma", gamma, "HAC", HAC,".csv", sep = ""))
  }
  }
  
  
  
  ################### 4. Accuracy Plot for Khat = 6: Figure S.4 ###################
  {if(HAC == 1 & gamma == 1){
    
    LassoA <-  vector(mode = 'list', length = 3)
    OracleLassoA <-  vector(mode = 'list', length = 3)
    LSA <- vector(mode = 'list', length = 3)
    OracleLSA <- vector(mode = 'list', length = 3)
    OracleOracleLSA <-  vector(mode = 'list', length = 3)
    PesaranCCE <-  vector(mode = 'list', length = 3)
    
    LassoB <-  vector(mode = 'list', length = 3)
    OracleLassoB <- vector(mode = 'list', length = 3)
    LSB <-  vector(mode = 'list', length = 3)
    OracleLSB <-  vector(mode = 'list', length = 3)
    OracleOracleLSB <-  vector(mode = 'list', length = 3)
    
    
    LassoC <-  vector(mode = 'list', length = 3)
    OracleLassoC <- vector(mode = 'list', length = 3)
    OracleOracleLSC <-  vector(mode = 'list', length = 3)
    
    
    tmpLassoA <- numeric(pA_vec[3])
    tmpOracleLassoA <- numeric(pA_vec[3])
    tmpLSA <- numeric(pA_vec[3])
    tmpOracleLSA <- numeric(pA_vec[3])
    tmpOracleOracleLSA <- numeric(10)
    tmpPesaranCCE <- numeric(pA_vec[3])
    
    tmpLassoB <- numeric(pB_vec[2])
    tmpOracleLassoB <- numeric(pB_vec[2])
    tmpLSB <- numeric(pB_vec[2])
    tmpOracleLSB<- numeric(pB_vec[2])
    tmpOracleOracleLSB <- numeric(10)
    
    
    tmpLassoC <- numeric(pC)
    tmpOracleLassoC <- numeric(pC)
    tmpOracleOracleLSC <- numeric(10)
    
    
    for (sim in 1:N_Sim) {
      tmpLassoA <- rbind(tmpLassoA, results_Factors[[sim]][[3]][[1]][["HighKLasso"]][["coefs"]])
      tmpOracleLassoA <- rbind(tmpOracleLassoA, results_Factors[[sim]][[3]][[1]][["OracleLasso"]])
      tmpLSA<-rbind(tmpLSA, results_Factors[[sim]][[3]][[1]][["HighKLS"]][["coefs"]])
      tmpOracleLSA <- rbind(tmpOracleLSA, results_Factors[[sim]][[3]][[1]][["OracleLS"]])
      tmpOracleOracleLSA<- rbind(tmpOracleOracleLSA, results_Factors[[sim]][[3]][[1]][["OracleOracleLS"]])
      tmpPesaranCCE <- rbind(tmpPesaranCCE, results_Factors[[sim]][[3]][[1]][["PesaranCCE"]])
      
      tmpLassoB <- rbind(tmpLassoB, results_Factors[[sim]][[2]][[2]][["HighKLasso"]][["coefs"]])
      tmpOracleLassoB <- rbind(tmpOracleLassoB, results_Factors[[sim]][[2]][[2]][["OracleLasso"]])
      tmpLSB <- rbind(tmpLSB, results_Factors[[sim]][[2]][[2]][["HighKLS"]][["coefs"]])
      tmpOracleLSB<- rbind(tmpOracleLSB,  results_Factors[[sim]][[2]][[2]][["OracleLS"]])
      tmpOracleOracleLSB <- rbind(tmpOracleOracleLSB, results_Factors[[sim]][[2]][[2]][["OracleOracleLS"]])
      
      tmpLassoC <- rbind(tmpLassoC, results_Factors[[sim]][[3]][[3]][["HighKLasso"]][["coefs"]])
      tmpOracleLassoC <- rbind(tmpOracleLassoC,results_Factors[[sim]][[3]][[3]][["OracleLasso"]])
      tmpOracleOracleLSC <- rbind(tmpOracleOracleLSC, results_Factors[[sim]][[3]][[3]][["OracleOracleLS"]])
      
    }
    LassoA[[3]] <-  tmpLassoA[-1,]
    OracleLassoA[[3]] <- tmpOracleLassoA[-1,]
    LSA[[3]] <- tmpLSA[-1,]
    OracleLSA[[3]] <- tmpOracleLSA[-1,]
    OracleOracleLSA[[3]] <- tmpOracleOracleLSA[-1,]
    PesaranCCE[[3]] <- tmpPesaranCCE[-1,]
    
    LassoB[[2]] <- tmpLassoB[-1,]
    OracleLassoB[[2]] <- tmpOracleLassoB[-1,]
    LSB[[2]] <- tmpLSB[-1,]
    OracleLSB[[2]] <- tmpOracleLSB[-1,]
    OracleOracleLSB[[2]] <- tmpOracleOracleLSB[-1,]
    
    LassoC[[3]] <- tmpLassoC[-1,]
    OracleLassoC[[3]] <- tmpOracleLassoC[-1,]
    OracleOracleLSC[[3]] <- tmpOracleOracleLSC[-1,]
    
    
    if(T_dim == 15){
      
      col_vec2 <- c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 4))
      col_vec_A <-  c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 5))
      
      
      pdf(file = paste0("./plots/FigS4/","Robustness", "Khat", 6, "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pC,".pdf", sep = ""), width = 6, height = 6)
      
      ## Scenario C
      boxplot(data.frame(OracleOracleLSC[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,3]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[3]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[3]][,6]-rep(c_val, N_Sim), OracleLassoC[[3]][,1]-rep(c_val, N_Sim), OracleLassoC[[3]][,2]-rep(c_val, N_Sim), OracleLassoC[[3]][,3]-rep(c_val, N_Sim), OracleLassoC[[3]][,5], 
                         OracleLassoC[[3]][,(2+300)]-rep(c_val, N_Sim), 
                         OracleLassoC[[3]][,(5+300)], OracleLassoC[[3]][,(2+2*300)]-rep(c_val, N_Sim), OracleLassoC[[3]][,(5+2*300)], 
                         LassoC[[3]][,1]-rep(c_val, N_Sim), LassoC[[3]][,2]-rep(c_val, N_Sim), LassoC[[3]][,3]-rep(c_val, N_Sim), LassoC[[3]][,5], LassoC[[3]][,(2+300)]-rep(c_val, N_Sim), 
                         LassoC[[3]][,(5+300)], LassoC[[3]][,(2+2*300)]-rep(c_val, N_Sim), LassoC[[3]][,(5+2*300)]),
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22)), col = col_vec2, cex = 0.2, boxwex = 0.4, lwd = 0.8, whisklty = 1)
      title("p =  901",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19), labels = c(TeX(r"(LS-O$^2$)"), "L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      ## Scenario A
      pdf(file = paste0("./plots/FigS4/","Robustness", "Khat", 6, "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pA_vec[3],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame( 
        OracleOracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,4]-rep(c_val, N_Sim), 
        OracleOracleLSA[[3]][,6]-rep(c_val, N_Sim), PesaranCCE[[3]][,1]-rep(c_val, N_Sim),  PesaranCCE[[3]][,2]-rep(c_val, N_Sim),  PesaranCCE[[3]][,3]-rep(c_val, N_Sim), PesaranCCE[[3]][,5], PesaranCCE[[3]][,6]-rep(c_val, N_Sim), 
        PesaranCCE[[3]][,9], PesaranCCE[[3]][,10]-rep(c_val, N_Sim), PesaranCCE[[3]][,13], OracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleLSA[[3]][,2]-rep(c_val, N_Sim),  OracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleLSA[[3]][,5], OracleLSA[[3]][,6]-rep(c_val, N_Sim), 
        OracleLSA[[3]][,9], OracleLSA[[3]][,10]-rep(c_val, N_Sim),OracleLSA[[3]][,13], LSA[[3]][,1]-rep(c_val, N_Sim), LSA[[3]][,2]-rep(c_val, N_Sim), LSA[[3]][,3]-rep(c_val, N_Sim), LSA[[3]][,5], LSA[[3]][,6]-rep(c_val, N_Sim), LSA[[3]][,9], LSA[[3]][,10]-rep(c_val, N_Sim), 
        LSA[[3]][,13], OracleLassoA[[3]][,1]-rep(c_val, N_Sim), OracleLassoA[[3]][,2]-rep(c_val, N_Sim), OracleLassoA[[3]][,3]-rep(c_val, N_Sim), OracleLassoA[[3]][,5], OracleLassoA[[3]][,6]-rep(c_val, N_Sim), OracleLassoA[[3]][,9], OracleLassoA[[3]][,10]-rep(c_val, N_Sim), OracleLassoA[[3]][,13], LassoA[[3]][,1]-rep(c_val, N_Sim), LassoA[[3]][,2]-rep(c_val, N_Sim), LassoA[[3]][,3]-rep(c_val, N_Sim), LassoA[[3]][,5], LassoA[[3]][,6]-rep(c_val, N_Sim), LassoA[[3]][,9], 
        LassoA[[3]][,10]-rep(c_val, N_Sim), LassoA[[3]][,13]),
        
        yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
        col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
      title("p =  13",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(h = 41, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      ## Scenario B
      pdf(file = paste0("./plots/FigS4/","Robustness", "Khat", 6, "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pB_vec[2],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame( OracleOracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,2]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,5]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[2]][,10]-rep(c_val, N_Sim), OracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleLSB[[2]][,2]-rep(c_val, N_Sim),  
                          OracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleLSB[[2]][,5], 
                          OracleLSB[[2]][,(2+50)]-rep(c_val,N_Sim),  
                          OracleLSB[[2]][,(5+50)], OracleLSB[[2]][,(2+2*50)]-rep(c_val, N_Sim), OracleLSB[[2]][,(5+2*50)], LSB[[2]][,1]-rep(c_val, N_Sim), LSB[[2]][,2]-rep(c_val, N_Sim), LSB[[2]][,3]-rep(c_val, N_Sim), 
                          LSB[[2]][,5], LSB[[2]][,(2+50)]-rep(c_val, N_Sim), LSB[[2]][,(5+50)], 
                          LSB[[2]][,(2+2*50)]-rep(c_val, N_Sim), 
                          LSB[[2]][,(5+2*50)], 
                          OracleLassoB[[2]][,1]-rep(c_val, N_Sim), OracleLassoB[[2]][,2]-rep(c_val, N_Sim), 
                          OracleLassoB[[2]][,3]-rep(c_val, N_Sim), OracleLassoB[[2]][,5], 
                          OracleLassoB[[2]][,(2+50)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+50)], 
                          OracleLassoB[[2]][,(2+2*50)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+2*50)],  LassoB[[2]][,1]-rep(c_val, N_Sim), LassoB[[2]][,2]-rep(c_val, N_Sim), LassoB[[2]][,3]-rep(c_val, N_Sim), 
                          LassoB[[2]][,5], LassoB[[2]][,(2+50)]-rep(c_val, N_Sim), LassoB[[2]][,(5+50)], 
                          LassoB[[2]][,(2+2*50)]-rep(c_val, N_Sim), LassoB[[2]][,(5+2*50)]), 
              
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title( "p = 151" , line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
    }
    
    if(T_dim == 50){
      col_vec2 <- c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 4))
      col_vec_A <-  c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 5))
      
      
      ## Scenario A
      # Boxplot for group representatives for each estimator
      pdf(file = paste0("./plots/FigS4/","Robustness", "Khat", 6, "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pA_vec[3],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame( 
        OracleOracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,4]-rep(c_val, N_Sim), 
        OracleOracleLSA[[3]][,6]-rep(c_val, N_Sim), PesaranCCE[[3]][,1]-rep(c_val, N_Sim),  PesaranCCE[[3]][,2]-rep(c_val, N_Sim),  PesaranCCE[[3]][,3]-rep(c_val, N_Sim), PesaranCCE[[3]][,5], PesaranCCE[[3]][,(2+15)]-rep(c_val, N_Sim), 
        PesaranCCE[[3]][,(5+15)], PesaranCCE[[3]][,(2+2*15)]-rep(c_val, N_Sim),PesaranCCE[[3]][,(5+2*15)], OracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleLSA[[3]][,2]-rep(c_val, N_Sim),  OracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleLSA[[3]][,5], OracleLSA[[3]][,(2+15)]-rep(c_val, N_Sim), 
        OracleLSA[[3]][,(5+15)], OracleLSA[[3]][,(2+2*15)]-rep(c_val, N_Sim),OracleLSA[[3]][,(5+2*15)], LSA[[3]][,1]-rep(c_val, N_Sim), LSA[[3]][,2]-rep(c_val, N_Sim), LSA[[3]][,3]-rep(c_val, N_Sim), LSA[[3]][,5], LSA[[3]][,(2+15)]-rep(c_val, N_Sim), LSA[[3]][,(5+15)], LSA[[3]][,(2+2*15)]-rep(c_val, N_Sim), 
        LSA[[3]][,(5+2*15)], OracleLassoA[[3]][,1]-rep(c_val, N_Sim), OracleLassoA[[3]][,2]-rep(c_val, N_Sim), OracleLassoA[[3]][,3]-rep(c_val, N_Sim), OracleLassoA[[3]][,15], OracleLassoA[[3]][,(2+15)]-rep(c_val, N_Sim), OracleLassoA[[3]][,(5+15)], OracleLassoA[[3]][,(2+2*15)]-rep(c_val, N_Sim), OracleLassoA[[3]][,(5+2*15)], LassoA[[3]][,1]-rep(c_val, N_Sim), LassoA[[3]][,2]-rep(c_val, N_Sim), LassoA[[3]][,3]-rep(c_val, N_Sim), LassoA[[3]][,5], LassoA[[3]][,(2+15)]-rep(c_val, N_Sim), LassoA[[3]][,(5+15)], 
        LassoA[[3]][,(2+2*15)]-rep(c_val, N_Sim), LassoA[[3]][,(5+2*15)]),
        
        yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
        col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
      title("p =  46",line=0.3, cex.main = 1)
      axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(h = 41, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      ## Scenario B
      pdf(file = paste0("./plots/FigS4/","Robustness", "Khat", 6, "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pB_vec[2],".pdf", sep = ""), width = 6, height = 6)
      boxplot(data.frame( OracleOracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,2]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,5]-rep(c_val, N_Sim), 
                          OracleOracleLSB[[2]][,10]-rep(c_val, N_Sim), OracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleLSB[[2]][,2]-rep(c_val, N_Sim),  
                          OracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleLSB[[2]][,5], 
                          OracleLSB[[2]][,(2+150)]-rep(c_val,N_Sim),  
                          OracleLSB[[2]][,(5+150)], OracleLSB[[2]][,(2+2*150)]-rep(c_val, N_Sim), OracleLSB[[2]][,(5+2*150)], LSB[[2]][,1]-rep(c_val, N_Sim), LSB[[2]][,2]-rep(c_val, N_Sim), LSB[[2]][,3]-rep(c_val, N_Sim), 
                          LSB[[2]][,5], LSB[[2]][,(2+150)]-rep(c_val, N_Sim), LSB[[2]][,(5+150)], 
                          LSB[[2]][,(2+2*150)]-rep(c_val, N_Sim), 
                          LSB[[2]][,(5+2*150)], 
                          OracleLassoB[[2]][,1]-rep(c_val, N_Sim), OracleLassoB[[2]][,2]-rep(c_val, N_Sim), 
                          OracleLassoB[[2]][,3]-rep(c_val, N_Sim), OracleLassoB[[2]][,5], 
                          OracleLassoB[[2]][,(2+150)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+150)], 
                          OracleLassoB[[2]][,(2+2*150)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+2*150)],  LassoB[[2]][,1]-rep(c_val, N_Sim), LassoB[[2]][,2]-rep(c_val, N_Sim), LassoB[[2]][,3]-rep(c_val, N_Sim), 
                          LassoB[[2]][,5], LassoB[[2]][,(2+150)]-rep(c_val, N_Sim), LassoB[[2]][,(5+150)], 
                          LassoB[[2]][,(2+2*150)]-rep(c_val, N_Sim), LassoB[[2]][,(5+2*150)]), 
              
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
              col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
      title("p =  451",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(h = 23, lwd = 1)
      abline(h = 32, lwd = 1)
      abline(v = 0)
      dev.off()
      
      
      
      
      ## Scenario C
      pdf(file = paste0("./plots/FigS4/","Robustness", "Khat", 6, "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pC,".pdf", sep = ""), width = 6, height = 6)
      
      boxplot(data.frame(OracleOracleLSC[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,3]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[3]][,4]-rep(c_val, N_Sim), 
                         OracleOracleLSC[[3]][,6]-rep(c_val, N_Sim), OracleLassoC[[3]][,1]-rep(c_val, N_Sim), OracleLassoC[[3]][,2]-rep(c_val, N_Sim), OracleLassoC[[3]][,3]-rep(c_val, N_Sim), OracleLassoC[[3]][,5], 
                         OracleLassoC[[3]][,(2+1000)]-rep(c_val, N_Sim), 
                         OracleLassoC[[3]][,(5+1000)], OracleLassoC[[3]][,(2+2*1000)]-rep(c_val, N_Sim), OracleLassoC[[3]][,(5+2*1000)], 
                         LassoC[[3]][,1]-rep(c_val, N_Sim), LassoC[[3]][,2]-rep(c_val, N_Sim), LassoC[[3]][,3]-rep(c_val, N_Sim), LassoC[[3]][,5], LassoC[[3]][,(2+1000)]-rep(c_val, N_Sim), 
                         LassoC[[3]][,(5+1000)], LassoC[[3]][,(2+2*1000)]-rep(c_val, N_Sim), LassoC[[3]][,(5+2*1000)]),
              yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22)), col = col_vec2, cex = 0.2, boxwex = 0.4, lwd = 0.8, whisklty = 1)
      title("p =  3001",line=0.3, cex.main = 1)
      axis(2, las = 0, at = c(2, 10, 19), labels = c(TeX(r"(LS-O$^2$)"), "L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
      abline(h = 5, lwd = 1)
      abline(h = 14, lwd = 1)
      abline(v = 0)
      dev.off()
    }
  }
    
  }
  
  ############# 5. Accuracy Plot for gamma = sqrt(10): Figure S.5 ###################
  {  
    if(HAC == 1 & gamma == sqrt(10)){
      
      LassoA <-  vector(mode = 'list', length = 3)
      OracleLassoA <-  vector(mode = 'list', length = 3)
      LSA <- vector(mode = 'list', length = 3)
      OracleLSA <- vector(mode = 'list', length = 3)
      OracleOracleLSA <-  vector(mode = 'list', length = 3)
      PesaranCCE <-  vector(mode = 'list', length = 3)
      
      LassoB <-  vector(mode = 'list', length = 3)
      OracleLassoB <- vector(mode = 'list', length = 3)
      LSB <-  vector(mode = 'list', length = 3)
      OracleLSB <-  vector(mode = 'list', length = 3)
      OracleOracleLSB <-  vector(mode = 'list', length = 3)
      
      
      LassoC <-  vector(mode = 'list', length = 3)
      OracleLassoC <- vector(mode = 'list', length = 3)
      OracleOracleLSC <-  vector(mode = 'list', length = 3)
      
      
      tmpLassoA <- numeric(pA_vec[3])
      tmpOracleLassoA <- numeric(pA_vec[3])
      tmpLSA <- numeric(pA_vec[3])
      tmpOracleLSA <- numeric(pA_vec[3])
      tmpOracleOracleLSA <- numeric(10)
      tmpPesaranCCE <- numeric(pA_vec[3])
      
      tmpLassoB <- numeric(pB_vec[2])
      tmpOracleLassoB <- numeric(pB_vec[2])
      tmpLSB <- numeric(pB_vec[2])
      tmpOracleLSB<- numeric(pB_vec[2])
      tmpOracleOracleLSB <- numeric(10)
      
      
      tmpLassoC <- numeric(pC)
      tmpOracleLassoC <- numeric(pC)
      tmpOracleOracleLSC <- numeric(10)
      
      
      for (sim in 1:N_Sim) {
        tmpLassoA <- rbind(tmpLassoA, results_Factors[[sim]][[3]][[1]][["Lasso"]]$coefs)
        tmpOracleLassoA <- rbind(tmpOracleLassoA, results_Factors[[sim]][[3]][[1]][["OracleLasso"]])
        tmpLSA<-rbind(tmpLSA, results_Factors[[sim]][[3]][[1]]$LS)
        tmpOracleLSA <- rbind(tmpOracleLSA, results_Factors[[sim]][[3]][[1]][["OracleLS"]])
        tmpOracleOracleLSA<- rbind(tmpOracleOracleLSA, results_Factors[[sim]][[3]][[1]][["OracleOracleLS"]])
        tmpPesaranCCE <- rbind(tmpPesaranCCE, results_Factors[[sim]][[3]][[1]][["PesaranCCE"]])
        
        tmpLassoB <- rbind(tmpLassoB, results_Factors[[sim]][[2]][[2]][["Lasso"]]$coefs)
        tmpOracleLassoB <- rbind(tmpOracleLassoB, results_Factors[[sim]][[2]][[2]][["OracleLasso"]])
        tmpLSB <- rbind(tmpLSB, results_Factors[[sim]][[2]][[2]]$LS)
        tmpOracleLSB<- rbind(tmpOracleLSB,  results_Factors[[sim]][[2]][[2]][["OracleLS"]])
        tmpOracleOracleLSB <- rbind(tmpOracleOracleLSB, results_Factors[[sim]][[2]][[2]][["OracleOracleLS"]])
        
        tmpLassoC <- rbind(tmpLassoC, results_Factors[[sim]][[3]][[3]][["Lasso"]]$coefs)
        tmpOracleLassoC <- rbind(tmpOracleLassoC,results_Factors[[sim]][[3]][[3]][["OracleLasso"]])
        tmpOracleOracleLSC <- rbind(tmpOracleOracleLSC, results_Factors[[sim]][[3]][[3]][["OracleOracleLS"]])
        
        
      }
      LassoA[[3]] <-  tmpLassoA[-1,]
      OracleLassoA[[3]] <- tmpOracleLassoA[-1,]
      LSA[[3]] <- tmpLSA[-1,]
      OracleLSA[[3]] <- tmpOracleLSA[-1,]
      OracleOracleLSA[[3]] <- tmpOracleOracleLSA[-1,]
      PesaranCCE[[3]] <- tmpPesaranCCE[-1,]
      
      LassoB[[2]] <- tmpLassoB[-1,]
      OracleLassoB[[2]] <- tmpOracleLassoB[-1,]
      LSB[[2]] <- tmpLSB[-1,]
      OracleLSB[[2]] <- tmpOracleLSB[-1,]
      OracleOracleLSB[[2]] <- tmpOracleOracleLSB[-1,]
      
      LassoC[[3]] <- tmpLassoC[-1,]
      OracleLassoC[[3]] <- tmpOracleLassoC[-1,]
      OracleOracleLSC[[3]] <- tmpOracleOracleLSC[-1,]
      
      
      if(T_dim == 15){
        
        col_vec2 <- c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 4))
        col_vec_A <-  c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 5))
        
        
        pdf(file = paste0("./plots/FigS5/","Robustness",  "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pC,".pdf", sep = ""), width = 6, height = 6)
        
        ## Scenario C
        boxplot(data.frame(OracleOracleLSC[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,3]-rep(c_val, N_Sim), 
                           OracleOracleLSC[[3]][,4]-rep(c_val, N_Sim), 
                           OracleOracleLSC[[3]][,6]-rep(c_val, N_Sim), OracleLassoC[[3]][,1]-rep(c_val, N_Sim), OracleLassoC[[3]][,2]-rep(c_val, N_Sim), OracleLassoC[[3]][,3]-rep(c_val, N_Sim), OracleLassoC[[3]][,5], 
                           OracleLassoC[[3]][,(2+300)]-rep(c_val, N_Sim), 
                           OracleLassoC[[3]][,(5+300)], OracleLassoC[[3]][,(2+2*300)]-rep(c_val, N_Sim), OracleLassoC[[3]][,(5+2*300)], 
                           LassoC[[3]][,1]-rep(c_val, N_Sim), LassoC[[3]][,2]-rep(c_val, N_Sim), LassoC[[3]][,3]-rep(c_val, N_Sim), LassoC[[3]][,5], LassoC[[3]][,(2+300)]-rep(c_val, N_Sim), 
                           LassoC[[3]][,(5+300)], LassoC[[3]][,(2+2*300)]-rep(c_val, N_Sim), LassoC[[3]][,(5+2*300)]),
                yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22)), col = col_vec2, cex = 0.2, boxwex = 0.4, lwd = 0.8, whisklty = 1)
        title("p =  901",line=0.3, cex.main = 1)
        axis(2, las = 0, at = c(2, 10, 19), labels = c(TeX(r"(LS-O$^2$)"), "L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
        abline(h = 5, lwd = 1)
        abline(h = 14, lwd = 1)
        abline(v = 0)
        dev.off()
        
        
        ## Scenario A
        pdf(file = paste0("./plots/FigS5/","Robustness", "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pA_vec[3],".pdf", sep = ""), width = 6, height = 6)
        boxplot(data.frame( 
          OracleOracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,4]-rep(c_val, N_Sim), 
          OracleOracleLSA[[3]][,6]-rep(c_val, N_Sim), PesaranCCE[[3]][,1]-rep(c_val, N_Sim),  PesaranCCE[[3]][,2]-rep(c_val, N_Sim),  PesaranCCE[[3]][,3]-rep(c_val, N_Sim), PesaranCCE[[3]][,5], PesaranCCE[[3]][,6]-rep(c_val, N_Sim), 
          PesaranCCE[[3]][,9], PesaranCCE[[3]][,10]-rep(c_val, N_Sim), PesaranCCE[[3]][,13], OracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleLSA[[3]][,2]-rep(c_val, N_Sim),  OracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleLSA[[3]][,5], OracleLSA[[3]][,6]-rep(c_val, N_Sim), 
          OracleLSA[[3]][,9], OracleLSA[[3]][,10]-rep(c_val, N_Sim),OracleLSA[[3]][,13], LSA[[3]][,1]-rep(c_val, N_Sim), LSA[[3]][,2]-rep(c_val, N_Sim), LSA[[3]][,3]-rep(c_val, N_Sim), LSA[[3]][,5], LSA[[3]][,6]-rep(c_val, N_Sim), LSA[[3]][,9], LSA[[3]][,10]-rep(c_val, N_Sim), 
          LSA[[3]][,13], OracleLassoA[[3]][,1]-rep(c_val, N_Sim), OracleLassoA[[3]][,2]-rep(c_val, N_Sim), OracleLassoA[[3]][,3]-rep(c_val, N_Sim), OracleLassoA[[3]][,5], OracleLassoA[[3]][,6]-rep(c_val, N_Sim), OracleLassoA[[3]][,9], OracleLassoA[[3]][,10]-rep(c_val, N_Sim), OracleLassoA[[3]][,13], LassoA[[3]][,1]-rep(c_val, N_Sim), LassoA[[3]][,2]-rep(c_val, N_Sim), LassoA[[3]][,3]-rep(c_val, N_Sim), LassoA[[3]][,5], LassoA[[3]][,6]-rep(c_val, N_Sim), LassoA[[3]][,9], 
          LassoA[[3]][,10]-rep(c_val, N_Sim), LassoA[[3]][,13]),
          
          yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
          col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
        title("p =  13",line=0.3, cex.main = 1)
        axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
        
        abline(h = 5, lwd = 1)
        abline(h = 14, lwd = 1)
        abline(h = 23, lwd = 1)
        abline(h = 32, lwd = 1)
        abline(h = 41, lwd = 1)
        abline(v = 0)
        dev.off()
        
        
        ## Scenario B
        pdf(file = paste0("./plots/FigS5/","Robustness", "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pB_vec[2],".pdf", sep = ""), width = 6, height = 6)
        boxplot(data.frame( OracleOracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,2]-rep(c_val, N_Sim), 
                            OracleOracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,5]-rep(c_val, N_Sim), 
                            OracleOracleLSB[[2]][,10]-rep(c_val, N_Sim), OracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleLSB[[2]][,2]-rep(c_val, N_Sim),  
                            OracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleLSB[[2]][,5], 
                            OracleLSB[[2]][,(2+50)]-rep(c_val,N_Sim),  
                            OracleLSB[[2]][,(5+50)], OracleLSB[[2]][,(2+2*50)]-rep(c_val, N_Sim), OracleLSB[[2]][,(5+2*50)], LSB[[2]][,1]-rep(c_val, N_Sim), LSB[[2]][,2]-rep(c_val, N_Sim), LSB[[2]][,3]-rep(c_val, N_Sim), 
                            LSB[[2]][,5], LSB[[2]][,(2+50)]-rep(c_val, N_Sim), LSB[[2]][,(5+50)], 
                            LSB[[2]][,(2+2*50)]-rep(c_val, N_Sim), 
                            LSB[[2]][,(5+2*50)], 
                            OracleLassoB[[2]][,1]-rep(c_val, N_Sim), OracleLassoB[[2]][,2]-rep(c_val, N_Sim), 
                            OracleLassoB[[2]][,3]-rep(c_val, N_Sim), OracleLassoB[[2]][,5], 
                            OracleLassoB[[2]][,(2+50)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+50)], 
                            OracleLassoB[[2]][,(2+2*50)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+2*50)],  LassoB[[2]][,1]-rep(c_val, N_Sim), LassoB[[2]][,2]-rep(c_val, N_Sim), LassoB[[2]][,3]-rep(c_val, N_Sim), 
                            LassoB[[2]][,5], LassoB[[2]][,(2+50)]-rep(c_val, N_Sim), LassoB[[2]][,(5+50)], 
                            LassoB[[2]][,(2+2*50)]-rep(c_val, N_Sim), LassoB[[2]][,(5+2*50)]), 
                
                yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
                col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
        title( "p = 151" , line=0.3, cex.main = 1)
        axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
        abline(h = 5, lwd = 1)
        abline(h = 14, lwd = 1)
        abline(h = 23, lwd = 1)
        abline(h = 32, lwd = 1)
        abline(v = 0)
        dev.off()
      }
      
      if(T_dim == 50){
        col_vec2 <- c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 4))
        col_vec_A <-  c("green", "blue", "red", "magenta", "bisque", rep(c("green", "blue", "red", "yellow", "magenta", "cyan", "bisque", "darkorange"), 5))
        
        
        ## Scenario A
        # Boxplot for group representatives for each estimator
        pdf(file = paste0("./plots/FigS5/","Robustness",  "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pA_vec[3],".pdf", sep = ""), width = 6, height = 6)
        boxplot(data.frame( 
          OracleOracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleOracleLSA[[3]][,4]-rep(c_val, N_Sim), 
          OracleOracleLSA[[3]][,6]-rep(c_val, N_Sim), PesaranCCE[[3]][,1]-rep(c_val, N_Sim),  PesaranCCE[[3]][,2]-rep(c_val, N_Sim),  PesaranCCE[[3]][,3]-rep(c_val, N_Sim), PesaranCCE[[3]][,5], PesaranCCE[[3]][,(2+15)]-rep(c_val, N_Sim), 
          PesaranCCE[[3]][,(5+15)], PesaranCCE[[3]][,(2+2*15)]-rep(c_val, N_Sim),PesaranCCE[[3]][,(5+2*15)], OracleLSA[[3]][,1]-rep(c_val, N_Sim), OracleLSA[[3]][,2]-rep(c_val, N_Sim),  OracleLSA[[3]][,3]-rep(c_val, N_Sim), OracleLSA[[3]][,5], OracleLSA[[3]][,(2+15)]-rep(c_val, N_Sim), 
          OracleLSA[[3]][,(5+15)], OracleLSA[[3]][,(2+2*15)]-rep(c_val, N_Sim),OracleLSA[[3]][,(5+2*15)], LSA[[3]][,1]-rep(c_val, N_Sim), LSA[[3]][,2]-rep(c_val, N_Sim), LSA[[3]][,3]-rep(c_val, N_Sim), LSA[[3]][,5], LSA[[3]][,(2+15)]-rep(c_val, N_Sim), LSA[[3]][,(5+15)], LSA[[3]][,(2+2*15)]-rep(c_val, N_Sim), 
          LSA[[3]][,(5+2*15)], OracleLassoA[[3]][,1]-rep(c_val, N_Sim), OracleLassoA[[3]][,2]-rep(c_val, N_Sim), OracleLassoA[[3]][,3]-rep(c_val, N_Sim), OracleLassoA[[3]][,15], OracleLassoA[[3]][,(2+15)]-rep(c_val, N_Sim), OracleLassoA[[3]][,(5+15)], OracleLassoA[[3]][,(2+2*15)]-rep(c_val, N_Sim), OracleLassoA[[3]][,(5+2*15)], LassoA[[3]][,1]-rep(c_val, N_Sim), LassoA[[3]][,2]-rep(c_val, N_Sim), LassoA[[3]][,3]-rep(c_val, N_Sim), LassoA[[3]][,5], LassoA[[3]][,(2+15)]-rep(c_val, N_Sim), LassoA[[3]][,(5+15)], 
          LassoA[[3]][,(2+2*15)]-rep(c_val, N_Sim), LassoA[[3]][,(5+2*15)]),
          
          yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40), seq(42, 49)),
          col = col_vec_A, cex = 0.2, boxwex = 0.7, lwd = 0.8, whisklty = 1)
        title("p =  46",line=0.3, cex.main = 1)
        axis(2, las =0, at = c(1, 10, 19, 27, 37, 46), labels = c(TeX(r"(LS-O$^2$)"), "CCE","LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
        
        abline(h = 5, lwd = 1)
        abline(h = 14, lwd = 1)
        abline(h = 23, lwd = 1)
        abline(h = 32, lwd = 1)
        abline(h = 41, lwd = 1)
        abline(v = 0)
        dev.off()
        
        
        ## Scenario B
        pdf(file = paste0("./plots/FigS5/","Robustness", "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pB_vec[2],".pdf", sep = ""), width = 6, height = 6)
        boxplot(data.frame( OracleOracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,2]-rep(c_val, N_Sim), 
                            OracleOracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleOracleLSB[[2]][,5]-rep(c_val, N_Sim), 
                            OracleOracleLSB[[2]][,10]-rep(c_val, N_Sim), OracleLSB[[2]][,1]-rep(c_val, N_Sim), OracleLSB[[2]][,2]-rep(c_val, N_Sim),  
                            OracleLSB[[2]][,3]-rep(c_val, N_Sim), OracleLSB[[2]][,5], 
                            OracleLSB[[2]][,(2+150)]-rep(c_val,N_Sim),  
                            OracleLSB[[2]][,(5+150)], OracleLSB[[2]][,(2+2*150)]-rep(c_val, N_Sim), OracleLSB[[2]][,(5+2*150)], LSB[[2]][,1]-rep(c_val, N_Sim), LSB[[2]][,2]-rep(c_val, N_Sim), LSB[[2]][,3]-rep(c_val, N_Sim), 
                            LSB[[2]][,5], LSB[[2]][,(2+150)]-rep(c_val, N_Sim), LSB[[2]][,(5+150)], 
                            LSB[[2]][,(2+2*150)]-rep(c_val, N_Sim), 
                            LSB[[2]][,(5+2*150)], 
                            OracleLassoB[[2]][,1]-rep(c_val, N_Sim), OracleLassoB[[2]][,2]-rep(c_val, N_Sim), 
                            OracleLassoB[[2]][,3]-rep(c_val, N_Sim), OracleLassoB[[2]][,5], 
                            OracleLassoB[[2]][,(2+150)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+150)], 
                            OracleLassoB[[2]][,(2+2*150)]-rep(c_val, N_Sim), OracleLassoB[[2]][,(5+2*150)],  LassoB[[2]][,1]-rep(c_val, N_Sim), LassoB[[2]][,2]-rep(c_val, N_Sim), LassoB[[2]][,3]-rep(c_val, N_Sim), 
                            LassoB[[2]][,5], LassoB[[2]][,(2+150)]-rep(c_val, N_Sim), LassoB[[2]][,(5+150)], 
                            LassoB[[2]][,(2+2*150)]-rep(c_val, N_Sim), LassoB[[2]][,(5+2*150)]), 
                
                yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22), seq(24,31), seq(33,40)),
                col = col_vec2, cex = 0.2, boxwex = 0.6, lwd = 0.8, whisklty = 1)
        title("p =  451",line=0.3, cex.main = 1)
        axis(2, las = 0, at = c(2, 10, 19, 28, 37), labels = c(TeX(r"(LS-O$^2$)"), "LS-O", "LS","L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
        abline(h = 5, lwd = 1)
        abline(h = 14, lwd = 1)
        abline(h = 23, lwd = 1)
        abline(h = 32, lwd = 1)
        abline(v = 0)
        dev.off()
        
        
        
        
        ## Scenario C
        pdf(file = paste0("./plots/FigS5/","Robustness",  "T_dim", T_dim, "gamma", gamma, "HAC", HAC, "p", pC,".pdf", sep = ""), width = 6, height = 6)
        
        boxplot(data.frame(OracleOracleLSC[[3]][,1]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,2]-rep(c_val, N_Sim), OracleOracleLSC[[3]][,3]-rep(c_val, N_Sim), 
                           OracleOracleLSC[[3]][,4]-rep(c_val, N_Sim), 
                           OracleOracleLSC[[3]][,6]-rep(c_val, N_Sim), OracleLassoC[[3]][,1]-rep(c_val, N_Sim), OracleLassoC[[3]][,2]-rep(c_val, N_Sim), OracleLassoC[[3]][,3]-rep(c_val, N_Sim), OracleLassoC[[3]][,5], 
                           OracleLassoC[[3]][,(2+1000)]-rep(c_val, N_Sim), 
                           OracleLassoC[[3]][,(5+1000)], OracleLassoC[[3]][,(2+2*1000)]-rep(c_val, N_Sim), OracleLassoC[[3]][,(5+2*1000)], 
                           LassoC[[3]][,1]-rep(c_val, N_Sim), LassoC[[3]][,2]-rep(c_val, N_Sim), LassoC[[3]][,3]-rep(c_val, N_Sim), LassoC[[3]][,5], LassoC[[3]][,(2+1000)]-rep(c_val, N_Sim), 
                           LassoC[[3]][,(5+1000)], LassoC[[3]][,(2+2*1000)]-rep(c_val, N_Sim), LassoC[[3]][,(5+2*1000)]),
                yaxt = "n",  horizontal = TRUE, at = c(seq(0,4), seq(6,13), seq(15,22)), col = col_vec2, cex = 0.2, boxwex = 0.4, lwd = 0.8, whisklty = 1)
        title("p =  3001",line=0.3, cex.main = 1)
        axis(2, las = 0, at = c(2, 10, 19), labels = c(TeX(r"(LS-O$^2$)"), "L-O", "L") , tick = FALSE, cex.axis = 0.8, mgp = c(3, 0.25, 1))
        abline(h = 5, lwd = 1)
        abline(h = 14, lwd = 1)
        abline(v = 0)
        dev.off()
      }
    }
  }
  
  
  return(NULL)
  
}



# Obtain Figures S.1, S.2, S.3, S.4 and Tables S.3, S.4, S.5 and S.6
res1 <- generate_results(N_Sim = 1000, T_dim = 15, HAC = 1, gamma = 1)
res2 <- generate_results(N_Sim = 1000, T_dim = 50, HAC = 1, gamma = 1)
# 1401.657 sec elapsed for 10 iterations -> approx 50 hours for 1000 iterations


# Obtain Figure S.5 and Table S.7   
res3 <- generate_results(N_Sim = 1000, T_dim = 15, HAC = 1, gamma = sqrt(10))
res4 <- generate_results(N_Sim = 1000, T_dim = 50, HAC = 1, gamma = sqrt(10))

# Obtain Tables S.8, S.9, S.10
res5 <- generate_results(N_Sim = 1000, T_dim = 15, HAC = 2, gamma = 1)
res6 <- generate_results(N_Sim = 1000, T_dim = 50, HAC = 2, gamma = 1)

# Obtain Tables S.11, S.12, S.13
res7 <- generate_results(N_Sim = 1000, T_dim = 15, HAC = 3, gamma = 1)
res8 <- generate_results(N_Sim = 1000, T_dim = 50, HAC = 3, gamma = 1)


