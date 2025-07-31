# Empirical Application for "Estimation and Inference in High-Dimensional Panel
# Data Models with Interactive Fixed Effects" (2025) by Maximilian Rücker, 
# Michael Vogt, Oliver Linton and Christopher Walsh.


########################### Install dependencies ###############################

# install.packages("glmnet")
# install.packages("matrixStats")
# install.packages("this.path")
setwd(this.path::here())

library(matrixStats)
library(glmnet)

########################### Read in data #######################################
data <- read.csv("./data/imputed_data.csv", header = TRUE)






########################### Data preparation ###################################
T_dim <- 60
p <- 90
n_dim <- 29

month_vec <- data$date
T_window <- unique(month_vec)[((507-(T_dim-1)):507)] # 507 consecutive months  
data_1722 <- data[data$date %in% T_window, ]

# Remove permnos 10333 and 24942 and select firm characteristics (column 7 until 108)
data_1722 <- data_1722[order(data_1722$permno),]
X_1722 <- data_1722[((data_1722$permno!= 10333) & (data_1722$permno!= 24942)),
                    7:108]

# Remove certain characteristics
excl_var_names <- c("ipo", "sin", "mom1m", "mom6m", "mom12m", "mom36m", 
                    "chmom", "ear", "maxret", "retvol", "ill",  "divi")



# Sort the columns alphabetically 
X_1722 <- X_1722[,order(names(X_1722))]
excl_var <- numeric(length(excl_var_names))

count <- 1
for(l in excl_var_names){
  excl_var[count] <- grep(l, colnames(X_1722)) 
  count <- count + 1
}


X_1722 <- X_1722[, -excl_var]

# Re-scale for alignment
X_1722$mve_ia <- X_1722$mve_ia*(10^(-4))
X_1722$zerotrade <- X_1722$zerotrade*(10^5)

Y_1722 <- data_1722$ret[((data_1722$permno!= 10333) & (data_1722$permno!= 24942))]

# Get the characteristics variance for better comparison
X_1722_var <- apply(X_1722, MARGIN = 2, FUN = var)


############################## Computations ###################################

###################### 1. Projection matrix ###################################
# Cross-sectional averages of the regressors

X_bar <- matrix(NA, ncol = p, nrow = T_dim)

for(t in 1:T_dim){
  indices <- seq(t,n_dim * T_dim, by = T_dim)
  X_bar[t,] <- colMeans(X_1722[indices,])
}

# Empirical covariance matrix and eigenstructure
Cov_X_bar <- (1/T_dim) * t(X_bar) %*% X_bar
Cov_X_bar_eigen <- eigen(Cov_X_bar, symmetric = TRUE)


#============================================================================#
#  Estimation of number of factors
#============================================================================#

# Normalize the eigenvalues with the largest one
eigen_values <- Cov_X_bar_eigen$values/Cov_X_bar_eigen$values[1]



# Estimated Factors
K_hat <- sum((0.01< eigen_values)) # K = 2 Factors estimated

W_hat_tmp <- X_bar %*% Cov_X_bar_eigen$vectors[,1:K_hat]

W_hat <- cbind((rep(1,(T_dim-1) )), W_hat_tmp[-T_dim,], W_hat_tmp[-1,])

Pi_hat <- diag((T_dim-1)) - W_hat %*% solve(t(W_hat) %*% W_hat)  %*% t(W_hat)

#============================================================================#
#  Project the data
#============================================================================#


Y_hat <- rep(NA, (T_dim-1) * n_dim)
X_hat <- matrix(NA, nrow = n_dim * (T_dim-1), ncol = p)


for(i in 1:n_dim){
  
  index1 <- ((i-1) * (T_dim-1) + 1):(i * (T_dim-1))
  index2 <- ((i-1) * T_dim + 1):(i * T_dim -1)
  
  Y_hat[index1] <- Pi_hat %*% t(t(Y_1722[index2]))
  X_hat[index1,] <- Pi_hat %*% t(t(X_1722[index2,]))
  
}



######################## 2. Lasso Estimation ###################################
# Leave one firm out CV
foldid <- rep(1:n_dim, each = (T_dim-1))

res <-  coef(
  glmnet::cv.glmnet(
    x = X_hat,
    y = Y_hat,
    foldid = foldid,
    standardize = TRUE,
    intercept = FALSE
  ),
  s = "lambda.min"
)[-1]

# Report the scaled (with characteristics standard deviation) coefficients
res_scaled <- res*sqrt(X_1722_var)
names(res_scaled) <- colnames(X_1722)


non_zero_coefs <- names(res_scaled[res_scaled!=0])
non_zero_coefs_index <- numeric(length(non_zero_coefs))

nbr.nonzerocoefs <- length(non_zero_coefs)
count <- 1
for (coef in non_zero_coefs) {
  non_zero_coefs_index[count] <- grep(paste("^",coef,"$", sep = ""), 
                                      colnames(X_1722))
  count <- count + 1
}
non_zero_coefs_val <- res_scaled[res_scaled!=0]

order_vec <- order(abs(res_scaled))



########################### 3. Inference ######################################

# Varianz estimator
significance_codes <- numeric(p)

# Coefficient for inference
alpha <- c(0.01, 0.05, 0.1)

#============================================================================#
#  Projection Matrix TILDE_Pi
#============================================================================#

# This line is for the Figure 8.2 
pdf(file = "./plots/CIandEstimates95.pdf", width = 7, height = 6)

par(mar=c(1.5,5,1,5), oma = c(0,0,0,0))
plot(res_scaled[order_vec], 1:p, cex = 0.75, xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
     pch = 18, xlim = c(-0.1, 0.1), tck = 0.02, type = "n")

y_tix <- 1:p
y_lab_boolean <- (y_tix %% 2 == 0)

for (i in 1:(p-1)) {
  y_bottom <- y_tix[i]
  y_top <- y_tix[i + 1]
  
  rect(xleft = -0.12, xright = 0.12, ybottom = y_bottom, ytop = y_top,
       col = ifelse(i %% 2 == 0, "grey90", "white"), border = NA)
}
rect(xleft = -0.12, xright = 0.12, ybottom = p, ytop = (p+1),
     col =  "grey90", border = NA)


CI_height <- 1 # This is for the CI height 

for (COEF_INDEX in order_vec) {
  
  # Empirical covariance matrix and eigenstructure
  Cov_X_bar_tilde <- (1/T_dim) * t(X_bar[,-COEF_INDEX]) %*% X_bar[,-COEF_INDEX]
  Cov_X_bar_tilde_eigen <- eigen(Cov_X_bar_tilde, symmetric = TRUE)
  
  
  W_tilde_tmp <- X_bar[,-COEF_INDEX] %*% Cov_X_bar_tilde_eigen$vectors[,1:K_hat]
  
  W_tilde <- cbind((rep(1,(T_dim-1) )), W_tilde_tmp[-1,], W_tilde_tmp[-T_dim,])
  
  Pi_tilde <- diag((T_dim-1)) -  W_tilde %*% solve(t(W_tilde) %*% W_tilde)  %*% t(W_tilde)
  
  #============================================================================#
  #  Project the data
  #============================================================================#
  
  
  Y_tilde <- rep(NA, (T_dim-1) * n_dim)
  X_tilde <- matrix(NA, nrow = n_dim * (T_dim-1), ncol = p)
  
  
  for(i in 1:n_dim){
    index1 <- ((i-1) * (T_dim-1) + 1):(i * (T_dim-1))
    index2 <- ((i-1) * T_dim + 1):(i * T_dim -1)
    
    Y_tilde[index1] <- Pi_tilde %*% t(t(Y_1722[index2]))
    X_tilde[index1,] <- Pi_tilde %*% t(t(X_1722[index2,]))
  }
  
  
  
  
  fit_Lasso <- glmnet::cv.glmnet(x = X_tilde, y = Y_tilde,
                                 foldid = foldid, family = "gaussian",
                                 alpha = 1, intercept = FALSE)
  
  coefs_Lasso <- stats::coef(fit_Lasso, s = "lambda.min")[-1]
  lambda_cv  <- fit_Lasso$lambda.min
  
  
  # Prediction Main regression and residuals
  yhat_Lasso <- stats::predict(fit_Lasso, newx = X_tilde,
                               type = "response", s = "lambda.min")
  resid_Lasso <- Y_tilde - yhat_Lasso
  
  
  
  fit_node_Lasso <- glmnet::cv.glmnet(x = X_tilde[,-COEF_INDEX], y = X_tilde[,COEF_INDEX],
                                      foldid = foldid, family = "gaussian",
                                      intercept = FALSE)
  
  coefs_node_Lasso <- stats::coef(fit_node_Lasso, s = "lambda.min")[-1]
  kappa_cv  <- fit_node_Lasso$lambda.min
  
  
  
  
  # Kappa grid that cv.glmnet uses for cross-validation
  kappa_grid <- fit_node_Lasso$lambda
  
  # Index for the CV chosen kappa
  kappa_cv_idx <- fit_node_Lasso$index[1]
  
  #============================================================================#
  # Calculate the de-sparsified estimator
  #============================================================================#
  
  #============================================================================#
  # Choose the nodewise penalty parameter
  #============================================================================#
  
  # Kappa grid length
  kappa_grid_len <- length(kappa_grid)
  
  # Calculate the "scaled asymptotic variance"
  var_scaled <- rep(0, len = kappa_grid_len)
  
  
  # Heteroscedastic robust estimator
  for(k in 1:kappa_grid_len){
    yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                                      type = "response", s = kappa_grid[k])
    resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso
    
    tmp <- numeric(n_dim)
    for (i in 1:n_dim) {
      sigma_eps_estimate <- (1/(T_dim-1))*((T_dim-1)/((T_dim-1)-K_hat))*sum(resid_Lasso[((i-1)*(T_dim-1)+1):(i*(T_dim-1))]^2)
      
      tmp[i] <- sigma_eps_estimate*sum(resid_node_Lasso[((i-1)*(T_dim-1)+1):(i*(T_dim-1))]^2)
    }
    var_scaled[k] <-  (1/(t(X_tilde[,COEF_INDEX]) %*% resid_node_Lasso)^2)*sum(tmp)
  }
  
  
  
  # 25% increase of the scaled variance for truncation
  V_TRUNC <- 1.25 * var_scaled[kappa_cv_idx]
  
  # Indicator whether the scaled variance is less than the threshold
  kappa_idx <- 1
  for (l in 1:kappa_grid_len) {
    if(var_scaled[l] <= V_TRUNC){
      kappa_idx <- l
    }
    if(var_scaled[l] > V_TRUNC){
      break
    }
  }
  
  
  
  #============================================================================#
  # Construction of the de-sparsified estimator
  #============================================================================#
  # Nodewise lasso residuals
  yhat_node_Lasso <- stats::predict(fit_node_Lasso, newx = X_tilde[,-COEF_INDEX],
                                    type = "response", s = kappa_grid[kappa_idx])
  resid_node_Lasso <- X_tilde[, COEF_INDEX] - yhat_node_Lasso
  
  despar_beta <- coefs_Lasso[COEF_INDEX] + t(resid_node_Lasso) %*% resid_Lasso / t(resid_node_Lasso) %*% X_tilde[, COEF_INDEX]
  
  # Collect results
  Avar <- sqrt(var_scaled[kappa_idx])
  
  # We use scaled (with characteristics standard deviation) confidence bands 
  conf_band_min <- sqrt(X_1722_var[COEF_INDEX])*(rep(despar_beta, length(alpha)) + Avar * qnorm(alpha/2))
  conf_band_max <- sqrt(X_1722_var[COEF_INDEX])*(rep(despar_beta, length(alpha)) + Avar * qnorm(1-alpha/2))
  
  confidence_band <- cbind(conf_band_min, conf_band_max)
  
  # This line is for the plot 
  
  upper <- (p-nbr.nonzerocoefs-1)
  if(CI_height <= upper){
    lines(c(conf_band_min[2], conf_band_max[2]), c(upper-(CI_height-1) + 0.5,
                                                   upper-(CI_height-1) + 0.5),
          col = "black")
  }else{
    lines(c(conf_band_min[2], conf_band_max[2]), c(CI_height+0.5,CI_height+0.5),
          col = "black")
  }
  CI_height <- CI_height + 1
  
  
  ########################## Generate nice x-axis labels ########################
  for(a in 1:3){
    if(conf_band_min[a] < 0 & 0 < conf_band_max[a])
      significance_codes[COEF_INDEX] <- a
  }
  
}######################## End iteration COEF INDEX ########################
significance_stars <- numeric(length(significance_codes))
for (j in 1:p) {
  significance_stars[j] <- paste(rep("*", (3-significance_codes[j])), collapse = "")
}

ylabel_complete <- paste(names(res_scaled), significance_stars, sep = "") 
ylabel_complete_sorted <- ylabel_complete[order(abs(res_scaled))]
ylabel_final <- c(rev(ylabel_complete_sorted[1:upper]),
                  ylabel_complete_sorted[(upper+1):p])

abline(v = 0, cex = 0.5)

points(res_scaled[order_vec], (1:p + 0.5), pch = 18, cex = 0.8)
axis(4, at = (y_tix[!(y_tix%%2 == 0)]+0.5), 
     labels = ylabel_final[!y_lab_boolean], cex.axis = 0.7, las = 1,
     mgp = c(1,0.1,0), tck = 0.02)
axis(2, at = (y_tix[(y_tix%%2 == 0)]+0.5), labels = ylabel_final[y_lab_boolean],
     cex.axis = 0.7, las = 1,
     mgp = c(1,0.1,0), tck = 0.02)
axis(1, at = seq(-0.1, 0.1, by = 0.01), cex.axis = 0.7,
     las = 1, mgp = c(1,0.1,0) , tck = 0.02)
############################ axis end #########################################

dev.off()



# Generate significance stars for non-selected coefficients
non_selected_significant_variables <- (!(significance_stars == "") & res_scaled == 0)
non_selected_significant_codes <- significance_stars[non_selected_significant_variables]
names(non_selected_significant_codes) <- names(res_scaled)[(!(significance_stars == "") & res_scaled == 0)]



# Generate Figure 8.1
ylabel <- paste(non_zero_coefs, significance_stars[non_zero_coefs_index], sep = "")
sorted_names_non_zero_coefs <- names(sort(abs(non_zero_coefs_val)))
order_names_non_zero_coefs <- order((abs(non_zero_coefs_val)))
sorted_ylabel <- ylabel[order_names_non_zero_coefs]

# Colors for sign of coefficients
color_vec_sign <- rep("steelblue1", length(non_zero_coefs))
count <- 1
for(j in order_names_non_zero_coefs){
  if(non_zero_coefs_val[j] > 0){
    color_vec_sign[count] <- "coral"
  }
  count <- count + 1
}


pdf(file = "./plots/VarianceImportanceMeasure.pdf", width = 6, height = 4)
par(mar=c(1.5,5.5,1,1), oma = c(0,0,0,0))
barplot(sort(abs(non_zero_coefs_val)), horiz = TRUE, names = sorted_ylabel, las = 1,
        cex.names = 0.7, cex.axis = 0.7, col = color_vec_sign, mgp = c(1,0,0), tck = 0.02)
dev.off()






