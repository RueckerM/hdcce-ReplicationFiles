# Clear memory
rm(list=ls())

# Load packages needed for later
#=============================
library(xtable)

# Functions needed for later
#==========================
source("helper_functions.R") 
# Contains the function "round_custom(...)", to prepare estimation results for 
# presentation in tables

source("hd_panel_estimator.r") 
# Contains the function "hd_panel_estimator(...)" to run the estimation 
# procedure

source("hd_panel_estimator_fixed_lambda.r") # variant of estimation procedure that allows to fix the penalty parameter
# Contains the function "hd_panel_estimator_fixed_pen(...)". This is the same 
# as "hd_panel_estimator(...)" with the extension of allowing a user-specified 
# penalty parameter, the results of which are additionally returned. 
# Careful: Currently this is only added if cross-validation is run, i.e. NFolds  
# (number of folds) and fold_vec (fold indicators) are also user supplied.

source("inference_lasso.r")
# Contains the function "inference_lasso(...)" to run the desparsified estimation 
# procedure in order to do inference.

#=====================#
#
# Step 1: Read in data 
#
#=====================

# Read in the data file
data <- read.csv("data_application_clean.csv", header = TRUE)

# Make sure the date is ordered by stocks (permno)
data_1722 <- data[order(data$permno),]

# Get the date of the observations which should be the last 5 years (=60 months) 
# of the data set from Apr 2017 - Mar 2022.
month_vec <- row.names(table(data_1722$date))

# Obtain the return vector for analysis
Y_1722 <- data_1722$ret
# Calculate the grand mean 
alpha_hat <- mean(Y_1722)
# De-mean the returns
Y_1722 <- Y_1722 - alpha_hat

# Obtain the design matrix for the analysis
X_1722 <- data_1722[, -(1:3)]

# Extract the number of regressors, and panel dimensions ready for estimation
dim_p <- dim(X_1722)[2]
dim_N <- unique(table(data_1722$date))
dim_T <- length(month_vec)

# Collect the data ready in the form needed for estimation 
data_design_1722 <- list(y = Y_1722, x = X_1722)


#=====================
#
# Step 2: Run the estimation 
#
#=====================

# Prepare the variables neede for leave one firm out cross-validation
nFolds <- dim_N
fold_vec <- vector("numeric",length=length(Y_1722))
for(i in 1:nFolds)
  fold_vec[((i-1)*dim_T + 1):(i*dim_T)] <- i

#=========================================================================
#
# Step 3a: Estimated the number of factors. 
#          Do one run to get the estimated eigenvalues and the scree plots
#          in order to estimate the number of factors
#
#=========================================================================

hd_cce_1722 <- hd_panel_estimator(data = data_design_1722, obs_N = dim_N, 
                                              obs_T = dim_T, TRUNC = 0.01,
                                              variant = "Lasso", penalty_choice = "NULL", 
                                              NFOLDS = nFolds, fold_vec = fold_vec, 
                                              NLAMBDA = NULL, scree_plot = FALSE) 

# Plot the scree plot and the scree plot for eigenvalues below the truncation
#==========================================
pdf(file = "Plots/Screeplots_a.pdf",width=4,height=3)
par(mar=c(2,4,1,1))
graphics::plot(hd_cce_1722$eigen_values,ylim=c(0,1),xlab="",ylab="",cex.axis=1.2)
title(ylab=expression(hat(psi)[j] / hat(psi)[1]), line=2, cex.lab=1.2)
dev.off()

pdf(file = "Plots/Screeplots_b.pdf",width=4,height=3)
par(mar=c(2,4,1,1))
graphics::plot(hd_cce_1722$eigen_values[-(1:2)],ylim=c(0,0.01),xlab="",ylab="",cex.axis=1.2)
title(ylab=expression(hat(psi)[j] / hat(psi)[1]), line=2, cex.lab=1.2)
dev.off()


# From the scree plots, there are K=2 clearly separated eigenvalues. 
# Zooming in there are an additional K=5 clearly separated eigenvalues. In the 
# paper, focus on the "cautious" choice of K=5.

# Set the number of factors equal to 5
#======================================
K <- as.integer(5)

# Retrieve the cross-validated penalty choice for K=5
Lambda <- hd_panel_estimator(data = data_design_1722, obs_N = dim_N, 
                             obs_T = dim_T, TRUNC = 0.01, NFACTORS = K,
                             variant = "Lasso", penalty_choice = "NULL", 
                             NFOLDS = nFolds, fold_vec = fold_vec, 
                             NLAMBDA = NULL, scree_plot = FALSE)$Lambda


#=========================================================================
#
# Step 3b: For this setting, see how the estimates change as we vary the 
#          number of factors from 1 to 6. 
#
#=========================================================================

# Range of factors for which to estimate the model
FACTORS <- 1:6

# Matrix to store the coefficients
coefs <- matrix(NA,nrow=dim(X_1722)[2],ncol=length(FACTORS)) 

for(Facs in FACTORS){
  NFacs <- as.integer(Facs)
  coefs[,Facs] <- hd_panel_estimator_fixed_pen(data = data_design_1722, obs_N = dim_N, 
                                                             obs_T = dim_T, TRUNC = 0.01, NFACTORS = NFacs,
                                                             variant = "Lasso", penalty_choice = "NULL", 
                                                             NFOLDS = nFolds, fold_vec = fold_vec, 
                                                             NLAMBDA = NULL, scree_plot = FALSE,
                                                             pen=Lambda)$coefs_pen
  }


# Construct the table of the results given in the appendix
#=================================
coefs_tab <- coefs
# Add the rownames
rownames(coefs_tab) <- colnames(X_1722)
# Only keep nonzero rows 
coefs_tab <- coefs_tab[(rowSums(abs(coefs_tab)) != 0),]
# Round the results for presentation
tab_results_K1_6 <- round_custom(coefs_tab)
# Add the penalty parameter
tab_results_K1_6   <- rbind(rep(Lambda,length(FACTORS)),tab_results_K1_6)
# Add the rownames so that they can be "easily" italized in Latex afterwards
var_names <- rownames(coefs_tab)
var_names <- paste0('\\textit{',var_names,'}}')
rownames(tab_results_K1_6) <- c("Lambda",var_names)

# Create a latex table
tab_results_K1_6 <- xtable(tab_results_K1_6)
print(tab_results_K1_6)


#=========================================================================
#
# Step 3c: Plot the estimated expected returns
#
#=========================================================================

# Calculate the estimated expected returns
exp_ret <- rep(alpha_hat,length(Y_1722)) + as.matrix(X_1722) %*% as.vector(coefs[,K])

# Calculate the "automatic" plot limits
plot_lims_Y <- c(min(Y_1722),max(Y_1722))
plot_lims_ret <- c(min(exp_ret),max(exp_ret)) 

# Plot the actual returns 
#============================================

# Select every third month
x_ticks <- seq(3,length(month_vec),by=3)
x_axis <- as.Date(month_vec[x_ticks],format="%Y-%m-%d") 


pdf(file = "Plots/Returns.pdf",width=6,height=4)
par(mar=c(4,4,1,1))
# Adjust the plot limits to get all tick labels showing 
plot(Y_1722[fold_vec==1],ylim=plot_lims_Y ,type="l",xaxt="n",col="darkgray",lwd=2,ylab="Excess returns",xlab="",cex.lab=1.2,cex.axis=1.2)
# Add dates to x-axis 
axis(1,x_ticks,format(x_axis,"%b %y"),las=2,cex.axis=1.2,cex.lab=1.2)
for(firm in 2:dim_N){
  lines(Y_1722[fold_vec==firm],col="darkgray",lwd=2)
}
dev.off()



# Plot the estimated expected returns (for K=5)
#============================================

pdf(file = "Plots/Exp_Ret_K5_Lasso.pdf",width=6,height=4)
par(mar=c(4,4,1,1))
plot(exp_ret[fold_vec==1],ylim=plot_lims_ret,type="l",xaxt="n",col="darkgray",lwd=2,ylab="Est. expected return",xlab="",cex.axis=1.2,cex.lab=1.2)
# Add dates to x-axis 
axis(1,x_ticks,format(x_axis,"%b %y"),las=2,cex.axis=1.2,cex.lab=1.2)
for(firm in 2:dim_N){
  lines(exp_ret[fold_vec==firm],col="darkgray",lwd=2)
}
dev.off()


#===============================================================================
#
# Step 4: Inference with debiased HD-CCE estimator: 
#                  Test whether beta_j = 0 or not
#
#===============================================================================

# Get desparsified estimates for 
#========================================================================

# Storage for the de-sparsified estimate 
coefs_despar <- vector(mode="double",length=dim_p) 

# Storage for the variance estimate 
Avar <- vector(mode="double",length=dim_p)  

# Get the de-sparsified estimates for inference
for(j in 1:dim_p){
  tmp <- inference_lasso(data = data_design_1722, obs_N = dim_N, nFactors = K, 
                         obs_T = dim_T, nFolds = nFolds, foldid = fold_vec, 
                         coef_index = j)
  coefs_despar[j] <- tmp$coef_despar
  Avar[j] <- tmp$Avar
}



# Calculate the Confidence Intervals fo the de-sparsified estimates
# pointwise
CI_95 <- matrix(NA,nrow = dim_p, ncol = 2)  
# with Bonferroni correction
CI_95_bc <- matrix(NA,nrow = dim_p, ncol = 2)  

# Set the desired Gaussian quantiles and for the bonferonni corrected procedure
alpha <- 0.05
c_95 <- c(qnorm(alpha/2),qnorm(1 - alpha/2))
bc_95 <- c(qnorm(alpha/dim_p/2),qnorm(1 - alpha/dim_p/2))

# Confidence intervals for K=5 case
CI_95[,1] <- coefs_despar + c_95[1] * Avar
CI_95[,2] <- coefs_despar + c_95[2] * Avar

CI_95_bc[,1] <- coefs_despar + bc_95[1] * Avar
CI_95_bc[,2] <- coefs_despar + bc_95[2] * Avar

# Add the variable names 
var_names <- colnames(X_1722)
rownames(CI_95) <- var_names
rownames(CI_95_bc) <- var_names

# Find the "significant" entries where the CIs do not cover 0.
sig_ind <- rep(0,dim_p) 
sig_ind_bc <- rep(0,dim_p)

for(j in 1:dim_p){
  sig_ind[j] <- (CI_95[j,1] > 0 | CI_95[j,2] < 0) 
  sig_ind_bc[j] <- (CI_95_bc[j,1] > 0 | CI_95_bc[j,2] < 0) 
}

# Store the "significant" entries 
CI_95_sig <- CI_95[(sig_ind==1),]
CI_95_bc_sig <- CI_95_bc[(sig_ind_bc==1),]

# Get the results for the paper
tab_CI_95_sig <- CI_95_sig
tab_CI_95_bc_sig <- CI_95_bc_sig

# Change the variable names to be able to easily adjust them in latex
# Pointwise CIs
var_names <- rownames(tab_CI_95_sig)
rownames(tab_CI_95_sig) <- paste0('\\textit{',var_names,'}}')

tab_CI_95_sig <- xtable(tab_CI_95_sig,digits=6)
print(tab_CI_95_sig)

# Bonferroni CIs
var_names <- rownames(tab_CI_95_bc_sig)
rownames(tab_CI_95_bc_sig) <- paste0('\\textit{',var_names,'}}')

tab_CI_95_bc_sig <- xtable(tab_CI_95_bc_sig,digits=6)
print(tab_CI_95_bc_sig)




#==============================================================================
#
#    Additional results using K=2 that are reported in the appendix
#
#=============================================================================

# Set the number of factors equal to 5
#======================================
K <- as.integer(2)

# Retrieve the cross-validated penalty choice for K=5
Lambda <- hd_panel_estimator(data = data_design_1722, obs_N = dim_N, 
                             obs_T = dim_T, TRUNC = 0.01, NFACTORS = K,
                             variant = "Lasso", penalty_choice = "NULL", 
                             NFOLDS = nFolds, fold_vec = fold_vec, 
                             NLAMBDA = NULL, scree_plot = FALSE)$Lambda


#=========================================================================
#
# Step 3b: For this setting, see how the estimates change as we vary the 
#          number of factors from 1 to 6. 
#
#=========================================================================

# Range of factors for which to estimate the model
FACTORS <- 1:6

# Matrix to store the coefficients
coefs <- matrix(NA,nrow=dim(X_1722)[2],ncol=length(FACTORS)) 

for(Facs in FACTORS){
  NFacs <- as.integer(Facs)
  coefs[,Facs] <- hd_panel_estimator_fixed_pen(data = data_design_1722, obs_N = dim_N, 
                                               obs_T = dim_T, TRUNC = 0.01, NFACTORS = NFacs,
                                               variant = "Lasso", penalty_choice = "NULL", 
                                               NFOLDS = nFolds, fold_vec = fold_vec, 
                                               NLAMBDA = NULL, scree_plot = FALSE,
                                               pen=Lambda)$coefs_pen
}


# Construct the table of the results given in the appendix
#=================================
coefs_tab <- coefs
# Add the rownames
rownames(coefs_tab) <- colnames(X_1722)
# Only keep nonzero rows 
coefs_tab <- coefs_tab[(rowSums(abs(coefs_tab)) != 0),]
# Round the results for presentation
tab_results_K1_6 <- round_custom(coefs_tab)
# Add the penalty parameter
tab_results_K1_6   <- rbind(rep(Lambda,length(FACTORS)),tab_results_K1_6)
# Add the rownames so that they can be "easily" italized in Latex afterwards
var_names <- rownames(coefs_tab)
var_names <- paste0('\\textit{',var_names,'}}')
rownames(tab_results_K1_6) <- c("Lambda",var_names)

# Create a latex table
tab_results_K1_6 <- xtable(tab_results_K1_6)
print(tab_results_K1_6)


#=========================================================================
#
# Step 3c: Plot the estimated expected returns
#
#=========================================================================

# Calculate the estimated expected returns
exp_ret <- rep(alpha_hat,length(Y_1722)) + as.matrix(X_1722) %*% as.vector(coefs[,K])

# Calculate the "automatic" plot limits
plot_lims_ret <- c(min(exp_ret),max(exp_ret)) 

# Plot the actual returns 
#============================================

# Select every third month
x_ticks <- seq(3,length(month_vec),by=3)
x_axis <- as.Date(month_vec[x_ticks],format="%Y-%m-%d") 

# Plot the estimated expected returns (for K=2)
#============================================

pdf(file = "Plots/Exp_Ret_K2_Lasso.pdf",width=6,height=4)
par(mar=c(4,4,1,1))
plot(exp_ret[fold_vec==1],ylim=plot_lims_ret,type="l",xaxt="n",col="darkgray",lwd=2,ylab="Est. expected return",xlab="",cex.axis=1.2,cex.lab=1.2)
# Add dates to x-axis 
axis(1,x_ticks,format(x_axis,"%b %y"),las=2,cex.axis=1.2,cex.lab=1.2)
for(firm in 2:dim_N){
  lines(exp_ret[fold_vec==firm],col="darkgray",lwd=2)
}
dev.off()


#===============================================================================
#
# Step 4: Inference with debiased HD-CCE estimator: 
#                  Test whether beta_j = 0 or not
#
#===============================================================================

# Get desparsified estimates for 
#========================================================================

# Storage for the de-sparsified estimate 
coefs_despar <- vector(mode="double",length=dim_p) 

# Storage for the variance estimate 
Avar <- vector(mode="double",length=dim_p)  

# Get the de-sparsified estimates for inference
for(j in 1:dim_p){
  tmp <- inference_lasso(data = data_design_1722, obs_N = dim_N, nFactors = K, 
                         obs_T = dim_T, nFolds = nFolds, foldid = fold_vec, 
                         coef_index = j)
  coefs_despar[j] <- tmp$coef_despar
  Avar[j] <- tmp$Avar
}



# Calculate the Confidence Intervals fo the de-sparsified estimates
# pointwise
CI_95 <- matrix(NA,nrow = dim_p, ncol = 2)  
# with Bonferroni correction
CI_95_bc <- matrix(NA,nrow = dim_p, ncol = 2)  

# Set the desired Gaussian quantiles and for the bonferonni corrected procedure
alpha <- 0.05
c_95 <- c(qnorm(alpha/2),qnorm(1 - alpha/2))
bc_95 <- c(qnorm(alpha/dim_p/2),qnorm(1 - alpha/dim_p/2))

# Confidence intervals for K=5 case
CI_95[,1] <- coefs_despar + c_95[1] * Avar
CI_95[,2] <- coefs_despar + c_95[2] * Avar

CI_95_bc[,1] <- coefs_despar + bc_95[1] * Avar
CI_95_bc[,2] <- coefs_despar + bc_95[2] * Avar

# Add the variable names 
var_names <- colnames(X_1722)
rownames(CI_95) <- var_names
rownames(CI_95_bc) <- var_names

# Find the "significant" entries where the CIs do not cover 0.
sig_ind <- rep(0,dim_p) 
sig_ind_bc <- rep(0,dim_p)

for(j in 1:dim_p){
  sig_ind[j] <- (CI_95[j,1] > 0 | CI_95[j,2] < 0) 
  sig_ind_bc[j] <- (CI_95_bc[j,1] > 0 | CI_95_bc[j,2] < 0) 
}

# Store the "significant" entries 
CI_95_sig <- CI_95[(sig_ind==1),]
CI_95_bc_sig <- CI_95_bc[(sig_ind_bc==1),]

# Get the results for the paper
tab_CI_95_sig <- CI_95_sig
tab_CI_95_bc_sig <- CI_95_bc_sig

# Change the variable names to be able to easily adjust them in latex
# Pointwise CIs
var_names <- rownames(tab_CI_95_sig)
rownames(tab_CI_95_sig) <- paste0('\\textit{',var_names,'}}')

tab_CI_95_sig <- xtable(tab_CI_95_sig,digits=6)
print(tab_CI_95_sig)

# Bonferroni CIs
var_names <- rownames(tab_CI_95_bc_sig)
rownames(tab_CI_95_bc_sig) <- paste0('\\textit{',var_names,'}}')

tab_CI_95_bc_sig <- xtable(tab_CI_95_bc_sig,digits=6)
print(tab_CI_95_bc_sig)
















