# Clear memory
rm(list=ls())

# Load packages needed for later
#=============================
library(xtable)

# Functions needed for later
#==========================
source("helper_functions.R") 
# Contains the functions: 
# "round_custom(...)", to prepare estimation results for presentation in tables
# "demean_panel_variables(...)", to calculate the individual specific means
# and get the corresponding demeaned variables

source("hd_panel_estimator.r") 
# Contains the function "hd_panel_estimator(...)" to run the estimation 
# procedure

source("inference_lasso.r")
# Contains the function "inference_lasso(...)" to run the desparsified estimation 
# procedure in order to do inference.

#=========================================================================#
#
# Step 1: Read in data and prepare data for estimation
#
#=========================================================================

# Read in the data file
data <- read.csv("data_application_clean.csv", header = TRUE)

# Make sure the date is ordered by stocks (permno)
data_1722 <- data[order(data$permno),]


# Get the date of the observations which should be the last 5 years (=60 months) 
# of the data set from Apr 2017 - Mar 2022.
month_vec <- row.names(table(data_1722$date))

dim_N <- unique(table(data_1722$date))

# Drop the date variable
data_1722 <-  data_1722[, !(names(data_1722) %in% "date")]

# Store the undemeaned variables
Y <- data_1722$ret
X <- data_1722[,-c(1,2)]
X_varnames <- names(X)


# Construct the stockwise demeaned variables that go into the estimation algorithm
#========================================================================

varnames <- names(data_1722)[-1]
demeaned_data_1722 <- demean_panel_variables(data_1722,id_column="permno",variables = varnames)

# Extract the data needed for the analysis
# Demeaned variables are all in the last p + columns
tmp <- demeaned_data_1722[,(dim(demeaned_data_1722)[2] - length(X_varnames)):(dim(demeaned_data_1722)[2])]
Y_1722 <- tmp$ret_demeaned
# Stockwise demeaned characteristics, which are in the last p columns
X_1722 <- tmp[,!(names(tmp) %in% "ret_demeaned")]
# Assign the characteristics their origianl variable names
names(X_1722) <- X_varnames

# Extract the number of regressors, and panel dimensions ready for estimation
dim_p <- dim(X_1722)[2]
dim_T <- length(month_vec)

# Collect the data ready in the form needed for estimation 
data_design_1722 <- list(y = Y_1722, x = X_1722)


#=====================
#
# Step 2: Data plot and preparations for leave one firm out 
#
#=====================

# Prepare the variables needed for leave one firm out cross-validation
nFolds <- dim_N
fold_vec <- vector("numeric",length=length(Y_1722))
for(i in 1:nFolds)
  fold_vec[((i-1)*dim_T + 1):(i*dim_T)] <- i


# Calculate the "automatic" plot limits
plot_lims_Y <- c(min(Y),max(Y))

# Select every third month
x_ticks <- seq(3,length(month_vec),by=3)
x_axis <- as.Date(month_vec[x_ticks],format="%Y-%m-%d") 


manualcolors<-c('black','gray12','gray33','dimgray','gray59','darkgray'
                ,'lightsteelblue3','slategray','gray','gainsboro', 'royalblue4'
                ,'darkblue','navyblue','mediumblue','slateblue','mediumslateblue'
                ,'cornflowerblue','dodgerblue','royalblue','steelblue1'
                ,'red','red3', 'indianred1', 'indianred3','firebrick4','firebrick2'
                ,'orangered','orange','coral2','darkred')

pdf(file = "Plots/Returns_cols_FE.pdf",width=6,height=4)
par(mar=c(4,4,1,1))
# Adjust the plot limits to get all tick labels showing 
plot(Y[fold_vec==1],ylim=plot_lims_Y ,type="l",xaxt="n",col="black",lwd=2,ylab="Excess returns",xlab="",cex.lab=1.2,cex.axis=1.2)
# Add dates to x-axis 
axis(1,x_ticks,format(x_axis,"%b %y"),las=2,cex.axis=1.2,cex.lab=1.2)
for(firm in 2:dim_N){
  lines(Y[fold_vec==firm],col=manualcolors[firm],lwd=2)
}
dev.off()



#=========================================================================
#
# Step 3: Run the HD CCE estimation procedure  
#         
#
#=========================================================================

hd_cce_1722 <- hd_panel_estimator(data = data_design_1722, obs_N = dim_N, 
                                              obs_T = dim_T, TRUNC = 0.01,
                                              variant = "Lasso", penalty_choice = "NULL", 
                                              NFOLDS = nFolds, fold_vec = fold_vec, 
                                              NLAMBDA = NULL, scree_plot = FALSE) 

# Plot the scree plot and the scree plot for eigenvalues below the truncation
#==========================================
pdf(file = "Plots/Screeplots_a_FE.pdf",width=4,height=3)
par(mar=c(2,4,1,1))
graphics::plot(hd_cce_1722$eigen_values,ylim=c(0,1),xlab="",ylab="",cex.axis=1.2)
title(ylab=expression(hat(psi)[j] / hat(psi)[1]), line=2, cex.lab=1.2)
dev.off()

pdf(file = "Plots/Screeplots_b_FE.pdf",width=4,height=3)
par(mar=c(2,4,1,1))
graphics::plot(hd_cce_1722$eigen_values,ylim=c(0,0.05),xlab="",ylab="",cex.axis=1.2)
abline(h=0.01,lty=2)
title(ylab=expression(hat(psi)[j] / hat(psi)[1]), line=2, cex.lab=1.2)
dev.off()

# From the truncation we choose K = 4 
K <- as.integer(hd_cce_1722$K_hat)

# Extract the penalty parameter
Lambda <- hd_cce_1722$Lambda

# Extract the coefficients
coefs <- hd_cce_1722$coefs

# Construct variable importance measure for K=4
#=================================================

# Get the squared estimated coefficients 
coefs_sq <- (coefs)^2

# Get the variances of the regressors
var_X <- apply(X_1722, 2, var, na.rm = TRUE)
  
# Construct "importance measure impacts" and put in a table with the coefficients
imp_var <- coefs_sq * var_X
tab_imp_var <- cbind(coefs,imp_var)

# Only keep the non-zero coefficients
tab_imp_var <- tab_imp_var[coefs_sq!=0,]

# Normalize to sum to one and add the column names
tab_imp_var[,2] <- tab_imp_var[,2]/sum(tab_imp_var[,2])
colnames(tab_imp_var) <- c("Estimates","Rel. Imp.")



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

# Set the desired Gaussian quantiles 
alpha <- 0.05
c_95 <- c(qnorm(alpha/2),qnorm(1 - alpha/2))

# Confidence intervals for K-hat case
CI_95[,1] <- coefs_despar + c_95[1] * Avar
CI_95[,2] <- coefs_despar + c_95[2] * Avar

# Add the variable names 
var_names <- colnames(X_1722)
rownames(CI_95) <- var_names

# Get the confidence intervals for the Lasso selected variables
CI_95_Lasso <- CI_95[rownames(CI_95) %in% rownames(tab_imp_var),]

# Add the results to one large table
tab_results <- cbind(tab_imp_var,CI_95_Lasso) 

# Sort the results by the relative importance measure
tab_results <- tab_results[order(tab_results[,2],decreasing=TRUE),]

# Adjust the rownames to make them easily modified in latex
rownames(tab_results) <- paste0('\\textit{',rownames(tab_results),'}}')

# Create a latex table
tab_results <- xtable(tab_results,digits=6)
print(tab_results)



