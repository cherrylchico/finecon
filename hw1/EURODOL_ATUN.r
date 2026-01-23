
# Monthly Average Exchange Rate of the Euro (ECU before 1999) against the Dollar

#### First we load the data by indicating the starting year and setting freq = 12 to denote monthly data. We view the loaded data in ts format and verify its attributes.
eurodol <- ts(read.table("eurodol.dat", header = F), start = 1995, freq = 12)
eurodol; str(eurodol)

#### Now, we plot the data and visualize the time series. I've included lines to denote the different years and the average of the time series as well.
plot(eurodol,
     main = "Monthly Average Exchange Rate of the Euro against the Dollar",
     lim = c(0, max(eurodol) + 0.1),
     ylab = "Average Exchange Rate")
abline(v = 1995:2020, col = 4, lty = 3)
abline(h = mean(eurodol), col = 2)

#### We see from the graph immediately a non-monotonic trend due to crashes from the late 90s to early 00s and a spike right before the 2010s. Non-constant variance is immediately noticeable. There are a few indicators of yearly seasonality where peaks are generally experienced in the middle of the year and drop or spike going into the next year.

#### Now we formally assess the variance with the boxplots and mean-variance plots. 

# Boxplot
boxplot(eurodol ~ floor(time(eurodol))) 

# Mean-Variance Plot
mn <- apply(matrix(eurodol,ncol=20),2,mean) 
var <- apply(matrix(eurodol,ncol=20),2,var)
plot(var ~ mn, main="Mean-Variance Plot")
abline(lm(var ~ mn))

#### We can clearly observe non-constant mean. Since we are dealing with financial data, we transform with the log-transform and compare graphs again.

lneurodol <- log(eurodol)
plot(lneurodol,
     main = "Monthly Log Average Exchange Rate of the Euro against the Dollar",
     lim = c(0, max(lneurodol) + 0.1),
     ylab = "Log Average Exchange Rate")
abline(v = 1995:2020, col = 4, lty = 3)
abline(h = mean(lneurodol), col = 2)


par(mfrow = c(2, 2))

# Boxplot
boxplot(eurodol ~ floor(time(eurodol)), main = "Boxplot (Log Transformed)") 
boxplot(lneurodol ~ floor(time(lneurodol)), main = "Boxplot")

# Mean-Variance Plot
lnmn <- apply(matrix(lneurodol,ncol=20),2,mean) 
lnvar <- apply(matrix(lneurodol,ncol=20),2,var)
plot(lnvar ~ lnmn, main="Mean-Variance Plot (Log Transformed)")
abline(lm(lnvar ~ lnmn))

mn <- apply(matrix(eurodol,ncol=20),2,mean) 
var <- apply(matrix(eurodol,ncol=20),2,var)
plot(var ~ mn, main="Mean-Variance Plot")
abline(lm(var ~ mn))

par(mfrow = c(1, 1))

#### Because the variance is not clearly constant, we double check with a box-cox transformation.

library(MASS)

# 1. Run the Box-Cox function
bc <- boxcox(eurodol + 1 ~ 1, plotit = TRUE)

# 2. Extract the exact lambda that maximizes the log-likelihood
optimal_lambda <- bc$x[which.max(bc$y)]
print(optimal_lambda)

#### Since the optimal lambda is 2, we can consider a squared transformation of our variable.

# 3. Transform the data using the exact lambda
bceurodol <- (eurodol^optimal_lambda - 1) / optimal_lambda

par(mfrow = c(2, 2))

# Boxplot
boxplot(bceurodol ~ floor(time(bceurodol)), main = "Boxplot (Square Transformed)") 
boxplot(lneurodol ~ floor(time(lneurodol)), main = "Boxplot (Log Transformed)")

# Mean-Variance Plot
bcmn <- apply(matrix(bceurodol,ncol=20),2,mean) 
bcvar <- apply(matrix(bceurodol,ncol=20),2,var)
plot(bcvar ~ bcmn, main="Mean-Variance Plot (Square Transformed)")
abline(lm(bcvar ~ bcmn))

lnmn <- apply(matrix(lneurodol,ncol=20),2,mean) 
lnvar <- apply(matrix(lneurodol,ncol=20),2,var)
plot(lnvar ~ lnmn, main="Mean-Variance Plot (Log Transformed)")
abline(lm(lnvar ~ lnmn))

par(mfrow = c(1, 1))

#### Comparing the box-cox transformation (squared transformation), we see that log transformation has a more constant variance based on the mean-variance plot. Hence, we go with log transform.

#### Now we formally check for seasonality.

monthplot(lneurodol)
plot(decompose(lneurodol))
acf(lneurodol, ylim = c(-1,1), lag.max = 60, col = c(2,rep(1,11)), lwd=2, main = "ACF Plot")

#### We see from both the monthplot and the seasonal section of the decomposition plots that seasonality clearly exists at a yearly level, hence we adjust for it. The ACF plot reinforces the non-stationarity of the current series as well.

d12lneurodol <- diff(lneurodol, lag = 12)
plot(d12lneurodol,
     main = "Seasonally Adjusted Monthly Log Average Exchange Rate of the Euro against the Dollar",
     lim = c(0, max(d12lneurodol) + 0.1),
     ylab = "Log Average Exchange Rate")
abline(v = 1995:2020, col = 4, lty = 3)
abline(h = mean(d12lneurodol), col = 2)


#### Seasonality is now clearly removed. From the resulting plot, a trend is clearly seen as the mean is not constant. We apply regular differencing until variance increases. For this, we use a function that applies differencing until the variance of the series increases. It returns the previously differenced series when true.

auto_diff <- function(x) {
  # Capture the original name of the object passed to the function
  orig_name <- deparse(substitute(x))
  
  # Initialize with the original data
  best_series <- x
  best_var <- var(best_series, na.rm = TRUE)
  diff_count <- 0
  
  cat("Initial Variance:", best_var, "\n")
  
  repeat {
    # Generate the next difference
    trial_series <- diff(best_series)
    trial_var <- var(trial_series, na.rm = TRUE)
    
    # If the variance increases, stop and keep the 'best_series'
    if (trial_var >= best_var) {
      if (diff_count == 0) {
        cat("Original series is optimal. No differencing applied.\n")
      } else {
        cat("Variance increased to", trial_var, "in additional differencing. Stopping at", diff_count, "difference(s).\n")
      }
      break
    }
    
    # Otherwise, update our best record and keep going
    best_series <- trial_series
    best_var <- trial_var
    diff_count <- diff_count + 1
    
    cat("Difference", diff_count, "applied. New Variance:", best_var, "\n")
    
    # Safety break to prevent infinite loops on weird data
    if (diff_count >= 20) break 
  }
  
  # Logic to rename and assign to Global Environment
  if (diff_count > 0) {
    # Create the name string (e.g., "d1d1eurodol")
    new_name <- paste0(paste(rep("d1", diff_count), collapse = ""), orig_name)
    
    # Assign the final series to the new name in the global environment
    assign(new_name, best_series, pos = 1)
    message("The variable '", new_name, "' has been created in your workspace.")
  }
  
  return(invisible(best_series))
}

# Example Usage:
auto_diff(d12lneurodol)

#### To access the new series, we open out environment and look for a series containing d1. After applying the function, we obtain the series d1d12lneurodol, meaning we applied one difference. 

#### We now see the detrended plot of the time series. We also check the ACF plot for non-stationarity.

par(mfrow = c(2, 1))

plot(d1d12lneurodol,
     main = "Detrended & Seasonally Adjusted Monthly 
             Log Average Exchange Rate of the Euro against the Dollar",
     lim = c(0, max(d1d12lneurodol) + 0.1),
     ylab = "Log Average Exchange Rate")
abline(v = 1995:2020, col = 4, lty = 3)
abline(h = mean(d1d12lneurodol), col = 2)

acf(d1d12lneurodol, ylim = c(-1,1), lag.max = 60, col = c(2,rep(1,11)), lwd=2, main = "ACF Plot")

par(mfrow = c(1, 1))

#### Now we stop as have obtained a stationary time series as seen from both the time series plot and ACF plot.






###########################################################################################





# Number of Individuals Registered as Unemployed at INEM Offices

#### First we load the data by indicating the starting year and setting freq = 12 to denote monthly data. We view the loaded data in ts format and verify its attributes.
atur <- ts(read.table("Atur.dat", header = F), start = 1996, freq = 12)
atur; str(atur)

#### Now, we plot the data and visualize the time series. I've included lines to denote the different years and the average of the time series as well.
plot(atur,
     main = "Monthly Number of Individuals Registered 
             as Unemployed at INEM Offices",
     lim = c(0, max(atur)),
     ylab = "Number of Individuals")
abline(v = 1996:2019, col = 4, lty = 3)
abline(h = mean(atur), col = 2)

#### We see from the graph immediately a trend due to abrupt changes from 2007 to 2010. Non-constant variance is not immediately noticeable. There are a few indicators of yearly seasonality where dips are generally experienced in the middle of the year and drop or spike going into the next year.

#### Now we formally assess the variance with the boxplots and mean-variance plots. 

# Boxplot
boxplot(atur ~ floor(time(atur))) 

# Mean-Variance Plot
mn <- apply(matrix(atur,ncol=20),2,mean) 
var <- apply(matrix(atur,ncol=20),2,var)
plot(var ~ mn, main="Mean-Variance Plot")
abline(lm(var ~ mn))

#### For the most part, we can see that variance is the similar for most years except in the period from 2007 to 2010 where variance jumps (as seen with the long IQRs and outliers in the mean-variance plot). To deal with these outliers possibly due to shocks, we apply log transformation to stabilize the variance and lessen the impact of the outliers.

lnatur <- log(atur)
plot(lnatur,
     main = "Monthly Log Number of Individuals Registered 
             as Unemployed at INEM Offices",
     lim = c(0, max(lnatur) + 0.1),
     ylab = "Log Number of Individuals")
abline(v = 1996:2019, col = 4, lty = 3)
abline(h = mean(lnatur), col = 2)


par(mfrow = c(2, 2))

# Boxplot
boxplot(atur ~ floor(time(atur)), main = "Boxplot (Log Transformed)") 
boxplot(lnatur ~ floor(time(lnatur)), main = "Boxplot")

# Mean-Variance Plot
lnmn <- apply(matrix(lnatur,ncol=20),2,mean) 
lnvar <- apply(matrix(lnatur,ncol=20),2,var)
plot(lnvar ~ lnmn, main="Mean-Variance Plot (Log Transformed)")
abline(lm(lnvar ~ lnmn))

mn <- apply(matrix(atur,ncol=20),2,mean) 
var <- apply(matrix(atur,ncol=20),2,var)
plot(var ~ mn, main="Mean-Variance Plot")
abline(lm(var ~ mn))

par(mfrow = c(1, 1))

#### We now observe, especially in the mean-variance plot, that variance is now more constant even in the presence of outliers.

#### Now we formally check for seasonality.

monthplot(lnatur)
plot(decompose(lnatur))
acf(lneurodol, ylim = c(-1,1), lag.max = 60, col = c(2,rep(1,11)), lwd=2, main = "ACF Plot")

#### We see from both the monthplot and the seasonal section of the decomposition plots that seasonality clearly exists at a yearly level, hence we adjust for it. The ACF plot reinforces the non-stationarity of the current series as well.

d12lnatur <- diff(lnatur, lag = 12)
plot(d12lnatur,
     main = "Seasonally Adjusted Log Monthly Number of Individuals Registered 
             as Unemployed at INEM Offices",
     lim = c(0, max(d12lnatur) + 0.1),
     ylab = "Log Monthly Number of Individuals")
abline(v = 1996:2019, col = 4, lty = 3)
abline(h = mean(d12lnatur), col = 2)


#### Seasonality is now clearly removed. From the resulting plot, a trend is clearly seen as the mean is not constant especially at the shock. We apply regular differencing until variance increases. For this, we use a function that applies differencing until the variance of the series increases. It returns the previously differenced series when true.

auto_diff(d12lnatur)

#### To access the new series, we open out environment and look for a series containing d1. After applying the function, we obtain the series d1d12lnatur, meaning we applied one difference. 

#### We now see the detrended plot of the time series. We also check the ACF plot for non-stationarity.

par(mfrow = c(2, 1))

plot(d1d12lnatur,
     main = "Detrended & Seasonally Adjusted Log Monthly Number of Individuals 
             Registered as Unemployed at INEM Offices",
     lim = c(0, max(d1d12lnatur) + 0.1),
     ylab = "Log Monthly Number of Individuals")
abline(v = 1996:2019, col = 4, lty = 3)
abline(h = mean(d1d12lnatur), col = 2)

acf(d1d12lnatur, ylim = c(-1,1), lag.max = 60, col = c(2,rep(1,11)), lwd=2, main = "ACF Plot")

par(mfrow = c(1, 1))

#### Now we stop as have obtained a stationary time series as seen from both the time series plot and ACF plot.
