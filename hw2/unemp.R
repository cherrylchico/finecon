library(forecast)
library(tseries)
library(lmtest)
library(nortest)
library(MASS)

# -- INDIVIDUAL ANALYSIS FUNCTION --

analyze_sarima <- function(series, order, seasonal) {
  
  # --- Setup ---
  s <- frequency(series)
  n <- length(series)
  model_name <- paste0("ARIMA(", paste(order, collapse=","), ")(", 
                       paste(seasonal$order, collapse=","), ")", s)
  
  cat("\n==================================================================\n")
  cat("  ANALYSIS FOR MODEL:", model_name, "\n")
  cat("==================================================================\n")
  
  # -----------------------------------------------------------------------
  # 1. Estimation & Significance Check
  # -----------------------------------------------------------------------
  fit_full <- arima(series, order = order, seasonal = seasonal)
  
  diag_var <- diag(fit_full$var.coef)
  if (any(diag_var <= 0)) {
    cat("\n⚠️  CRITICAL WARNING: The Hessian matrix is unstable (Negative Variances).\n")
    se_vec <- ifelse(diag_var > 0, sqrt(diag_var), NA)
  } else {
    se_vec <- sqrt(diag_var)
  }
  
  t_vec  <- fit_full$coef / se_vec
  p_vec  <- 2 * (1 - pnorm(abs(t_vec)))
  
  coeftest <- data.frame(
    Estimate = fit_full$coef,
    Std.Error = se_vec,
    t_stat = t_vec,
    p_value = p_vec
  )
  
  coeftest$Significant <- ifelse(!is.na(coeftest$p_value) & 
                                   (coeftest$p_value < 0.05 | abs(coeftest$t_stat) > 2), 
                                 "Yes", "No")
  
  cat("\n[1] COEFFICIENT SIGNIFICANCE (Alpha = 0.05 or |t| > 2):\n")
  print(round(coeftest[,1:4], 4))
  cat("\nSignificance Flags:\n")
  print(coeftest[,5, drop=FALSE])
  
  # -----------------------------------------------------------------------
  # 2. Residual Analysis (Diagnostics)
  # -----------------------------------------------------------------------
  resid <- residuals(fit_full)
  
  # A. Causality / Invertibility Roots (Min Modulus)
  # AR Roots
  ar_poly <- fit_full$model$phi
  if (length(ar_poly) > 0 && any(ar_poly != 0)) {
    ar_roots_mod <- Mod(polyroot(c(1, -ar_poly)))
    min_ar_root  <- min(ar_roots_mod)
  } else {
    min_ar_root  <- NA
  }
  
  # MA Roots
  ma_poly <- fit_full$model$theta
  if (length(ma_poly) > 0 && any(ma_poly != 0)) {
    ma_roots_mod <- Mod(polyroot(c(1, ma_poly)))
    min_ma_root  <- min(ma_roots_mod)
  } else {
    min_ma_root  <- NA
  }
  
  is_causal     <- ifelse(!is.na(min_ar_root), min_ar_root > 1, TRUE)
  is_invertible <- ifelse(!is.na(min_ma_root), min_ma_root > 1, TRUE)
  
  cat("\n[2] DIAGNOSTICS:\n")
  cat("Min AR Root Modulus: ", min_ar_root, "\n")
  cat("Min MA Root Modulus: ", min_ma_root, "\n")
  
  # B. Plots
  par(mfrow=c(2,2), mar=c(4,4,3,2))
  plot(resid, main="Residuals", ylab="Error")
  abline(h=0, col="gray"); abline(h=c(-3*sd(resid), 3*sd(resid)), lty=2, col="red")
  scatter.smooth(sqrt(abs(resid)), main="Variance Check", ylab="sqrt(|Resid|)", lpars=list(col="red"))
  qqnorm(resid); qqline(resid, col="red")
  hist(resid, breaks=20, freq=FALSE, main="Histogram"); curve(dnorm(x, mean=mean(resid), sd=sd(resid)), col="red", add=TRUE)
  
  # C. Tests
  sw_test <- shapiro.test(resid)
  lb_test <- Box.test(resid, lag = 20, type = "Ljung-Box")
  
  # -----------------------------------------------------------------------
  # 3. Forecasting & Cross-Validation
  # -----------------------------------------------------------------------
  train_end <- c(floor(time(series)[n-12]), cycle(series)[n-12])
  train_set <- window(series, end = train_end)
  test_set  <- window(series, start = time(series)[n-11]) 
  
  fit_train <- arima(train_set, order = order, seasonal = seasonal)
  fc <- predict(fit_train, n.ahead = 12)
  
  fc_orig <- exp(fc$pred); test_orig <- exp(test_set)
  lower_bound <- exp(fc$pred - 1.96 * fc$se); upper_bound <- exp(fc$pred + 1.96 * fc$se)
  
  # -----------------------------------------------------------------------
  # 4. Summary Table Generation
  # -----------------------------------------------------------------------
  summary_row <- data.frame(
    Model_Name      = model_name,
    AIC             = fit_full$aic,
    BIC             = BIC(fit_full),
    Sigma2          = fit_full$sigma2,
    
    Resid_Norm_PVal = sw_test$p.value,
    Normal       = ifelse(sw_test$p.value > 0.05, "Yes", "No"),
    
    Resid_Ind_PVal  = lb_test$p.value,
    Independent  = ifelse(lb_test$p.value > 0.05, "Yes", "No"),
    
    # Positioned as requested: Root value to the left of the Boolean status
    Min_AR_Root     = min_ar_root,
    Causal          = ifelse(is_causal, "Yes", "No"),
    
    Min_MA_Root     = min_ma_root,
    Invertible      = ifelse(is_invertible, "Yes", "No"),
    
    RMSE            = sqrt(mean((test_orig - fc_orig)^2)),
    MAE             = mean(abs(test_orig - fc_orig)),
    MAPE_Pct        = mean(abs((test_orig - fc_orig) / test_orig)) * 100,
    Mean_CI_Width   = mean(upper_bound - lower_bound)
  )
  
  cat("\n[3] PERFORMANCE SUMMARY:\n")
  print(t(summary_row))
  
  return(summary_row)
}




# -- COMPARATIVE ANALYSIS FUNCTION --


compare_sarima_models <- function(series, models_list) {
  
  cat("\n--- Starting Automated Model Comparison ---\n")
  
  # This will store the summary row for each successful model
  comparison_results <- list()
  
  # Loop through each model specification in the list
  for (i in seq_along(models_list)) {
    spec <- models_list[[i]]
    
    cat(paste0("\nProcessing Model ", i, " of ", length(models_list), "...\n"))
    
    # Use tryCatch to prevent the loop from stopping if a model fails to converge
    result <- tryCatch({
      analyze_sarima(series = series, 
                     order = spec$order, 
                     seasonal = spec$seasonal)
    }, error = function(e) {
      cat(paste0("⚠️ Error in Model ", i, ": ", e$message, "\n"))
      return(NULL) # Return NULL so we can skip it
    })
    
    # If the model was successful, add it to our list
    if (!is.null(result)) {
      comparison_results[[length(comparison_results) + 1]] <- result
    }
  }
  
  # Combine all the rows into one master table
  final_table <- do.call(rbind, comparison_results)
  
  # Sort by BIC by default (it penalizes complexity more than AIC)
  final_table <- final_table[order(final_table$BIC), ]
  
  cat("\n--- Comparison Complete ---\n")
  return(final_table)
}


# -- IMPLEMENTATION --


# Load Data
data <- window(ts(read.table("Unemploy.dat",dec=",")/1e6,start=1996,freq=12),start=c(2000),end=c(2025,12))

# Plot Data
plot(data,main="Unemploy Spain")
abline(v=1996:2026,col=4,lty=3)

# Insert Other Stuff


# Specify Models (for now I just made up some)
my_models <- list(
  list(order = c(1,1,1), 
       seasonal = list(order = c(2,0,2), period = 12)),
  
  list(order = c(1,1,2), 
       seasonal = list(order = c(2,0,2), period = 12)),
  
  list(order = c(0,1,9), 
       seasonal = list(order = c(2,0,2), period = 12)),
  
  list(order = c(7,1,0), 
       seasonal = list(order = c(2,0,2), period = 12))
)

# Run Analysis
comparison_table <- compare_sarima_models(data, my_models)
print(comparison_table)
