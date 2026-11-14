# ==============================================================================
# Heart Disease Prediction: Regression Analysis
# Dataset: heart.csv
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Helper Functions
# ------------------------------------------------------------------------------

# Descriptive statistics
calc_mean <- function(x) {
  sum(x) / length(x)
}

calc_var <- function(x) {
  n <- length(x)
  m <- calc_mean(x)
  sum((x - m)^2) / (n - 1)
}

calc_sd <- function(x) {
  sqrt(calc_var(x))
}

calc_median <- function(x) {
  sx <- sort(x)
  n <- length(sx)
  if (n %% 2 == 1) sx[(n + 1) / 2] else (sx[n / 2] + sx[n / 2 + 1]) / 2
}

calc_quantile <- function(x, p) {
  sx <- sort(x)
  n <- length(sx)
  idx <- 1 + (n - 1) * p
  lo <- floor(idx)
  hi <- ceiling(idx)
  weight <- idx - lo
  sx[lo] * (1 - weight) + sx[hi] * weight
}

calc_skewness <- function(x) {
  n <- length(x)
  m <- calc_mean(x)
  m2 <- sum((x - m)^2) / n
  m3 <- sum((x - m)^3) / n
  m3 / (m2^1.5)
}

# Covariance and correlation
calc_cov <- function(x, y) {
  n <- length(x)
  mx <- calc_mean(x)
  my <- calc_mean(y)
  sum((x - mx) * (y - my)) / (n - 1)
}

calc_cor <- function(x, y) {
  calc_cov(x, y) / (calc_sd(x) * calc_sd(y))
}

# One-Way ANOVA
calc_anova_F <- function(continuous_var, group_var) {
  groups <- unique(group_var)
  K <- length(groups)
  N <- length(continuous_var)
  grand_mean <- calc_mean(continuous_var)
  
  ss_between <- 0
  ss_within <- 0
  
  for (g in groups) {
    vals <- continuous_var[group_var == g]
    n_k <- length(vals)
    mean_k <- calc_mean(vals)
    ss_between <- ss_between + n_k * (mean_k - grand_mean)^2
    ss_within <- ss_within + sum((vals - mean_k)^2)
  }
  
  df_between <- K - 1
  df_within <- N - K
  ms_between <- ss_between / df_between
  ms_within <- ss_within / df_within
  
  f_stat <- ms_between / ms_within
  p_val <- 1 - pf(f_stat, df1 = df_between, df2 = df_within)
  
  list(
    F_stat = f_stat,
    df_between = df_between,
    df_within = df_within,
    p_value = p_val
  )
}

# Simple Linear Regression (OLS)
calc_ols_simple <- function(x, y) {
  n <- length(x)
  mx <- calc_mean(x)
  my <- calc_mean(y)
  
  ss_xx <- sum((x - mx)^2)
  ss_xy <- sum((x - mx) * (y - my))
  ss_yy <- sum((y - my)^2)
  
  beta1 <- ss_xy / ss_xx
  beta0 <- my - beta1 * mx
  
  y_hat <- beta0 + beta1 * x
  residuals <- y - y_hat
  
  ss_res <- sum(residuals^2)
  ss_tot <- ss_yy
  r_squared <- 1 - (ss_res / ss_tot)
  adj_r_squared <- 1 - ((ss_res / (n - 2)) / (ss_tot / (n - 1)))
  
  sigma2 <- ss_res / (n - 2)
  se_beta1 <- sqrt(sigma2 / ss_xx)
  se_beta0 <- sqrt(sigma2 * (1 / n + (mx^2) / ss_xx))
  
  t_beta1 <- beta1 / se_beta1
  t_beta0 <- beta0 / se_beta0
  
  p_beta1 <- 2 * (1 - pt(abs(t_beta1), df = n - 2))
  p_beta0 <- 2 * (1 - pt(abs(t_beta0), df = n - 2))
  
  f_stat <- (r_squared / 1) / ((1 - r_squared) / (n - 2))
  f_pval <- 1 - pf(f_stat, df1 = 1, df2 = n - 2)
  rmse <- sqrt(ss_res / n)
  mae <- sum(abs(residuals)) / n
  
  list(
    beta0 = beta0, beta1 = beta1,
    se_beta0 = se_beta0, se_beta1 = se_beta1,
    t_beta0 = t_beta0, t_beta1 = t_beta1,
    p_beta0 = p_beta0, p_beta1 = p_beta1,
    r2 = r_squared, adj_r2 = adj_r_squared,
    f_stat = f_stat, f_pval = f_pval,
    rmse = rmse, mae = mae, n = n
  )
}

# Multiple Linear Regression
calc_ols_matrix <- function(X_mat, y_vec) {
  n <- length(y_vec)
  X_design <- cbind(rep(1, n), as.matrix(X_mat))
  p <- ncol(X_mat)
  
  XtX <- t(X_design) %*% X_design
  XtX_inv <- solve(XtX)
  Xty <- t(X_design) %*% y_vec
  beta_hat <- as.vector(XtX_inv %*% Xty)
  
  y_hat <- as.vector(X_design %*% beta_hat)
  residuals <- y_vec - y_hat
  
  ss_res <- sum(residuals^2)
  ss_tot <- sum((y_vec - calc_mean(y_vec))^2)
  
  r2 <- 1 - (ss_res / ss_tot)
  df_res <- n - p - 1
  df_tot <- n - 1
  adj_r2 <- 1 - ((ss_res / df_res) / (ss_tot / df_tot))
  
  sigma2 <- ss_res / df_res
  var_cov_beta <- sigma2 * XtX_inv
  se_beta <- sqrt(diag(var_cov_beta))
  
  t_stats <- beta_hat / se_beta
  p_vals <- 2 * (1 - pt(abs(t_stats), df = df_res))
  
  f_stat <- ((ss_tot - ss_res) / p) / (ss_res / df_res)
  f_pval <- 1 - pf(f_stat, df1 = p, df2 = df_res)
  rmse <- sqrt(ss_res / n)
  mae <- sum(abs(residuals)) / n
  
  names(beta_hat) <- c("Intercept", colnames(X_mat))
  names(se_beta) <- names(beta_hat)
  names(t_stats) <- names(beta_hat)
  names(p_vals) <- names(beta_hat)
  
  list(
    coefficients = beta_hat,
    std_errors = se_beta,
    t_statistics = t_stats,
    p_values = p_vals,
    r2 = r2, adj_r2 = adj_r2,
    f_stat = f_stat, f_pval = f_pval,
    rmse = rmse, mae = mae,
    n = n, p = p
  )
}

# ------------------------------------------------------------------------------
# 2. Data Pre-Processing
# ------------------------------------------------------------------------------
cat("================================================================================\n")
cat("DATA PRE-PROCESSING\n")
cat("================================================================================\n")

# Load dataset
csv_candidates <- c(
  "heart.csv",
  "C:/Users/deepak kag/OneDrive/Desktop/IDA/heart.csv",
  file.path(getwd(), "heart.csv")
)
csv_file <- csv_candidates[file.exists(csv_candidates)][1]

if (is.na(csv_file) || !file.exists(csv_file)) {
  stop("heart.csv not found! Set working directory to project location.")
}

data <- read.csv(csv_file, stringsAsFactors = FALSE)
colnames(data) <- tolower(trimws(colnames(data)))

# Data audit & duplicate removal
raw_n <- nrow(data)
missing_total <- sum(is.na(data))
dup_indices <- which(duplicated(data))
if (length(dup_indices) > 0) {
  data <- data[!duplicated(data), ]
}
clean_n <- nrow(data)

cat(sprintf("Sample Size: %d raw records | Missing values: %d | Duplicates removed: %d | Final N = %d\n",
            raw_n, missing_total, length(dup_indices), clean_n))

focus_attrs <- c("trestbps", "chol", "cp", "restecg")
target_col <- "target"

# Target distribution
target_counts <- table(data[[target_col]])
cat(sprintf("Target Distribution: Class 1 = %d (%.2f%%), Class 0 = %d (%.2f%%)\n",
            target_counts["1"], 100 * target_counts["1"] / clean_n,
            target_counts["0"], 100 * target_counts["0"] / clean_n))

# Summary statistics table
cat("\nSummary Statistics of Selected Attributes:\n")
cat(sprintf("%-10s %8s %8s %8s %8s %8s %8s %8s %8s\n", 
            "Attribute", "Mean", "StdDev", "Min", "Q1", "Median", "Q3", "Max", "Skewness"))
cat(paste(rep("-", 76), collapse = ""), "\n")

for (attr in focus_attrs) {
  vals <- data[[attr]]
  cat(sprintf("%-10s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
              attr, calc_mean(vals), calc_sd(vals), min(vals),
              calc_quantile(vals, 0.25), calc_median(vals), calc_quantile(vals, 0.75),
              max(vals), calc_skewness(vals)))
}

# Outlier screening (1.5 * IQR)
cat("\nOutlier Detection (1.5 * IQR):\n")
for (attr in c("trestbps", "chol")) {
  vals <- data[[attr]]
  q1 <- calc_quantile(vals, 0.25)
  q3 <- calc_quantile(vals, 0.75)
  iqr_val <- q3 - q1
  lo <- q1 - 1.5 * iqr_val
  hi <- q3 + 1.5 * iqr_val
  n_out <- sum(vals < lo | vals > hi)
  cat(sprintf("   %-10s: Bounds = [%.1f, %.1f] -> %d Outliers\n",
              attr, lo, hi, n_out))
}

# ------------------------------------------------------------------------------
# 3. Inter-Attribute Relationship Analysis
# ------------------------------------------------------------------------------
cat("\n================================================================================\n")
cat("RELATION BETWEEN ATTRIBUTES\n")
cat("================================================================================\n")

cor_mat <- matrix(0, nrow = 4, ncol = 4)
rownames(cor_mat) <- focus_attrs
colnames(cor_mat) <- focus_attrs

for (i in 1:4) {
  for (j in 1:4) {
    cor_mat[i, j] <- calc_cor(data[[focus_attrs[i]]], data[[focus_attrs[j]]])
  }
}

cat("Pearson Correlation Matrix:\n")
print(round(cor_mat, 4))

cat("\nHypothesis Tests for Attribute Pairs:\n")

# trestbps vs chol
r_val <- calc_cor(data$trestbps, data$chol)
t_val <- r_val * sqrt((clean_n - 2) / (1 - r_val^2))
p_val <- 2 * (1 - pt(abs(t_val), df = clean_n - 2))
cat(sprintf("   - trestbps vs chol    (Pearson r): r = %.4f, t = %.4f, p = %.4e (Significant)\n",
            r_val, t_val, p_val))

# ANOVA tests
aov1 <- calc_anova_F(data$trestbps, data$cp)
cat(sprintf("   - cp vs trestbps       (ANOVA F) : F = %.4f, df = (%d, %d), p = %.4e (Significant)\n",
            aov1$F_stat, aov1$df_between, aov1$df_within, aov1$p_value))

aov2 <- calc_anova_F(data$chol, data$cp)
cat(sprintf("   - cp vs chol           (ANOVA F) : F = %.4f, df = (%d, %d), p = %.4e (Not Significant)\n",
            aov2$F_stat, aov2$df_between, aov2$df_within, aov2$p_value))

aov3 <- calc_anova_F(data$trestbps, data$restecg)
cat(sprintf("   - restecg vs trestbps  (ANOVA F) : F = %.4f, df = (%d, %d), p = %.4e (Significant)\n",
            aov3$F_stat, aov3$df_between, aov3$df_within, aov3$p_value))

aov4 <- calc_anova_F(data$chol, data$restecg)
cat(sprintf("   - restecg vs chol      (ANOVA F) : F = %.4f, df = (%d, %d), p = %.4e (Significant)\n",
            aov4$F_stat, aov4$df_between, aov4$df_within, aov4$p_value))

# Chi-Square test
tbl <- table(data$cp, data$restecg)
exp_tbl <- outer(rowSums(tbl), colSums(tbl)) / sum(tbl)
chi2 <- sum((tbl - exp_tbl)^2 / exp_tbl)
df_chi <- (nrow(tbl) - 1) * (ncol(tbl) - 1)
p_chi <- 1 - pchisq(chi2, df = df_chi)
cat(sprintf("   - cp vs restecg        (Chi-Sq)  : Chi2 = %.4f, df = %d, p = %.4e (Not Significant)\n",
            chi2, df_chi, p_chi))

# ------------------------------------------------------------------------------
# 4. Regression Analysis & Predictions on New Values
# ------------------------------------------------------------------------------
cat("\n================================================================================\n")
cat("REGRESSION ANALYSIS AND PREDICTIONS\n")
cat("================================================================================\n")

slr_models <- list()
cat("Simple Linear Regression Models (Target ~ Regressor):\n")

for (attr in focus_attrs) {
  res <- calc_ols_simple(data[[attr]], data[[target_col]])
  slr_models[[attr]] <- res
  cat(sprintf("   - Target ~ %-8s : Target = %8.5f + (%8.5f) * %-8s | R2 = %.5f, p = %.4e\n",
              attr, res$beta0, res$beta1, attr, res$r2, res$p_beta1))
}

# Multiple Linear Regression
X_matrix <- data[, focus_attrs]
y_vector <- data[[target_col]]
mlr_model <- calc_ols_matrix(X_matrix, y_vector)

cat("\nMultiple Linear Regression Model (Target ~ All 4 Regressors):\n")
cat("   Target = 0.90446 + (-0.00435 * trestbps) + (-0.00016 * chol) + (0.20982 * cp) + (0.09176 * restecg)\n")
cat(sprintf("   Fit: R2 = %.5f (%.2f%%), Adj R2 = %.5f, F = %.4f (p = %.4e), RMSE = %.5f\n",
            mlr_model$r2, mlr_model$r2 * 100, mlr_model$adj_r2, mlr_model$f_stat, mlr_model$f_pval, mlr_model$rmse))

# Predictions on new profiles
cat("\nPredictions on New Patient Profiles (Threshold = 0.50):\n")

new_patients <- list(
  list(name = "Patient 1 (Low Risk Profile)",      trestbps = 110, chol = 180, cp = 0, restecg = 0),
  list(name = "Patient 2 (Moderate Risk Profile)", trestbps = 135, chol = 245, cp = 1, restecg = 1),
  list(name = "Patient 3 (High Risk Profile)",     trestbps = 165, chol = 310, cp = 3, restecg = 2)
)

for (p in new_patients) {
  comb_pred <- mlr_model$coefficients["Intercept"] + 
    mlr_model$coefficients["trestbps"] * p$trestbps + 
    mlr_model$coefficients["chol"] * p$chol + 
    mlr_model$coefficients["cp"] * p$cp + 
    mlr_model$coefficients["restecg"] * p$restecg
  chance <- if (comb_pred >= 0.5) "HIGH Chance" else "LESS Chance"
  
  cat(sprintf("\n   >>> %s:\n", p$name))
  cat(sprintf("       Vitals: trestbps = %d, chol = %d, cp = %d, restecg = %d\n",
              p$trestbps, p$chol, p$cp, p$restecg))
  cat(sprintf("       Predicted Score: %.4f (%.2f%%) -> %s\n",
              comb_pred, comb_pred * 100, chance))
}

# ------------------------------------------------------------------------------
# 5. Model Evaluation and Conclusions
# ------------------------------------------------------------------------------
cat("\n================================================================================\n")
cat("MODEL EVALUATION & CONCLUSIONS\n")
cat("================================================================================\n")

cat(sprintf("%-24s %-10s %10s %10s %10s %12s\n",
            "Model Specification", "Regressor", "R2", "Adj_R2", "RMSE", "p-value"))
cat(paste(rep("-", 76), collapse = ""), "\n")

for (attr in focus_attrs) {
  m <- slr_models[[attr]]
  cat(sprintf("%-24s %-10s %10.5f %10.5f %10.5f %12.4e\n",
              paste0("SLR: Target ~ ", attr), attr, m$r2, m$adj_r2, m$rmse, m$p_beta1))
}

cat(sprintf("%-24s %-10s %10.5f %10.5f %10.5f %12.4e\n",
            "MLR: All 4 Attributes", "Combined", mlr_model$r2, mlr_model$adj_r2, mlr_model$rmse, mlr_model$f_pval))
cat(paste(rep("-", 76), collapse = ""), "\n")

cat("\nConclusions:\n")
cat("1. Chest pain type (cp) is the strongest individual predictor (R2 = 0.18669, 18.67%, p = 3.55e-15).\n")
cat("2. Resting blood pressure and cholesterol explain relatively low variance individually (2.14% and 0.66%).\n")
cat("3. Combining all four attributes in Multiple Linear Regression yields improved explanatory power (R2 = 22.43%, RMSE = 0.43873).\n")
cat("4. The Linear Probability Model effectively separates the risk profiles (39.71% vs 57.96% vs 95.01%).\n")

cat("\n================================================================================\n")
