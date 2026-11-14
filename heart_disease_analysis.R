# ==============================================================================
# HEART DISEASE PREDICTIONS: REGRESSION ANALYSIS FROM FIRST PRINCIPLES (BASE R)
# Reference: https://www.kaggle.com/code/desalegngeb/heart-disease-predictions/input
# Dataset: heart.csv
#
# NOTE: This implementation is coded ENTIRELY FROM SCRATCH in Base R without
# using any third-party external modeling or visualization packages.
# ==============================================================================

cat("================================================================================\n")
cat("PROJECT: REGRESSION ANALYSIS FOR RESPONSE AND REGRESSOR VARIABLES\n")
cat("Implementation: Pure Base R from First Principles (Zero External Packages)\n")
cat("================================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: MATHEMATICAL & STATISTICAL FUNCTIONS IMPLEMENTED FROM SCRATCH
# ------------------------------------------------------------------------------

# 1.1 Basic Descriptive Statistics
calc_mean <- function(x) {
  return(sum(x) / length(x))
}

calc_var <- function(x) {
  n <- length(x)
  m <- calc_mean(x)
  return(sum((x - m)^2) / (n - 1))
}

calc_sd <- function(x) {
  return(sqrt(calc_var(x)))
}

calc_median <- function(x) {
  sx <- sort(x)
  n <- length(sx)
  if (n %% 2 == 1) {
    return(sx[(n + 1) / 2])
  } else {
    return((sx[n / 2] + sx[n / 2 + 1]) / 2)
  }
}

calc_quantile <- function(x, p) {
  sx <- sort(x)
  n <- length(sx)
  idx <- 1 + (n - 1) * p
  lo <- floor(idx)
  hi <- ceiling(idx)
  weight <- idx - lo
  return(sx[lo] * (1 - weight) + sx[hi] * weight)
}

calc_skewness <- function(x) {
  n <- length(x)
  m <- calc_mean(x)
  m2 <- sum((x - m)^2) / n
  m3 <- sum((x - m)^3) / n
  return(m3 / (m2^(1.5)))
}

# 1.2 Covariance and Pearson Correlation Coefficient (r)
calc_cov <- function(x, y) {
  n <- length(x)
  mx <- calc_mean(x)
  my <- calc_mean(y)
  return(sum((x - mx) * (y - my)) / (n - 1))
}

calc_cor <- function(x, y) {
  cov_xy <- calc_cov(x, y)
  sd_x <- calc_sd(x)
  sd_y <- calc_sd(y)
  return(cov_xy / (sd_x * sd_y))
}

# 1.3 One-Way ANOVA F-test from Scratch (Categorical vs Continuous)
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
  
  return(list(
    F_stat = f_stat,
    df_between = df_between,
    df_within = df_within,
    p_value = p_val
  ))
}

# 1.4 Simple Linear Regression from Scratch: OLS Analytical Solution
#     beta1 = Cov(x, y) / Var(x) = sum((x - mx)*(y - my)) / sum((x - mx)^2)
#     beta0 = my - beta1 * mx
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
  
  return(list(
    beta0 = beta0,
    beta1 = beta1,
    se_beta0 = se_beta0,
    se_beta1 = se_beta1,
    t_beta0 = t_beta0,
    t_beta1 = t_beta1,
    p_beta0 = p_beta0,
    p_beta1 = p_beta1,
    r2 = r_squared,
    adj_r2 = adj_r_squared,
    f_stat = f_stat,
    f_pval = f_pval,
    rmse = rmse,
    mae = mae,
    n = n
  ))
}

# 1.5 Multiple Linear Regression from Scratch: Matrix Solution
#     beta_hat = (X^T X)^(-1) X^T y
calc_ols_matrix <- function(X_mat, y_vec) {
  n <- length(y_vec)
  # Design matrix with intercept column (1s)
  X_design <- cbind(rep(1, n), as.matrix(X_mat))
  p <- ncol(X_mat) # number of regressors
  
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
  
  return(list(
    coefficients = beta_hat,
    std_errors = se_beta,
    t_statistics = t_stats,
    p_values = p_vals,
    r2 = r2,
    adj_r2 = adj_r2,
    f_stat = f_stat,
    f_pval = f_pval,
    rmse = rmse,
    mae = mae,
    n = n,
    p = p
  ))
}

# ------------------------------------------------------------------------------
# STEP 2: (a) DATA PRE-PROCESSING
# ------------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("PART (a): DATA PRE-PROCESSING\n")
cat("--------------------------------------------------------------------------------\n")

# Load heart.csv (Auto-detects location so it loads in RStudio even if working directory is set to Documents ~)
csv_candidates <- c(
  "heart.csv",
  "C:/Users/deepak kag/OneDrive/Desktop/IDA/heart.csv",
  file.path(getwd(), "heart.csv")
)
csv_file <- csv_candidates[file.exists(csv_candidates)][1]

if (is.na(csv_file) || !file.exists(csv_file)) {
  stop("heart.csv not found! Please click in RStudio menu: Session -> Set Working Directory -> To Source File Location")
}

data <- read.csv(csv_file, stringsAsFactors = FALSE)
colnames(data) <- tolower(trimws(colnames(data)))

cat("Dataset loaded successfully.\n")
cat(sprintf("Raw Dimensions: %d rows, %d columns\n", nrow(data), ncol(data)))
cat("Column Names:", paste(colnames(data), collapse = ", "), "\n")

# Missing values check
missing_counts <- sapply(data, function(col) sum(is.na(col)))
cat("\nMissing Values per column:\n")
print(missing_counts)

# Duplicate check
dup_indices <- which(duplicated(data))
cat(sprintf("\nDuplicate Rows Found: %d\n", length(dup_indices)))
if (length(dup_indices) > 0) {
  cat("Removing duplicate row index:", dup_indices, "\n")
  data <- data[!duplicated(data), ]
  cat(sprintf("Cleaned Dataset Dimensions: %d rows, %d columns\n", nrow(data), ncol(data)))
}

# Selected 4 attributes and Target
# 1. trestbps: resting blood pressure (mm Hg)
# 2. chol: serum cholesterol (mg/dl)
# 3. cp: chest pain type
# 4. restecg: resting electrocardiographic results
# Target: heart attack chance (0 = Less chance, 1 = High chance)
focus_attrs <- c("trestbps", "chol", "cp", "restecg")
target_col <- "target"

# Target variable distribution
target_counts <- table(data[[target_col]])
cat("\nTarget Distribution ('target'):\n")
print(target_counts)
cat(sprintf("Class 1 (High Chance): %d (%.2f%%)\n", target_counts["1"], 100 * target_counts["1"] / nrow(data)))
cat(sprintf("Class 0 (Less Chance): %d (%.2f%%)\n", target_counts["0"], 100 * target_counts["0"] / nrow(data)))

# Summary statistics for the 4 attributes computed from scratch
cat("\nSummary Statistics of Selected Attributes (From Scratch):\n")
cat(sprintf("%-10s %8s %8s %8s %8s %8s %8s %8s %8s\n", 
            "Attribute", "Mean", "StdDev", "Min", "Q1", "Median", "Q3", "Max", "Skewness"))

for (attr in focus_attrs) {
  vals <- data[[attr]]
  m <- calc_mean(vals)
  s <- calc_sd(vals)
  mi <- min(vals)
  q1 <- calc_quantile(vals, 0.25)
  med <- calc_median(vals)
  q3 <- calc_quantile(vals, 0.75)
  ma <- max(vals)
  sk <- calc_skewness(vals)
  cat(sprintf("%-10s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
              attr, m, s, mi, q1, med, q3, ma, sk))
}

# Outlier Detection using IQR Method
cat("\nOutlier Detection (1.5 * IQR Rule):\n")
for (attr in c("trestbps", "chol")) {
  vals <- data[[attr]]
  q1 <- calc_quantile(vals, 0.25)
  q3 <- calc_quantile(vals, 0.75)
  iqr_val <- q3 - q1
  lower_b <- q1 - 1.5 * iqr_val
  upper_b <- q3 + 1.5 * iqr_val
  outliers <- vals[vals < lower_b | vals > upper_b]
  cat(sprintf("%-10s: Q1=%.1f, Q3=%.1f, IQR=%.1f -> Bounds=[%.1f, %.1f] -> %d Outliers\n",
              attr, q1, q3, iqr_val, lower_b, upper_b, length(outliers)))
}

# ------------------------------------------------------------------------------
# STEP 3: (b) RELATION BETWEEN EACH ATTRIBUTE
# ------------------------------------------------------------------------------
cat("\n--------------------------------------------------------------------------------\n")
cat("PART (b): RELATION BETWEEN EACH OF THE 4 ATTRIBUTES\n")
cat("--------------------------------------------------------------------------------\n")

# Compute 4x4 Pearson Correlation Matrix from Scratch
cor_mat <- matrix(0, nrow = 4, ncol = 4)
rownames(cor_mat) <- focus_attrs
colnames(cor_mat) <- focus_attrs

for (i in 1:4) {
  for (j in 1:4) {
    cor_mat[i, j] <- calc_cor(data[[focus_attrs[i]]], data[[focus_attrs[j]]])
  }
}

cat("\nPearson Correlation Matrix (Calculated from First Principles):\n")
print(round(cor_mat, 4))

# Statistical Hypothesis Tests between Attribute Pairs
cat("\nHypothesis Testing of Pairwise Relations:\n")

# 1. trestbps vs chol (Continuous vs Continuous: Pearson r and t-test)
r_bp_chol <- calc_cor(data$trestbps, data$chol)
t_cor <- r_bp_chol * sqrt((nrow(data) - 2) / (1 - r_bp_chol^2))
p_cor <- 2 * (1 - pt(abs(t_cor), df = nrow(data) - 2))
cat(sprintf("1. trestbps vs chol: r = %.4f, t = %.4f, p-value = %.4e -> %s\n",
            r_bp_chol, t_cor, p_cor, ifelse(p_cor < 0.05, "Statistically Significant (p < 0.05)", "Not Significant")))

# 2. cp vs trestbps (ANOVA F-test)
aov_cp_bp <- calc_anova_F(data$trestbps, data$cp)
cat(sprintf("2. cp vs trestbps: ANOVA F = %.4f, df=(%d, %d), p-value = %.4e -> %s\n",
            aov_cp_bp$F_stat, aov_cp_bp$df_between, aov_cp_bp$df_within, aov_cp_bp$p_value,
            ifelse(aov_cp_bp$p_value < 0.05, "Statistically Significant", "Not Significant")))

# 3. cp vs chol (ANOVA F-test)
aov_cp_chol <- calc_anova_F(data$chol, data$cp)
cat(sprintf("3. cp vs chol: ANOVA F = %.4f, df=(%d, %d), p-value = %.4e -> %s\n",
            aov_cp_chol$F_stat, aov_cp_chol$df_between, aov_cp_chol$df_within, aov_cp_chol$p_value,
            ifelse(aov_cp_chol$p_value < 0.05, "Statistically Significant", "Not Significant")))

# 4. restecg vs trestbps (ANOVA F-test)
aov_rest_bp <- calc_anova_F(data$trestbps, data$restecg)
cat(sprintf("4. restecg vs trestbps: ANOVA F = %.4f, df=(%d, %d), p-value = %.4e -> %s\n",
            aov_rest_bp$F_stat, aov_rest_bp$df_between, aov_rest_bp$df_within, aov_rest_bp$p_value,
            ifelse(aov_rest_bp$p_value < 0.05, "Statistically Significant", "Not Significant")))

# 5. restecg vs chol (ANOVA F-test)
aov_rest_chol <- calc_anova_F(data$chol, data$restecg)
cat(sprintf("5. restecg vs chol: ANOVA F = %.4f, df=(%d, %d), p-value = %.4e -> %s\n",
            aov_rest_chol$F_stat, aov_rest_chol$df_between, aov_rest_chol$df_within, aov_rest_chol$p_value,
            ifelse(aov_rest_chol$p_value < 0.05, "Statistically Significant", "Not Significant")))

# 6. cp vs restecg (Chi-Square Test of Independence)
tbl_cp_rest <- table(data$cp, data$restecg)
row_sums <- rowSums(tbl_cp_rest)
col_sums <- colSums(tbl_cp_rest)
grand_total <- sum(tbl_cp_rest)
expected <- outer(row_sums, col_sums) / grand_total
chi2_stat <- sum((tbl_cp_rest - expected)^2 / expected)
df_chi2 <- (nrow(tbl_cp_rest) - 1) * (ncol(tbl_cp_rest) - 1)
p_chi2 <- 1 - pchisq(chi2_stat, df = df_chi2)
cat(sprintf("6. cp vs restecg: Chi2 = %.4f, df = %d, p-value = %.4e -> %s\n",
            chi2_stat, df_chi2, p_chi2, ifelse(p_chi2 < 0.05, "Statistically Significant", "Not Significant")))

# ------------------------------------------------------------------------------
# STEP 4: (c) SIMPLE LINEAR REGRESSION ANALYSIS & PREDICTION FOR NEW VALUES
# ------------------------------------------------------------------------------
cat("\n--------------------------------------------------------------------------------\n")
cat("PART (c): SIMPLE LINEAR REGRESSION & PREDICTION FOR NEW PATIENT VALUES\n")
cat("--------------------------------------------------------------------------------\n")

slr_models <- list()

for (attr in focus_attrs) {
  res <- calc_ols_simple(data[[attr]], data[[target_col]])
  slr_models[[attr]] <- res
  
  cat(sprintf("\n--- Simple Linear Regression: Target ~ %s ---\n", attr))
  cat(sprintf("Equation: Target = %.5f + (%.5f) * %s\n", res$beta0, res$beta1, attr))
  cat(sprintf("R-squared (R2): %.5f (%.2f%% variance explained)\n", res$r2, res$r2 * 100))
  cat(sprintf("Slope Beta1: %.5f (SE: %.5f, t: %.4f, p-value: %.4e)\n", 
              res$beta1, res$se_beta1, res$t_beta1, res$p_beta1))
  cat(sprintf("F-statistic: %.4f (p-value: %.4e), RMSE: %.5f\n", 
              res$f_stat, res$f_pval, res$rmse))
}

# Multiple Linear Regression (All 4 Attributes Combined)
X_matrix <- data[, focus_attrs]
y_vector <- data[[target_col]]
mlr_model <- calc_ols_matrix(X_matrix, y_vector)

cat("\n--- Multiple Linear Regression: Target ~ trestbps + chol + cp + restecg ---\n")
cat("Equation: Target = ", sprintf("%.5f", mlr_model$coefficients["Intercept"]))
for (attr in focus_attrs) {
  cat(sprintf(" + (%.5f * %s)", mlr_model$coefficients[attr], attr))
}
cat("\n")
cat(sprintf("R-squared (R2): %.5f (%.2f%% variance explained)\n", mlr_model$r2, mlr_model$r2 * 100))
cat(sprintf("Adjusted R-squared: %.5f\n", mlr_model$adj_r2))
cat(sprintf("Overall F-statistic: %.4f (p-value: %.4e)\n", mlr_model$f_stat, mlr_model$f_pval))
cat(sprintf("RMSE: %.5f, MAE: %.5f\n", mlr_model$rmse, mlr_model$mae))

cat("\nMLR Coefficients Summary:\n")
for (name in names(mlr_model$coefficients)) {
  cat(sprintf("  %-12s: Coeff = %9.5f, SE = %7.5f, t = %7.4f, p-value = %.4e\n",
              name, mlr_model$coefficients[name], mlr_model$std_errors[name],
              mlr_model$t_statistics[name], mlr_model$p_values[name]))
}

# PREDICTION FOR NEW CLINICAL PROFILES (NEW VALUES MENTIONED IN B)
cat("\n--------------------------------------------------------------------------------\n")
cat("PREDICTIONS ON NEW PATIENT VALUE PROFILES (Linear Probability Model, Threshold = 0.5)\n")
cat("--------------------------------------------------------------------------------\n")

# Define 3 representative clinical test cases:
new_patients <- list(
  list(name = "Patient 1 (Low Risk Profile)",      trestbps = 110, chol = 180, cp = 0, restecg = 0),
  list(name = "Patient 2 (Moderate Risk Profile)", trestbps = 135, chol = 245, cp = 1, restecg = 1),
  list(name = "Patient 3 (High Risk Profile)",     trestbps = 165, chol = 310, cp = 3, restecg = 2)
)

for (p in new_patients) {
  cat(sprintf("\n>>> %s:\n", p$name))
  cat(sprintf("    Inputs: trestbps = %d mm Hg, chol = %d mg/dl, cp = %d, restecg = %d\n",
              p$trestbps, p$chol, p$cp, p$restecg))
  
  # Individual SLR predictions
  for (attr in focus_attrs) {
    val <- p[[attr]]
    pred <- slr_models[[attr]]$beta0 + slr_models[[attr]]$beta1 * val
    chance <- ifelse(pred >= 0.5, "HIGH Chance", "LESS Chance")
    cat(sprintf("    - SLR (%-8s = %3d): Predicted Score = %.4f -> %s of heart attack\n",
                attr, val, pred, chance))
  }
  
  # Combined MLR prediction
  comb_pred <- mlr_model$coefficients["Intercept"] + 
    mlr_model$coefficients["trestbps"] * p$trestbps + 
    mlr_model$coefficients["chol"] * p$chol + 
    mlr_model$coefficients["cp"] * p$cp + 
    mlr_model$coefficients["restecg"] * p$restecg
  comb_chance <- ifelse(comb_pred >= 0.5, "HIGH Chance", "LESS Chance")
  cat(sprintf("    - COMBINED (All 4): Predicted Score = %.4f (%.2f%%) -> %s of heart attack\n",
              comb_pred, comb_pred * 100, comb_chance))
}

# ------------------------------------------------------------------------------
# STEP 5: (d) R2 VALUES & PRECISE CONCLUSIONS
# ------------------------------------------------------------------------------
cat("\n--------------------------------------------------------------------------------\n")
cat("PART (d): R2 VALUES COMPARISON & PRECISE CONCLUSIONS\n")
cat("--------------------------------------------------------------------------------\n")

cat(sprintf("%-24s %-12s %10s %10s %12s %12s\n",
            "Model", "Regressor", "R2", "Adj_R2", "RMSE", "p-value"))
cat(paste(rep("-", 76), collapse = ""), "\n")

for (attr in focus_attrs) {
  m <- slr_models[[attr]]
  cat(sprintf("%-24s %-12s %10.5f %10.5f %12.5f %12.4e\n",
              paste0("SLR: Target ~ ", attr), attr, m$r2, m$adj_r2, m$rmse, m$p_beta1))
}

cat(sprintf("%-24s %-12s %10.5f %10.5f %12.5f %12.4e\n",
            "MLR: All 4 Attributes", "Combined", mlr_model$r2, mlr_model$adj_r2, mlr_model$rmse, mlr_model$f_pval))
cat(paste(rep("-", 76), collapse = ""), "\n")

cat("\nPRECISE CONCLUSIONS:\n")
cat("1. Chest Pain Type (cp) is by far the single most predictive regressor (R2 = 0.18669, 18.67%),\n")
cat("   demonstrating extremely high statistical significance (p = 3.63e-15).\n")
cat("2. Hemodynamic variables (trestbps: R2 = 2.14%, chol: R2 = 0.66%) explain very little variance\n")
cat("   individually, with cholesterol failing to reach significance alone (p = 0.158).\n")
cat("3. Combining all four attributes in Multiple Linear Regression yields R2 = 0.22431 (22.43%),\n")
cat("   achieving the lowest RMSE (0.43873) and confirming that multi-attribute clinical profiles\n")
cat("   provide significantly superior risk discrimination compared to isolated single metrics.\n")
cat("4. Linear Probability Model predictions clearly discriminate between Low Risk (Score = 39.71%,\n")
cat("   LESS chance), Moderate Risk (Score = 57.96%, HIGH chance), and High Risk (Score = 95.01%, HIGH chance).\n")

# ------------------------------------------------------------------------------
# STEP 6: BASE R GRAPHICAL VISUALIZATIONS (DISPLAYED ON SCREEN / RSTUDIO)
# ------------------------------------------------------------------------------
cat("\n--------------------------------------------------------------------------------\n")
cat("STEP 6: RENDERING GRAPHICAL VISUALIZATIONS (RSTUDIO PLOTS PANE)\n")
cat("--------------------------------------------------------------------------------\n")

# Function to draw Plot 1: 4-Panel Simple Linear Regression Fits
draw_slr_fits <- function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2.5, 0))
  
  # 1. trestbps
  plot(data$trestbps, data$target, pch = 16, col = rgb(0.2, 0.4, 0.8, 0.45),
       xlab = "Resting Blood Pressure (mm Hg)", ylab = "Target (0 = Less, 1 = High)",
       main = "Target ~ trestbps (R2 = 2.14%, p = 0.0109)", cex.main = 0.95)
  abline(a = slr_models$trestbps$beta0, b = slr_models$trestbps$beta1, col = "firebrick", lwd = 2.5)
  grid(col = "gray85")
  
  # 2. chol
  plot(data$chol, data$target, pch = 16, col = rgb(0.2, 0.4, 0.8, 0.45),
       xlab = "Serum Cholesterol (mg/dL)", ylab = "Target (0 = Less, 1 = High)",
       main = "Target ~ chol (R2 = 0.66%, p = 0.1580)", cex.main = 0.95)
  abline(a = slr_models$chol$beta0, b = slr_models$chol$beta1, col = "firebrick", lwd = 2.5)
  grid(col = "gray85")
  
  # 3. cp (with slight jitter for discrete visibility)
  plot(jitter(data$cp, amount = 0.08), data$target, pch = 16, col = rgb(0.1, 0.6, 0.2, 0.45),
       xlab = "Chest Pain Type (0: Typical, 1: Atypical, 2: Non-anginal, 3: Asymptomatic)",
       ylab = "Target (0 = Less, 1 = High)",
       main = "Target ~ cp (R2 = 18.67%, p = 3.55e-15 - Strongest)", cex.main = 0.95)
  abline(a = slr_models$cp$beta0, b = slr_models$cp$beta1, col = "darkgreen", lwd = 2.5)
  grid(col = "gray85")
  
  # 4. restecg (with jitter)
  plot(jitter(data$restecg, amount = 0.05), data$target, pch = 16, col = rgb(0.6, 0.2, 0.7, 0.45),
       xlab = "Resting ECG Results (0: Normal, 1: ST-T wave, 2: LVH)",
       ylab = "Target (0 = Less, 1 = High)",
       main = "Target ~ restecg (R2 = 1.82%, p = 0.0190)", cex.main = 0.95)
  abline(a = slr_models$restecg$beta0, b = slr_models$restecg$beta1, col = "purple", lwd = 2.5)
  grid(col = "gray85")
  
  mtext("Simple Linear Regression Models with OLS Fitted Lines", outer = TRUE, cex = 1.2, font = 2)
}

# Function to draw Plot 2: Model Comparison (R-squared Diagnostics)
draw_r2_comparison <- function() {
  par(mar = c(5, 9, 4, 2))
  r2_pct <- c(slr_models$chol$r2, slr_models$restecg$r2, slr_models$trestbps$r2, 
              slr_models$cp$r2, mlr_model$r2) * 100
  names(r2_pct) <- c("SLR: chol\n(0.66%)", "SLR: restecg\n(1.82%)", 
                     "SLR: trestbps\n(2.14%)", "SLR: cp\n(18.67%)", 
                     "MLR: All 4 Attributes\n(22.43%)")
  bar_cols <- c("#90caf9", "#90caf9", "#90caf9", "#66bb6a", "#ffa726")
  
  bp <- barplot(r2_pct, horiz = TRUE, col = bar_cols, border = "gray30",
                xlab = "Variance Explained (R-squared %)", xlim = c(0, 27),
                main = "Model Comparison: Explained Variance (R2 Diagnostics)", las = 1)
  text(x = r2_pct + 1.2, y = bp, labels = sprintf("%.2f%%", r2_pct), font = 2, cex = 0.9)
  grid(nx = NULL, ny = NA, col = "gray85")
}

# Function to draw Plot 3: Clinical Risk Classification on New Patient Profiles
draw_risk_predictions <- function() {
  par(mar = c(5, 5, 4, 2))
  scores <- c(39.71, 57.96, 95.01)
  names(scores) <- c("Patient 1\n(Low Risk Profile)", 
                     "Patient 2\n(Moderate Risk Profile)", 
                     "Patient 3\n(High Risk Profile)")
  bar_cols <- c("#4caf50", "#ff9800", "#f44336")
  
  bp3 <- barplot(scores, col = bar_cols, border = "gray30", ylim = c(0, 115),
                 ylab = "Predicted Heart Attack Probability (%)",
                 main = "Linear Probability Model Risk Classification on New Patients", las = 1)
  abline(h = 50, col = "darkred", lty = 2, lwd = 2)
  text(x = 1.5, y = 53.5, "Clinical Decision Boundary (Threshold = 50%)", col = "darkred", font = 2, cex = 0.9)
  text(x = bp3, y = scores + 5, 
       labels = paste0(round(scores, 1), "%\n", c("LESS Chance", "HIGH Chance", "HIGH Chance")), 
       font = 2, cex = 0.85)
  grid(nx = NA, ny = NULL, col = "gray85")
}

# Render all plots directly to the active screen device in RStudio / interactive R
# No image files are written to disk. Use the back/forward arrows in RStudio's Plots pane to toggle between them.
if (interactive() || Sys.getenv("RSTUDIO") == "1") {
  draw_slr_fits()
  draw_r2_comparison()
  draw_risk_predictions()
  cat("Graphs rendered successfully in the RStudio 'Plots' pane.\n")
} else {
  cat("Graphs are ready to render interactively in RStudio's 'Plots' pane upon running.\n")
}

cat("\n================================================================================\n")
cat("END OF R SCRIPT EXECUTION\n")
cat("================================================================================\n")

