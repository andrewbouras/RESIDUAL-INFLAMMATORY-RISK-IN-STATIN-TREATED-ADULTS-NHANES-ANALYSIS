# =============================================================================
# RIR Manuscript - Table Generation
# =============================================================================
# Generates publication-quality tables for the RIR manuscript
# Table 1: Characteristics by RIR status
# Table 2: 6-group LDL×Statin prevalence
# Table 3: Regression results
# =============================================================================

library(survey)
library(tidyverse)
library(haven)

# Handle lonely PSUs
options(survey.lonely.psu = "adjust")

# -----------------------------------------------------------------------------
# Setup paths
# -----------------------------------------------------------------------------
project_root <- tryCatch({
  dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
}, error = function(e) {
  getwd()
})

if (is.null(project_root) || project_root == "" || project_root == ".") {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) {
    project_root <- dirname(dirname(normalizePath(script_path)))
  } else {
    project_root <- getwd()
  }
}

data_path <- file.path(project_root, "data", "processed", "rir_analytic_cohort.csv")
output_dir <- file.path(project_root, "manuscript", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("="
, strrep("=", 69), "\n")
cat("RIR MANUSCRIPT - TABLE GENERATION\n")
cat(strrep("=", 70), "\n\n")

# -----------------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------------
cat("[Loading data]\n")
df <- read_csv(data_path, show_col_types = FALSE)
cat("  Loaded", nrow(df), "participants\n")

# Prepare variables
df <- df %>%
  mutate(
    sex_cat = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    race_cat = factor(race_eth, 
                      levels = c(1, 2, 3, 4, 5),
                      labels = c("Mexican American", "Other Hispanic", 
                                "Non-Hispanic White", "Non-Hispanic Black", 
                                "Other/Multi-racial")),
    age_group = cut(age, 
                    breaks = c(20, 40, 50, 60, 70, Inf),
                    labels = c("20-39", "40-49", "50-59", "60-69", "70+"),
                    right = FALSE),
    bmi_cat = cut(bmi,
                  breaks = c(0, 25, 30, Inf),
                  labels = c("Normal/Underweight", "Overweight", "Obese"),
                  right = FALSE),
    glycemic_cat = factor(glycemic_status,
                          levels = c(0, 1, 2),
                          labels = c("Normal", "Prediabetes", "Diabetes")),
    rir_cat = factor(rir, levels = c(0, 1), labels = c("No RIR", "RIR")),
    diabetes_bin = as.numeric(diabetes == 1),
    hypertension_bin = as.numeric(hypertension == 1),
    obese_bin = as.numeric(obese == 1),
    current_smoker_bin = as.numeric(current_smoker == 1),
    prediabetes_bin = as.numeric(prediabetes == 1)
  )

# Adjust weights for multi-cycle
n_cycles <- length(unique(df$cycle))
df$adj_weight <- df$fasting_weight / n_cycles

# Create survey design - full cohort
design_full <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~adj_weight,
  data = df,
  nest = TRUE
)

# Primary cohort: statin users with LDL<70
design_primary <- subset(design_full, statin_user == 1 & ldl_under_70 == 1)

cat("  Primary cohort N:", sum(df$statin_user == 1 & df$ldl_under_70 == 1, na.rm = TRUE), "\n\n")

# =============================================================================
# TABLE 1: Characteristics by RIR Status
# =============================================================================
cat(strrep("=", 70), "\n")
cat("TABLE 1: CHARACTERISTICS BY RIR STATUS\n")
cat(strrep("=", 70), "\n\n")

# Function to calculate weighted mean and SD by group
calc_continuous <- function(design, varname) {
  formula <- as.formula(paste0("~", varname))
  
  # Means by RIR status
  means <- tryCatch({
    svyby(formula, ~rir_cat, design = design, FUN = svymean, na.rm = TRUE)
  }, error = function(e) NULL)
  
  # SDs (using svyvar)
  vars <- tryCatch({
    svyby(formula, ~rir_cat, design = design, FUN = svyvar, na.rm = TRUE)
  }, error = function(e) NULL)
  
  # P-value from svyttest
  pval <- tryCatch({
    test <- svyttest(as.formula(paste0(varname, " ~ rir")), design = design)
    test$p.value
  }, error = function(e) NA)
  
  if (is.null(means) || is.null(vars)) return(NULL)
  
  # Extract values
  no_rir_mean <- means[1, varname]
  no_rir_sd <- sqrt(vars[1, varname])
  rir_mean <- means[2, varname]
  rir_sd <- sqrt(vars[2, varname])
  
  data.frame(
    Variable = varname,
    No_RIR = sprintf("%.1f (%.1f)", no_rir_mean, no_rir_sd),
    RIR = sprintf("%.1f (%.1f)", rir_mean, rir_sd),
    P_value = ifelse(is.na(pval), "—", sprintf("%.3f", pval))
  )
}

# Function to calculate weighted proportions by group
calc_categorical <- function(design, varname, level_label = NULL) {
  formula <- as.formula(paste0("~", varname))
  
  # Proportions by RIR status
  props <- tryCatch({
    svyby(formula, ~rir_cat, design = design, FUN = svymean, na.rm = TRUE)
  }, error = function(e) NULL)
  
  # Chi-square test (skip for indicator variables)
  pval <- tryCatch({
    test <- svyttest(as.formula(paste0(varname, " ~ rir")), design = design)
    test$p.value
  }, error = function(e) NA)
  
  if (is.null(props)) return(NULL)
  
  # For indicator variables like I(sex == 2), the column name is the full expression
  # Find the column that contains TRUE or the varname
  col_idx <- which(grepl("TRUE|^[^s]", names(props)) & !grepl("rir_cat|^se", names(props)))
  if (length(col_idx) == 0) {
    col_idx <- 2  # Default to second column

  } else {
    col_idx <- col_idx[1]
  }
  
  no_rir_pct <- props[1, col_idx] * 100
  rir_pct <- props[2, col_idx] * 100
  
  label <- ifelse(is.null(level_label), varname, level_label)
  
  data.frame(
    Variable = label,
    No_RIR = sprintf("%.1f%%", no_rir_pct),
    RIR = sprintf("%.1f%%", rir_pct),
    P_value = ifelse(is.na(pval), "—", sprintf("%.3f", pval))
  )
}

# Build Table 1
table1_rows <- list()

# Sample sizes
n_no_rir <- sum(df$statin_user == 1 & df$ldl_under_70 == 1 & df$rir == 0, na.rm = TRUE)
n_rir <- sum(df$statin_user == 1 & df$ldl_under_70 == 1 & df$rir == 1, na.rm = TRUE)

table1_rows[[1]] <- data.frame(
  Variable = "N (unweighted)",
  No_RIR = as.character(n_no_rir),
  RIR = as.character(n_rir),
  P_value = ""
)

# Continuous variables
cat("Calculating continuous variables...\n")
continuous_vars <- c("age", "bmi", "ldl_calc", "hdl", "trigly", "hscrp", 
                     "sbp_mean", "dbp_mean", "hba1c")
var_labels <- c("Age, years", "BMI, kg/m²", "LDL-C, mg/dL", "HDL-C, mg/dL",
                "Triglycerides, mg/dL", "hs-CRP, mg/L", "Systolic BP, mmHg",
                "Diastolic BP, mmHg", "HbA1c, %")

for (i in seq_along(continuous_vars)) {
  result <- calc_continuous(design_primary, continuous_vars[i])
  if (!is.null(result)) {
    result$Variable <- var_labels[i]
    table1_rows[[length(table1_rows) + 1]] <- result
  }
}

# Categorical variables
cat("Calculating categorical variables...\n")

# Female
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "I(sex == 2)", "Female, %"
)

# Race/ethnicity - Non-Hispanic White
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "I(race_eth == 3)", "Non-Hispanic White, %"
)

# Race/ethnicity - Non-Hispanic Black
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "I(race_eth == 4)", "Non-Hispanic Black, %"
)

# Diabetes
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "diabetes_bin", "Diabetes, %"
)

# Prediabetes
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "prediabetes_bin", "Prediabetes, %"
)

# Hypertension
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "hypertension_bin", "Hypertension, %"
)

# Obesity
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "obese_bin", "Obesity (BMI ≥30), %"
)

# Current smoker
table1_rows[[length(table1_rows) + 1]] <- calc_categorical(
  design_primary, "current_smoker_bin", "Current smoker, %"
)

# Combine all rows
table1 <- bind_rows(table1_rows)

# Print table
cat("\n")
cat("TABLE 1: Baseline Characteristics of Statin Users with LDL-C <70 mg/dL\n")
cat("        by Residual Inflammatory Risk Status\n")
cat(strrep("-", 70), "\n")
print(table1, row.names = FALSE)

# Save to CSV
write_csv(table1, file.path(output_dir, "Table1_Characteristics.csv"))
cat("\nSaved: Table1_Characteristics.csv\n")

# =============================================================================
# TABLE 2: 6-Group LDL×Statin Prevalence
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE 2: hs-CRP ≥2 mg/L PREVALENCE BY LDL×STATIN GROUP\n")
cat(strrep("=", 70), "\n\n")

# Create group labels
df <- df %>%
  mutate(
    ldl_statin_group_f = factor(ldl_statin_group,
                                levels = c(1, 2, 3, 4, 5, 6),
                                labels = c("LDL ≤70, Statin", "LDL ≤70, No Statin",
                                          "LDL 70-130, Statin", "LDL 70-130, No Statin",
                                          "LDL >130, Statin", "LDL >130, No Statin"))
  )

# Recreate design with new variable
design_full <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~adj_weight,
  data = df,
  nest = TRUE
)

# Calculate prevalence by group
table2_rows <- list()

groups <- c("LDL ≤70, Statin", "LDL ≤70, No Statin",
            "LDL 70-130, Statin", "LDL 70-130, No Statin",
            "LDL >130, Statin", "LDL >130, No Statin")

for (grp in groups) {
  grp_design <- subset(design_full, ldl_statin_group_f == grp)
  n_grp <- nrow(grp_design)
  
  if (n_grp > 0) {
    prev <- svymean(~hscrp_high, design = grp_design, na.rm = TRUE)
    ci <- confint(prev)
    
    table2_rows[[length(table2_rows) + 1]] <- data.frame(
      Group = grp,
      N = n_grp,
      Prevalence = sprintf("%.1f%%", coef(prev) * 100),
      CI_95 = sprintf("(%.1f%% - %.1f%%)", ci[1] * 100, ci[2] * 100),
      SE = sprintf("%.1f%%", sqrt(attr(prev, "var")) * 100)
    )
  }
}

table2 <- bind_rows(table2_rows)

cat("hs-CRP ≥2 mg/L Prevalence Across LDL-C and Statin Status Groups\n")
cat(strrep("-", 70), "\n")
print(table2, row.names = FALSE)

# Save
write_csv(table2, file.path(output_dir, "Table2_LDL_Statin_Groups.csv"))
cat("\nSaved: Table2_LDL_Statin_Groups.csv\n")

# =============================================================================
# TABLE 2B: Stratified by Glycemic Status
# =============================================================================
cat("\n", strrep("-", 50), "\n")
cat("Table 2B: Stratified by Glycemic Status\n")
cat(strrep("-", 50), "\n\n")

table2b_rows <- list()

for (grp in groups) {
  for (glyc in c("Normal", "Prediabetes", "Diabetes")) {
    grp_design <- subset(design_full, ldl_statin_group_f == grp & glycemic_cat == glyc)
    n_grp <- nrow(grp_design)
    
    if (n_grp >= 10) {
      prev <- tryCatch({
        svymean(~hscrp_high, design = grp_design, na.rm = TRUE)
      }, error = function(e) NULL)
      
      if (!is.null(prev)) {
        table2b_rows[[length(table2b_rows) + 1]] <- data.frame(
          Group = grp,
          Glycemic_Status = glyc,
          N = n_grp,
          Prevalence = sprintf("%.1f%%", coef(prev) * 100),
          SE = sprintf("%.1f%%", sqrt(attr(prev, "var")) * 100)
        )
      }
    }
  }
}

table2b <- bind_rows(table2b_rows)
print(table2b, row.names = FALSE)

write_csv(table2b, file.path(output_dir, "Table2B_Glycemic_Stratified.csv"))
cat("\nSaved: Table2B_Glycemic_Stratified.csv\n")

# =============================================================================
# TABLE 2C: Stratified by Sex
# =============================================================================
cat("\n", strrep("-", 50), "\n")
cat("Table 2C: Stratified by Sex\n")
cat(strrep("-", 50), "\n\n")

table2c_rows <- list()

for (grp in groups) {
  for (sx in c("Male", "Female")) {
    grp_design <- subset(design_full, ldl_statin_group_f == grp & sex_cat == sx)
    n_grp <- nrow(grp_design)
    
    if (n_grp >= 10) {
      prev <- tryCatch({
        svymean(~hscrp_high, design = grp_design, na.rm = TRUE)
      }, error = function(e) NULL)
      
      if (!is.null(prev)) {
        table2c_rows[[length(table2c_rows) + 1]] <- data.frame(
          Group = grp,
          Sex = sx,
          N = n_grp,
          Prevalence = sprintf("%.1f%%", coef(prev) * 100),
          SE = sprintf("%.1f%%", sqrt(attr(prev, "var")) * 100)
        )
      }
    }
  }
}

table2c <- bind_rows(table2c_rows)
print(table2c, row.names = FALSE)

write_csv(table2c, file.path(output_dir, "Table2C_Sex_Stratified.csv"))
cat("\nSaved: Table2C_Sex_Stratified.csv\n")

# =============================================================================
# TABLE 3: Regression Results
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE 3: PREDICTORS OF RIR (LOGISTIC REGRESSION)\n")
cat(strrep("=", 70), "\n\n")

design_primary <- subset(design_full, statin_user == 1 & ldl_under_70 == 1)

# Fit model
model <- tryCatch({
  svyglm(rir ~ age + sex_cat + race_cat + bmi + diabetes_bin + 
           hypertension_bin + current_smoker_bin + trigly,
         design = design_primary,
         family = quasibinomial())
}, error = function(e) {
  cat("Error fitting model:", e$message, "\n")
  NULL
})

if (!is.null(model)) {
  # Extract results
  coefs <- summary(model)$coefficients
  or <- exp(coef(model))
  ci <- exp(confint(model))
  
  table3 <- data.frame(
    Predictor = rownames(coefs),
    OR = sprintf("%.2f", or),
    CI_95 = sprintf("(%.2f - %.2f)", ci[, 1], ci[, 2]),
    P_value = sprintf("%.3f", coefs[, "Pr(>|t|)"])
  )
  
  # Clean up predictor names
  table3$Predictor <- gsub("sex_catFemale", "Female (vs Male)", table3$Predictor)
  table3$Predictor <- gsub("race_cat", "", table3$Predictor)
  table3$Predictor <- gsub("diabetes_bin", "Diabetes", table3$Predictor)
  table3$Predictor <- gsub("hypertension_bin", "Hypertension", table3$Predictor)
  table3$Predictor <- gsub("current_smoker_bin", "Current Smoker", table3$Predictor)
  table3$Predictor <- gsub("trigly", "Triglycerides (per mg/dL)", table3$Predictor)
  table3$Predictor <- gsub("bmi", "BMI (per kg/m²)", table3$Predictor)
  table3$Predictor <- gsub("age", "Age (per year)", table3$Predictor)
  
  # Remove intercept
  table3 <- table3[table3$Predictor != "(Intercept)", ]
  
  cat("Survey-Weighted Logistic Regression for RIR\n")
  cat("Population: Statin users with LDL-C <70 mg/dL (N=577, 141 events)\n")
  cat(strrep("-", 70), "\n")
  print(table3, row.names = FALSE)
  
  write_csv(table3, file.path(output_dir, "Table3_Regression.csv"))
  cat("\nSaved: Table3_Regression.csv\n")
}

# =============================================================================
# SUPPLEMENTAL TABLE: Sensitivity Analyses
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("SUPPLEMENTAL TABLE: SENSITIVITY ANALYSES\n")
cat(strrep("=", 70), "\n\n")

# hs-CRP thresholds
cat("A. hs-CRP Threshold Sensitivity\n")
cat(strrep("-", 50), "\n")

hscrp_thresholds <- c(1.5, 2.0, 2.5, 3.0, 4.0, 5.0)
sens_hscrp <- list()

for (thresh in hscrp_thresholds) {
  df_temp <- df %>%
    filter(statin_user == 1, ldl_under_70 == 1) %>%
    mutate(rir_temp = as.numeric(hscrp >= thresh))
  
  design_temp <- svydesign(
    id = ~psu, strata = ~strata, weights = ~adj_weight,
    data = df_temp, nest = TRUE
  )
  
  prev <- svymean(~rir_temp, design = design_temp, na.rm = TRUE)
  ci <- confint(prev)
  
  sens_hscrp[[length(sens_hscrp) + 1]] <- data.frame(
    Threshold = paste0("≥", thresh, " mg/L"),
    Prevalence = sprintf("%.1f%%", coef(prev) * 100),
    CI_95 = sprintf("(%.1f%% - %.1f%%)", ci[1] * 100, ci[2] * 100)
  )
}

sens_hscrp_df <- bind_rows(sens_hscrp)
print(sens_hscrp_df, row.names = FALSE)

# LDL thresholds
cat("\nB. LDL-C Threshold Sensitivity\n")
cat(strrep("-", 50), "\n")

ldl_thresholds <- c(55, 60, 65, 70, 80, 100)
sens_ldl <- list()

for (thresh in ldl_thresholds) {
  df_temp <- df %>%
    filter(statin_user == 1, ldl_calc < thresh) %>%
    mutate(rir_temp = as.numeric(hscrp >= 2))
  
  if (nrow(df_temp) >= 30) {
    design_temp <- svydesign(
      id = ~psu, strata = ~strata, weights = ~adj_weight,
      data = df_temp, nest = TRUE
    )
    
    prev <- svymean(~rir_temp, design = design_temp, na.rm = TRUE)
    ci <- confint(prev)
    
    sens_ldl[[length(sens_ldl) + 1]] <- data.frame(
      Threshold = paste0("<", thresh, " mg/dL"),
      N = nrow(df_temp),
      Prevalence = sprintf("%.1f%%", coef(prev) * 100),
      CI_95 = sprintf("(%.1f%% - %.1f%%)", ci[1] * 100, ci[2] * 100)
    )
  }
}

sens_ldl_df <- bind_rows(sens_ldl)
print(sens_ldl_df, row.names = FALSE)

# Combine and save
suppl_table <- list(
  hscrp_sensitivity = sens_hscrp_df,
  ldl_sensitivity = sens_ldl_df
)

write_csv(sens_hscrp_df, file.path(output_dir, "SuppTable_Sensitivity_hsCRP.csv"))
write_csv(sens_ldl_df, file.path(output_dir, "SuppTable_Sensitivity_LDL.csv"))
cat("\nSaved: SuppTable_Sensitivity_hsCRP.csv\n")
cat("Saved: SuppTable_Sensitivity_LDL.csv\n")

# =============================================================================
# DONE
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE GENERATION COMPLETE\n")
cat(strrep("=", 70), "\n\n")
cat("Files saved to:", output_dir, "\n")

