# =============================================================================
# RIR Study - Survey-Weighted Analysis
# =============================================================================
# 
# Performs all weighted analyses for the Residual Inflammatory Risk study:
# 1. Weighted prevalence of RIR among statin users with LDL<70
# 2. Table 1: Characteristics by RIR status
# 3. Subgroup prevalence (age, sex, race/ethnicity)
# 4. Survey-weighted logistic regression for RIR predictors
# 5. Sensitivity analyses (hs-CRP ≥3, LDL <55)
#
# Uses NHANES fasting subsample weights (WTSAF2YR)
# =============================================================================

library(survey)
library(tidyverse)
library(haven)
# library(tableone)  # Optional - for Table 1 creation

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
# Get project root - works in RStudio or command line
project_root <- tryCatch({
  dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
}, error = function(e) {
  getwd()
})

# Fallback if not in RStudio or result is empty
if (is.null(project_root) || project_root == "" || project_root == ".") {
  # Try to detect from script path
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) {
    project_root <- dirname(dirname(normalizePath(script_path)))
  } else {
    project_root <- getwd()
  }
}

data_path <- file.path(project_root, "data", "processed", "rir_analytic_cohort.csv")
output_dir <- file.path(project_root, "output", "tables")

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Handle lonely PSUs (common in NHANES subgroup analyses)
options(survey.lonely.psu = "adjust")

cat(strrep("=", 70), "\n")
cat("RIR STUDY - SURVEY-WEIGHTED ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
cat("[Loading data]\n")

if (!file.exists(data_path)) {
  stop("Data file not found. Run 02_process_data.py first.\n", data_path)
}

df <- read_csv(data_path, show_col_types = FALSE)
cat("  Loaded", nrow(df), "participants\n")

# -----------------------------------------------------------------------------
# Prepare variables for analysis
# -----------------------------------------------------------------------------
cat("[Preparing variables]\n")

# Create factor variables with labels
df <- df %>%
  mutate(
    # Sex
    sex_cat = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    
    # Race/ethnicity (RIDRETH1: 1=Mexican Am, 2=Other Hisp, 3=NH White, 4=NH Black, 5=Other)
    race_cat = factor(race_eth, 
                      levels = c(1, 2, 3, 4, 5),
                      labels = c("Mexican American", "Other Hispanic", 
                                "Non-Hispanic White", "Non-Hispanic Black", 
                                "Other/Multi-racial")),
    
    # Age groups
    age_group = cut(age, 
                    breaks = c(20, 40, 50, 60, 70, Inf),
                    labels = c("20-39", "40-49", "50-59", "60-69", "70+"),
                    right = FALSE),
    
    # BMI categories
    bmi_cat = cut(bmi,
                  breaks = c(0, 18.5, 25, 30, Inf),
                  labels = c("Underweight", "Normal", "Overweight", "Obese"),
                  right = FALSE),
    
    # Smoking status
    smoking_cat = factor(smoking_status,
                         levels = c(1, 2, 3),
                         labels = c("Current", "Former", "Never")),
    
    # Education (1=<9th, 2=9-11th, 3=HS/GED, 4=Some college, 5=College grad+)
    education_cat = factor(
      case_when(
        education %in% c(1, 2) ~ 1,
        education == 3 ~ 2,
        education == 4 ~ 3,
        education == 5 ~ 4,
        TRUE ~ NA_real_
      ),
      levels = c(1, 2, 3, 4),
      labels = c("Less than HS", "HS/GED", "Some college", "College graduate")
    ),
    
    # Binary versions for analysis
    diabetes_bin = as.numeric(diabetes == 1),
    hypertension_bin = as.numeric(hypertension == 1),
    obese_bin = as.numeric(obese == 1),
    current_smoker_bin = as.numeric(current_smoker == 1),
    
    # RIR status factor
    rir_cat = factor(rir, levels = c(0, 1), labels = c("No RIR", "RIR"))
  )

# -----------------------------------------------------------------------------
# Create Survey Design Objects
# -----------------------------------------------------------------------------
cat("[Creating survey design]\n")

# Full cohort - using fasting weights for lipid-based analyses
# For multi-cycle analysis, we adjust weights
# NHANES recommends dividing by number of cycles when combining

n_cycles <- length(unique(df$cycle))
df$adj_fasting_weight <- df$fasting_weight / n_cycles

# Full cohort design
design_full <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~adj_fasting_weight,
  data = df,
  nest = TRUE
)

# Subset: Statin users only
design_statin <- subset(design_full, statin_user == 1)

# Subset: Statin users with LDL <70 (primary analysis population)
design_primary <- subset(design_full, statin_user == 1 & ldl_under_70 == 1)

# Subset: Statin users with LDL <55 (sensitivity)
design_sensitivity <- subset(design_full, statin_user == 1 & ldl_under_55 == 1)

cat("  Full cohort N:", nrow(df), "\n")
cat("  Statin users N:", sum(df$statin_user == 1, na.rm = TRUE), "\n")
cat("  Statin + LDL<70 N:", sum(df$statin_user == 1 & df$ldl_under_70 == 1, na.rm = TRUE), "\n")

# =============================================================================
# ANALYSIS 1: Overall RIR Prevalence
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 1: RIR PREVALENCE\n")
cat(strrep("=", 70), "\n\n")

# Primary analysis: RIR prevalence among statin users with LDL<70
primary_cohort <- df %>% filter(statin_user == 1, ldl_under_70 == 1)

if (nrow(primary_cohort) > 0) {
  # Weighted prevalence
  rir_prev <- svymean(~rir, design = design_primary, na.rm = TRUE)
  rir_ci <- confint(rir_prev)
  
  cat("PRIMARY ANALYSIS: Statin users with LDL-C <70 mg/dL\n")
  cat("  N =", nrow(primary_cohort), "\n")
  cat("  RIR prevalence (hs-CRP ≥2 mg/L):\n")
  cat(sprintf("    %.1f%% (95%% CI: %.1f%% - %.1f%%)\n", 
              coef(rir_prev)["rir"] * 100,
              rir_ci["rir", 1] * 100,
              rir_ci["rir", 2] * 100))
  
  # Sensitivity: hs-CRP ≥3 mg/L
  rir_strict_prev <- svymean(~rir_strict, design = design_primary, na.rm = TRUE)
  rir_strict_ci <- confint(rir_strict_prev)
  
  cat("\n  RIR prevalence (hs-CRP ≥3 mg/L, sensitivity):\n")
  cat(sprintf("    %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
              coef(rir_strict_prev)["rir_strict"] * 100,
              rir_strict_ci["rir_strict", 1] * 100,
              rir_strict_ci["rir_strict", 2] * 100))
}

# Secondary: LDL <55 threshold
secondary_cohort <- df %>% filter(statin_user == 1, ldl_under_55 == 1)
if (nrow(secondary_cohort) > 10) {
  rir_ldl55_prev <- svymean(~rir_ldl55, design = design_sensitivity, na.rm = TRUE)
  rir_ldl55_ci <- confint(rir_ldl55_prev)
  
  cat("\nSENSITIVITY: Statin users with LDL-C <55 mg/dL\n")
  cat("  N =", nrow(secondary_cohort), "\n")
  cat(sprintf("  RIR prevalence: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
              coef(rir_ldl55_prev)["rir_ldl55"] * 100,
              rir_ldl55_ci["rir_ldl55", 1] * 100,
              rir_ldl55_ci["rir_ldl55", 2] * 100))
}

# =============================================================================
# ANALYSIS 2: RIR Prevalence by Subgroups
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 2: RIR PREVALENCE BY SUBGROUPS\n")
cat(strrep("=", 70), "\n\n")

# Function to calculate subgroup prevalence
calc_subgroup_prev <- function(design, formula, by_var) {
  # Get prevalence by subgroup
  prev <- svyby(formula, by_var, design = design, FUN = svymean, na.rm = TRUE)
  return(prev)
}

# By sex
cat("BY SEX:\n")
rir_by_sex <- svyby(~rir, ~sex_cat, design = design_primary, FUN = svymean, na.rm = TRUE)
print(rir_by_sex)

# By age group
cat("\nBY AGE GROUP:\n")
rir_by_age <- svyby(~rir, ~age_group, design = design_primary, FUN = svymean, na.rm = TRUE)
print(rir_by_age)

# By race/ethnicity
cat("\nBY RACE/ETHNICITY:\n")
rir_by_race <- svyby(~rir, ~race_cat, design = design_primary, FUN = svymean, na.rm = TRUE)
print(rir_by_race)

# =============================================================================
# ANALYSIS 3: Table 1 - Characteristics by RIR Status
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 3: CHARACTERISTICS BY RIR STATUS\n")
cat(strrep("=", 70), "\n\n")

# Variables for Table 1
table1_vars <- c("age", "sex_cat", "race_cat", "bmi", "bmi_cat",
                 "ldl_calc", "hdl", "trigly", "total_chol", "hscrp",
                 "diabetes_bin", "hypertension_bin", "current_smoker_bin",
                 "sbp_mean", "dbp_mean", "hba1c", "fasting_glucose")

# Create Table 1 comparing RIR vs No RIR
# Use only primary cohort
primary_df <- df %>% filter(statin_user == 1, ldl_under_70 == 1)

if (nrow(primary_df) > 0) {
  # Weighted means/proportions
  cat("Weighted characteristics (statin users with LDL<70):\n\n")
  
  # Continuous variables - means by RIR status
  continuous_vars <- c("age", "bmi", "ldl_calc", "hdl", "trigly", "hscrp",
                       "sbp_mean", "dbp_mean", "hba1c")
  
  for (var in continuous_vars) {
    if (var %in% names(primary_df)) {
      formula <- as.formula(paste0("~", var))
      by_formula <- ~rir_cat
      
      tryCatch({
        means <- svyby(formula, by_formula, design = design_primary, 
                       FUN = svymean, na.rm = TRUE)
        cat(sprintf("%s:\n", var))
        cat(sprintf("  No RIR: %.2f (SE: %.2f)\n", 
                    means[1, var], means[1, paste0("se.", var)]))
        cat(sprintf("  RIR:    %.2f (SE: %.2f)\n",
                    means[2, var], means[2, paste0("se.", var)]))
        
        # Survey-adjusted t-test
        test <- svyttest(formula, by = ~rir, design = design_primary)
        cat(sprintf("  p-value: %.4f\n\n", test$p.value))
      }, error = function(e) {
        cat(sprintf("%s: Could not calculate - %s\n\n", var, e$message))
      })
    }
  }
  
  # Categorical variables - proportions by RIR status
  cat("\nCategorical variables:\n")
  
  cat("\nSex (% Female):\n")
  tryCatch({
    sex_tab <- svyby(~I(sex == 2), ~rir_cat, design = design_primary, 
                     FUN = svymean, na.rm = TRUE)
    print(sex_tab)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  cat("\nDiabetes (%):\n")
  tryCatch({
    dm_tab <- svyby(~diabetes_bin, ~rir_cat, design = design_primary, 
                    FUN = svymean, na.rm = TRUE)
    print(dm_tab)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  cat("\nHypertension (%):\n")
  tryCatch({
    htn_tab <- svyby(~hypertension_bin, ~rir_cat, design = design_primary, 
                     FUN = svymean, na.rm = TRUE)
    print(htn_tab)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  cat("\nCurrent Smoker (%):\n")
  tryCatch({
    smoke_tab <- svyby(~current_smoker_bin, ~rir_cat, design = design_primary, 
                       FUN = svymean, na.rm = TRUE)
    print(smoke_tab)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  cat("\nObesity (BMI ≥30) (%):\n")
  tryCatch({
    obese_tab <- svyby(~obese_bin, ~rir_cat, design = design_primary, 
                       FUN = svymean, na.rm = TRUE)
    print(obese_tab)
  }, error = function(e) cat("  Error:", e$message, "\n"))
}

# =============================================================================
# ANALYSIS 4: Logistic Regression - Predictors of RIR
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 4: PREDICTORS OF RIR (Logistic Regression)\n")
cat(strrep("=", 70), "\n\n")

# Fit survey-weighted logistic regression
cat("Survey-weighted logistic regression for RIR\n")
cat("Population: Statin users with LDL<70 mg/dL\n\n")

tryCatch({
  # Model 1: Demographics only
  model1 <- svyglm(rir ~ age + sex_cat + race_cat, 
                   design = design_primary, 
                   family = quasibinomial())
  
  cat("MODEL 1: Demographics\n")
  cat(strrep("-", 50), "\n")
  summary_m1 <- summary(model1)
  print(summary_m1$coefficients)
  
  # Odds ratios
  cat("\nOdds Ratios:\n")
  or_m1 <- exp(coef(model1))
  ci_m1 <- exp(confint(model1))
  or_table1 <- data.frame(
    OR = or_m1,
    CI_lower = ci_m1[, 1],
    CI_upper = ci_m1[, 2]
  )
  print(round(or_table1, 3))
  
  # Model 2: Add clinical factors
  cat("\n\nMODEL 2: Demographics + Clinical Factors\n")
  cat(strrep("-", 50), "\n")
  
  model2 <- svyglm(rir ~ age + sex_cat + race_cat + bmi + 
                     diabetes_bin + hypertension_bin + current_smoker_bin + trigly,
                   design = design_primary,
                   family = quasibinomial())
  
  summary_m2 <- summary(model2)
  print(summary_m2$coefficients)
  
  cat("\nOdds Ratios:\n")
  or_m2 <- exp(coef(model2))
  ci_m2 <- exp(confint(model2))
  or_table2 <- data.frame(
    OR = or_m2,
    CI_lower = ci_m2[, 1],
    CI_upper = ci_m2[, 2]
  )
  print(round(or_table2, 3))
  
}, error = function(e) {
  cat("Error fitting models:", e$message, "\n")
})

# =============================================================================
# ANALYSIS 5: Sensitivity Analyses
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 5: SENSITIVITY ANALYSES\n")
cat(strrep("=", 70), "\n\n")

# A. Different hs-CRP thresholds
cat("A. hs-CRP THRESHOLD SENSITIVITY\n")
cat(strrep("-", 40), "\n")

thresholds <- c(1.5, 2, 2.5, 3, 4, 5)
for (thresh in thresholds) {
  df_temp <- df %>% 
    filter(statin_user == 1, ldl_under_70 == 1) %>%
    mutate(rir_temp = as.numeric(hscrp >= thresh))
  
  if (sum(df_temp$rir_temp, na.rm = TRUE) > 0) {
    design_temp <- svydesign(
      id = ~psu, strata = ~strata, weights = ~adj_fasting_weight,
      data = df_temp, nest = TRUE
    )
    prev <- svymean(~rir_temp, design = design_temp, na.rm = TRUE)
    ci <- confint(prev)
    cat(sprintf("hs-CRP ≥%.1f mg/L: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
                thresh, coef(prev) * 100, ci[1] * 100, ci[2] * 100))
  }
}

# B. Different LDL thresholds
cat("\nB. LDL-C THRESHOLD SENSITIVITY\n")
cat(strrep("-", 40), "\n")

ldl_thresholds <- c(55, 60, 65, 70, 80, 100)
for (ldl_thresh in ldl_thresholds) {
  df_temp <- df %>% 
    filter(statin_user == 1, ldl_calc < ldl_thresh) %>%
    mutate(rir_temp = as.numeric(hscrp >= 2))
  
  if (nrow(df_temp) > 30 && sum(df_temp$rir_temp, na.rm = TRUE) > 0) {
    design_temp <- svydesign(
      id = ~psu, strata = ~strata, weights = ~adj_fasting_weight,
      data = df_temp, nest = TRUE
    )
    n_cohort <- nrow(df_temp)
    prev <- svymean(~rir_temp, design = design_temp, na.rm = TRUE)
    ci <- confint(prev)
    cat(sprintf("LDL-C <%d mg/dL (N=%d): %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
                ldl_thresh, n_cohort, coef(prev) * 100, ci[1] * 100, ci[2] * 100))
  }
}

# =============================================================================
# SAVE RESULTS
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("SAVING RESULTS\n")
cat(strrep("=", 70), "\n\n")

# Create summary results dataframe
results_summary <- data.frame(
  Analysis = character(),
  N = numeric(),
  Prevalence_pct = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  stringsAsFactors = FALSE
)

# Add primary result
if (nrow(primary_cohort) > 0) {
  rir_prev <- svymean(~rir, design = design_primary, na.rm = TRUE)
  rir_ci <- confint(rir_prev)
  
  results_summary <- rbind(results_summary, data.frame(
    Analysis = "Primary: RIR (statin + LDL<70 + hs-CRP≥2)",
    N = nrow(primary_cohort),
    Prevalence_pct = round(coef(rir_prev)["rir"] * 100, 1),
    CI_lower = round(rir_ci["rir", 1] * 100, 1),
    CI_upper = round(rir_ci["rir", 2] * 100, 1)
  ))
  
  # Strict definition
  rir_strict_prev <- svymean(~rir_strict, design = design_primary, na.rm = TRUE)
  rir_strict_ci <- confint(rir_strict_prev)
  
  results_summary <- rbind(results_summary, data.frame(
    Analysis = "Sensitivity: RIR (hs-CRP≥3)",
    N = nrow(primary_cohort),
    Prevalence_pct = round(coef(rir_strict_prev)["rir_strict"] * 100, 1),
    CI_lower = round(rir_strict_ci["rir_strict", 1] * 100, 1),
    CI_upper = round(rir_strict_ci["rir_strict", 2] * 100, 1)
  ))
}

# Save results
write_csv(results_summary, file.path(output_dir, "rir_prevalence_results.csv"))
cat("Saved: rir_prevalence_results.csv\n")

# Save subgroup results
if (exists("rir_by_sex")) {
  rir_by_sex$subgroup <- "Sex"
  rir_by_age$subgroup <- "Age"
  rir_by_race$subgroup <- "Race"
  
  subgroup_results <- bind_rows(
    rir_by_sex %>% rename(category = sex_cat),
    rir_by_age %>% rename(category = age_group),
    rir_by_race %>% rename(category = race_cat)
  )
  
  write_csv(subgroup_results, file.path(output_dir, "rir_prevalence_by_subgroup.csv"))
  cat("Saved: rir_prevalence_by_subgroup.csv\n")
}

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")

