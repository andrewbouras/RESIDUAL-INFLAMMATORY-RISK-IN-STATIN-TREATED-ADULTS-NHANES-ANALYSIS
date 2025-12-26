# =============================================================================
# RIR Study - Statistical Analysis
# Survey-weighted prevalence and logistic regression
# =============================================================================

# Suppress package startup messages
suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
  library(tableone)
  library(broom)
})

# Handle singleton PSUs (common in NHANES subsamples)
options(survey.lonely.psu = "adjust")

# Set working directory - handle both RStudio and command line
if (interactive() && exists("rstudioapi") && rstudioapi::isAvailable()) {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # When run from command line
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) {
    project_root <- dirname(dirname(normalizePath(script_path)))
  } else {
    project_root <- getwd()
  }
}
setwd(project_root)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("RIR STUDY - STATISTICAL ANALYSIS\n")
cat(strrep("=", 60), "\n")

data_path <- file.path(project_root, "data", "processed", "rir_cohort.csv")

if (!file.exists(data_path)) {
  stop("Data file not found. Run Python scripts first.")
}

df <- read_csv(data_path, show_col_types = FALSE)
cat("Loaded", nrow(df), "participants\n")

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Filter to primary cohort (statin users with LDL <70)
df_primary <- df %>%
  filter(primary_cohort == 1) %>%
  drop_na(fasting_weight, psu, strata)

cat("Primary cohort (statin + LDL<70):", nrow(df_primary), "\n")
cat("RIR cases:", sum(df_primary$rir), "\n")

# Create factors
df_primary <- df_primary %>%
  mutate(
    sex_f = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    race_eth_f = case_when(
      race_eth == 1 ~ "Mexican American",
      race_eth == 2 ~ "Other Hispanic",
      race_eth == 3 ~ "NH White",
      race_eth == 4 ~ "NH Black",
      race_eth == 6 ~ "NH Asian",
      race_eth == 7 ~ "Other/Multi",
      TRUE ~ "Other"
    ) %>% factor(),
    smoking_f = factor(smoking_status, 
                       levels = c(0, 1, 2),
                       labels = c("Never", "Former", "Current")),
    bmi_cat_f = factor(bmi_cat),
    age_group = cut(age, breaks = c(20, 45, 65, 100),
                    labels = c("20-44", "45-64", "≥65")),
    rir_f = factor(rir, levels = c(0, 1), labels = c("No RIR", "RIR")),
    # Continuous versions for regression
    hdl_low = (hdl < 40) %>% as.integer(),
    tg_high = (trigly >= 150) %>% as.integer()
  )

# =============================================================================
# SURVEY DESIGN
# =============================================================================

# Number of cycles pooled (for weight adjustment)
# 2005-2010 = 3 cycles, 2015-2020 = 2 cycles... but not contiguous
# For simplicity: divide by 5 (number of 2-year cycles represented)
n_cycles <- 5

# Create survey design for all eligible statin users
df_statin <- df %>%
  filter(statin_eligible == 1) %>%
  drop_na(fasting_weight, psu, strata) %>%
  mutate(
    sex_f = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    race_eth_f = case_when(
      race_eth == 1 ~ "Mexican American",
      race_eth == 2 ~ "Other Hispanic", 
      race_eth == 3 ~ "NH White",
      race_eth == 4 ~ "NH Black",
      race_eth == 6 ~ "NH Asian",
      race_eth == 7 ~ "Other/Multi",
      TRUE ~ "Other"
    ) %>% factor(),
    smoking_f = factor(smoking_status,
                       levels = c(0, 1, 2),
                       labels = c("Never", "Former", "Current")),
    age_group = cut(age, breaks = c(20, 45, 65, 100),
                    labels = c("20-44", "45-64", "≥65")),
    rir_f = factor(rir, levels = c(0, 1), labels = c("No RIR", "RIR"))
  )

design_statin <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~I(fasting_weight / n_cycles),
  data = df_statin,
  nest = TRUE
)

# Design for primary cohort
design_primary <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~I(fasting_weight / n_cycles),
  data = df_primary,
  nest = TRUE
)

cat("\nSurvey designs created\n")
cat("Statin users weighted N:", round(sum(weights(design_statin)), 0), "\n")
cat("Primary cohort weighted N:", round(sum(weights(design_primary)), 0), "\n")

# =============================================================================
# TABLE 1: BASELINE CHARACTERISTICS BY RIR STATUS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("TABLE 1: CHARACTERISTICS BY RIR STATUS\n")
cat(strrep("=", 60), "\n")

table1_vars <- c(
  "age", "female", "race_eth_f", "bmi", "obesity",
  "hypertension", "diabetes", "smoking_f", "cvd_history",
  "ldl", "hdl", "trigly", "hscrp"
)

cat_vars <- c("female", "race_eth_f", "obesity", "hypertension", 
              "diabetes", "smoking_f", "cvd_history")

table1 <- svyCreateTableOne(
  vars = table1_vars,
  strata = "rir_f",
  data = design_primary,
  factorVars = cat_vars,
  test = TRUE
)

print(table1, smd = TRUE)

# Export
dir.create(file.path(project_root, "output", "tables"), 
           showWarnings = FALSE, recursive = TRUE)

table1_df <- print(table1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_df, 
          file.path(project_root, "output", "tables", "table1_rir_characteristics.csv"))

# =============================================================================
# TABLE 2: PREVALENCE OF RIR BY SUBGROUPS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("TABLE 2: RIR PREVALENCE BY SUBGROUPS\n")
cat(strrep("=", 60), "\n")

# Overall prevalence in primary cohort
overall_prev <- svymean(~rir, design_primary, na.rm = TRUE)
cat("\nOverall RIR prevalence:", 
    sprintf("%.1f%% (%.1f%%-%.1f%%)",
            100*coef(overall_prev),
            100*(coef(overall_prev) - 1.96*SE(overall_prev)),
            100*(coef(overall_prev) + 1.96*SE(overall_prev))), "\n")

# Function to calculate subgroup prevalence
calc_subgroup_prev <- function(design, var_name) {
  result <- svyby(~rir, as.formula(paste("~", var_name)), 
                  design, svymean, na.rm = TRUE)
  
  result %>%
    mutate(
      prev_pct = sprintf("%.1f", 100 * rir),
      ci_low = 100 * (rir - 1.96 * se),
      ci_high = 100 * (rir + 1.96 * se),
      ci = sprintf("(%.1f-%.1f)", ci_low, ci_high),
      subgroup_var = var_name,
      level = as.character(result[[1]])  # Convert first column to character
    ) %>%
    select(subgroup_var, level, prev_pct, ci)
}

# Calculate by subgroups
prev_age <- calc_subgroup_prev(design_primary, "age_group")
prev_sex <- calc_subgroup_prev(design_primary, "sex_f")
prev_race <- calc_subgroup_prev(design_primary, "race_eth_f")
prev_obesity <- calc_subgroup_prev(design_primary, "obesity")
prev_diabetes <- calc_subgroup_prev(design_primary, "diabetes")
prev_smoking <- calc_subgroup_prev(design_primary, "smoking_f")
prev_cvd <- calc_subgroup_prev(design_primary, "cvd_history")

# Combine
table2 <- bind_rows(
  tibble(subgroup_var = "Overall", level = "All", 
         prev_pct = sprintf("%.1f", 100*coef(overall_prev)),
         ci = sprintf("(%.1f-%.1f)", 
                      100*(coef(overall_prev) - 1.96*SE(overall_prev)),
                      100*(coef(overall_prev) + 1.96*SE(overall_prev)))),
  prev_age, prev_sex, prev_race, prev_obesity, 
  prev_diabetes, prev_smoking, prev_cvd
)

print(table2)

write_csv(table2, 
          file.path(project_root, "output", "tables", "table2_rir_prevalence_subgroups.csv"))

# =============================================================================
# TABLE 3: LOGISTIC REGRESSION FOR RIR PREDICTORS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("TABLE 3: PREDICTORS OF RIR (LOGISTIC REGRESSION)\n")
cat(strrep("=", 60), "\n")

# Model 1: Demographics only
model1 <- svyglm(rir ~ age + sex_f + race_eth_f,
                 design = design_primary,
                 family = quasibinomial())

cat("\nModel 1: Demographics\n")
print(tidy(model1, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term != "(Intercept)") %>%
        select(term, estimate, conf.low, conf.high, p.value))

# Model 2: Add metabolic factors
model2 <- svyglm(rir ~ age + sex_f + race_eth_f + bmi + diabetes + hypertension,
                 design = design_primary,
                 family = quasibinomial())

cat("\nModel 2: + Metabolic factors\n")
print(tidy(model2, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term != "(Intercept)") %>%
        select(term, estimate, conf.low, conf.high, p.value))

# Model 3: Full model
model3 <- svyglm(rir ~ age + sex_f + race_eth_f + bmi + diabetes + 
                   hypertension + smoking_f + trigly + hdl + cvd_history,
                 design = design_primary,
                 family = quasibinomial())

cat("\nModel 3: Full model\n")
tidy_model3 <- tidy(model3, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value)
print(tidy_model3)

# Export regression results
table3 <- bind_rows(
  tidy(model1, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(model = "Model 1"),
  tidy(model2, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(model = "Model 2"),
  tidy(model3, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(model = "Model 3")
) %>%
  mutate(
    or_ci = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    p_value = sprintf("%.3f", p.value)
  ) %>%
  select(model, term, or_ci, p_value)

write_csv(table3,
          file.path(project_root, "output", "tables", "table3_rir_predictors.csv"))

# =============================================================================
# SENSITIVITY ANALYSES
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SENSITIVITY ANALYSES\n")
cat(strrep("=", 60), "\n")

# Sensitivity 1: RIR with hs-CRP ≥3 mg/L
cat("\n1. RIR defined as hs-CRP ≥3 mg/L:\n")
prev_crp3 <- svymean(~rir_crp3, design_primary, na.rm = TRUE)
cat("   Prevalence:", sprintf("%.1f%% (%.1f%%-%.1f%%)",
                              100*coef(prev_crp3),
                              100*(coef(prev_crp3) - 1.96*SE(prev_crp3)),
                              100*(coef(prev_crp3) + 1.96*SE(prev_crp3))), "\n")

# Sensitivity 2: LDL <55 mg/dL cohort
cat("\n2. LDL <55 mg/dL cohort:\n")
df_ldl55 <- df %>%
  filter(cohort_ldl55 == 1) %>%
  drop_na(fasting_weight, psu, strata)

if (nrow(df_ldl55) >= 30) {
  design_ldl55 <- svydesign(
    id = ~psu,
    strata = ~strata,
    weights = ~I(fasting_weight / n_cycles),
    data = df_ldl55,
    nest = TRUE
  )
  
  prev_ldl55 <- svymean(~rir_ldl55, design_ldl55, na.rm = TRUE)
  cat("   N:", nrow(df_ldl55), "\n")
  cat("   RIR prevalence:", sprintf("%.1f%% (%.1f%%-%.1f%%)",
                                    100*coef(prev_ldl55),
                                    100*(coef(prev_ldl55) - 1.96*SE(prev_ldl55)),
                                    100*(coef(prev_ldl55) + 1.96*SE(prev_ldl55))), "\n")
} else {
  cat("   Insufficient sample size (n=", nrow(df_ldl55), ")\n")
}

# Sensitivity 3: Restrict to clinical ASCVD
cat("\n3. ASCVD subgroup only:\n")
df_ascvd <- df %>%
  filter(primary_cohort == 1, cvd_history == 1) %>%
  drop_na(fasting_weight, psu, strata)

if (nrow(df_ascvd) >= 30) {
  design_ascvd <- svydesign(
    id = ~psu,
    strata = ~strata,
    weights = ~I(fasting_weight / n_cycles),
    data = df_ascvd,
    nest = TRUE
  )
  
  prev_ascvd <- svymean(~rir, design_ascvd, na.rm = TRUE)
  cat("   N:", nrow(df_ascvd), "\n")
  cat("   RIR prevalence:", sprintf("%.1f%% (%.1f%%-%.1f%%)",
                                    100*coef(prev_ascvd),
                                    100*(coef(prev_ascvd) - 1.96*SE(prev_ascvd)),
                                    100*(coef(prev_ascvd) + 1.96*SE(prev_ascvd))), "\n")
} else {
  cat("   Insufficient sample size (n=", nrow(df_ascvd), ")\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 60), "\n")

cat("\nKey findings:\n")
cat("1. Primary cohort (statin + LDL<70): n =", nrow(df_primary), "\n")
cat("2. RIR cases (hs-CRP ≥2): n =", sum(df_primary$rir), "\n")
cat("3. Weighted RIR prevalence:", sprintf("%.1f%%", 100*coef(overall_prev)), "\n")

cat("\nOutput files:\n")
cat("  - table1_rir_characteristics.csv\n")
cat("  - table2_rir_prevalence_subgroups.csv\n")
cat("  - table3_rir_predictors.csv\n")

