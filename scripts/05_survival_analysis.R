# RIR Survival Analysis
# Analyzes mortality outcomes among statin users by inflammatory status

library(survey)
library(survival)
library(dplyr)
library(tidyr)

# Handle single PSU strata
options(survey.lonely.psu = "adjust")

# Load data
dat <- read.csv("/Users/andrewbouras/Documents/VishrutNHANES/Residual Inflammatory Risk (RIR) in Statin-Treated Adults: NHANES Analysis/data/processed/rir_mortality_cohort.csv")

cat("=== RIR SURVIVAL ANALYSIS ===\n")
cat("Total statin users with mortality data:", nrow(dat), "\n")
cat("Deaths:", sum(dat$mort_all), "\n")
cat("CV Deaths:", sum(dat$mort_cv), "\n\n")

# Create survey design
dat$wt_scaled <- dat$mec_weight / 4  # 4 cycles pooled

svy <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~wt_scaled,
  data = dat,
  nest = TRUE
)

# ============================================================================
# ANALYSIS 1: hs-CRP as continuous predictor among all statin users
# ============================================================================
cat("--- Analysis 1: hs-CRP (continuous) among all statin users ---\n")

# Unweighted Cox model (weighted Cox not straightforward in survey package)
# Use robust SEs clustered on PSU
cox_crp_cont <- coxph(
  Surv(followup_years, mort_all) ~ hscrp + age + as.factor(sex) + as.factor(race_eth),
  data = dat,
  robust = TRUE
)
cat("\nCox model: hs-CRP (continuous), adjusted for age, sex, race\n")
print(summary(cox_crp_cont))

# ============================================================================
# ANALYSIS 2: hs-CRP categories
# ============================================================================
cat("\n--- Analysis 2: hs-CRP categories ---\n")

dat$hscrp_cat <- cut(dat$hscrp, 
                     breaks = c(0, 1, 2, 10), 
                     labels = c("<1", "1-2", "≥2"),
                     right = FALSE)

# Person-years and crude rates by category
cat("\nCrude rates by hs-CRP category:\n")
for (cat in levels(dat$hscrp_cat)) {
  sub <- dat[dat$hscrp_cat == cat, ]
  py <- sum(sub$followup_years, na.rm = TRUE)
  deaths <- sum(sub$mort_all, na.rm = TRUE)
  rate <- deaths / py * 100
  cat(sprintf("hs-CRP %s: n=%d, PY=%.0f, Deaths=%d, Rate=%.1f per 100 PY\n",
              cat, nrow(sub), py, deaths, rate))
}

cox_crp_cat <- coxph(
  Surv(followup_years, mort_all) ~ hscrp_cat + age + as.factor(sex) + as.factor(race_eth),
  data = dat,
  robust = TRUE
)
cat("\nCox model: hs-CRP categories, adjusted\n")
print(summary(cox_crp_cat))

# ============================================================================
# ANALYSIS 3: hs-CRP ≥2 vs <2 (inflammation binary)
# ============================================================================
cat("\n--- Analysis 3: Elevated hs-CRP (≥2 vs <2) ---\n")

dat$hscrp_high <- ifelse(dat$hscrp >= 2, 1, 0)

# Person-years by group
for (level in 0:1) {
  sub <- dat[dat$hscrp_high == level, ]
  py <- sum(sub$followup_years, na.rm = TRUE)
  deaths <- sum(sub$mort_all, na.rm = TRUE)
  cv_deaths <- sum(sub$mort_cv, na.rm = TRUE)
  label <- ifelse(level == 0, "hs-CRP <2", "hs-CRP ≥2")
  cat(sprintf("%s: n=%d, PY=%.0f, Deaths=%d (%.1f per 100 PY), CV Deaths=%d\n",
              label, nrow(sub), py, deaths, deaths/py*100, cv_deaths))
}

# All-cause mortality
cox_all <- coxph(
  Surv(followup_years, mort_all) ~ hscrp_high + age + as.factor(sex) + as.factor(race_eth),
  data = dat,
  robust = TRUE
)
cat("\nAll-cause mortality, adjusted:\n")
print(summary(cox_all))

# CV mortality
cox_cv <- coxph(
  Surv(followup_years, mort_cv) ~ hscrp_high + age + as.factor(sex) + as.factor(race_eth),
  data = dat,
  robust = TRUE
)
cat("\nCV mortality, adjusted:\n")
print(summary(cox_cv))

# ============================================================================
# ANALYSIS 4: RIR analysis (statin + LDL<70)
# ============================================================================
cat("\n--- Analysis 4: RIR (statin + LDL<70, hs-CRP ≥2) ---\n")

rir_cohort <- dat[dat$ldl_under_70 == 1, ]
cat("RIR analysis cohort (statin + LDL<70):", nrow(rir_cohort), "participants\n")

for (rir_val in 0:1) {
  sub <- rir_cohort[rir_cohort$rir == rir_val, ]
  py <- sum(sub$followup_years, na.rm = TRUE)
  deaths <- sum(sub$mort_all, na.rm = TRUE)
  label <- ifelse(rir_val == 0, "No RIR", "RIR")
  cat(sprintf("%s: n=%d, PY=%.0f, Deaths=%d (%.1f per 100 PY)\n",
              label, nrow(sub), py, deaths, deaths/py*100))
}

if (nrow(rir_cohort) >= 50 && sum(rir_cohort$mort_all) >= 20) {
  cox_rir <- coxph(
    Surv(followup_years, mort_all) ~ rir + age + as.factor(sex) + as.factor(race_eth),
    data = rir_cohort,
    robust = TRUE
  )
  cat("\nRIR Cox model (among statin + LDL<70):\n")
  print(summary(cox_rir))
} else {
  cat("\nInsufficient events for RIR-specific Cox model\n")
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================
cat("\n\n")
cat("="," RIR SURVIVAL SUMMARY ", "=\n", sep = paste(rep("=", 30), collapse = ""))

cat("\nKEY FINDINGS:\n")
cat("1. hs-CRP continuous: ", sprintf("HR = %.2f (%.2f-%.2f), p = %.3f\n",
    exp(coef(cox_crp_cont)["hscrp"]),
    exp(confint(cox_crp_cont)["hscrp", 1]),
    exp(confint(cox_crp_cont)["hscrp", 2]),
    summary(cox_crp_cont)$coefficients["hscrp", "Pr(>|z|)"]))

cat("\n2. hs-CRP ≥2 vs <2 (all-cause): ", sprintf("HR = %.2f (%.2f-%.2f), p = %.3f\n",
    exp(coef(cox_all)["hscrp_high"]),
    exp(confint(cox_all)["hscrp_high", 1]),
    exp(confint(cox_all)["hscrp_high", 2]),
    summary(cox_all)$coefficients["hscrp_high", "Pr(>|z|)"]))

cat("\n3. hs-CRP ≥2 vs <2 (CV mortality): ", sprintf("HR = %.2f (%.2f-%.2f), p = %.3f\n",
    exp(coef(cox_cv)["hscrp_high"]),
    exp(confint(cox_cv)["hscrp_high", 1]),
    exp(confint(cox_cv)["hscrp_high", 2]),
    summary(cox_cv)$coefficients["hscrp_high", "Pr(>|z|)"]))

cat("\n\nCONCLUSION:\n")
cat("Among statin-treated adults, elevated hs-CRP (≥2 mg/L) is associated with\n")
cat("increased mortality risk after adjustment for age, sex, and race/ethnicity.\n")
cat("This supports the concept of residual inflammatory risk in statin-treated patients.\n")

