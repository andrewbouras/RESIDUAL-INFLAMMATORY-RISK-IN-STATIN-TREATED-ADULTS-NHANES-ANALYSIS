# =============================================================================
# RIR: KAPLAN-MEIER SURVIVAL CURVES
# =============================================================================
# Shows all-cause mortality by hs-CRP status among statin users
# N=1,714 statin users, 488 deaths, mean 8.3 years follow-up
# =============================================================================

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# Color palette (consistent with RIR manuscript - navy blue)
rir_colors <- list(
  primary = "#1E4A6E",      # Deep navy blue (hs-CRP >=2)
  secondary = "#5B8FA8",    # Medium steel blue (hs-CRP <2)
  light = "#A3C4D9"
)

# Paths
data_path <- file.path("data", "processed", "rir_mortality_cohort.csv")
output_dir <- file.path("manuscript", "figures")

cat("Loading data...\n")
df <- read.csv(data_path)

# Prepare data
df_km <- df %>%
  mutate(
    hscrp_elevated = factor(ifelse(hscrp >= 2, "hs-CRP >=2 mg/L", "hs-CRP <2 mg/L"),
                            levels = c("hs-CRP <2 mg/L", "hs-CRP >=2 mg/L")),
    event = as.numeric(mort_all),
    time = followup_years
  ) %>%
  filter(!is.na(time) & !is.na(event) & !is.na(hscrp_elevated))

cat("Analytic sample N:", nrow(df_km), "\n")
cat("Deaths:", sum(df_km$event), "\n")
cat("By hs-CRP status:\n")
print(table(df_km$hscrp_elevated, df_km$event))

# =============================================================================
# FIT KAPLAN-MEIER MODEL
# =============================================================================

fit <- survfit(Surv(time, event) ~ hscrp_elevated, data = df_km)

cat("\nSurvival summary:\n")
print(fit)

# Cox model for HR
cox_model <- coxph(Surv(time, event) ~ hscrp_elevated + age + sex + race_eth, data = df_km)
cat("\nCox model (adjusted for age, sex, race):\n")
print(summary(cox_model)$conf.int)

# =============================================================================
# CREATE KAPLAN-MEIER PLOT
# =============================================================================

cat("\nCreating Kaplan-Meier figure...\n")

# Custom theme
km_theme <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0),
    plot.caption = element_text(size = 8, color = "gray50", hjust = 0),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Colors for the two groups
plot_colors <- c(rir_colors$secondary, rir_colors$primary)

# Use ggsurvplot
p_km <- ggsurvplot(
  fit,
  data = df_km,
  palette = plot_colors,
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  pval = TRUE,
  pval.coord = c(1, 0.25),
  pval.size = 4,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  ncensor.plot = FALSE,
  legend.title = "",
  legend.labs = c("hs-CRP <2 mg/L", "hs-CRP >=2 mg/L"),
  xlab = "Follow-up Time (Years)",
  ylab = "Survival Probability",
  title = "All-Cause Mortality by hs-CRP Status Among Statin Users",
  subtitle = "NHANES 2005-2016 with mortality follow-up | Adjusted HR = 1.49 (1.03-2.15), P=0.032",
  ggtheme = km_theme,
  surv.median.line = "none",
  break.time.by = 2
)

# Add caption
p_km$plot <- p_km$plot +
  labs(caption = "Statin users (N=1,714) with up to 14 years mortality follow-up. 488 deaths occurred.\nElevated hs-CRP (>=2 mg/L) associated with 49% higher mortality risk after adjustment for age, sex, and race/ethnicity.")

# Save
ggsave(file.path(output_dir, "Figure6_KM_Survival.png"), 
       print(p_km), width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6_KM_Survival.pdf"), 
       print(p_km), width = 8, height = 7, bg = "white")

cat("\n✅ Saved: Figure6_KM_Survival.png/pdf\n")

# =============================================================================
# ADDITIONAL: RIR SUBGROUP (LDL<70) SURVIVAL
# =============================================================================

cat("\n--- RIR Subgroup Analysis (LDL<70) ---\n")

df_rir <- df_km %>%
  filter(ldl_calc < 70) %>%
  mutate(
    rir_status = factor(ifelse(hscrp >= 2, "RIR (hs-CRP >=2)", "No RIR (hs-CRP <2)"),
                        levels = c("No RIR (hs-CRP <2)", "RIR (hs-CRP >=2)"))
  )

cat("LDL<70 subgroup N:", nrow(df_rir), "\n")
cat("Deaths:", sum(df_rir$event), "\n")
print(table(df_rir$rir_status, df_rir$event))

if (nrow(df_rir) >= 50) {
  fit_rir <- survfit(Surv(time, event) ~ rir_status, data = df_rir)
  
  p_rir <- ggsurvplot(
    fit_rir,
    data = df_rir,
    palette = plot_colors,
    conf.int = TRUE,
    conf.int.alpha = 0.15,
    pval = TRUE,
    pval.coord = c(1, 0.4),
    pval.size = 4,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    legend.title = "",
    legend.labs = c("No RIR (hs-CRP <2)", "RIR (hs-CRP >=2)"),
    xlab = "Follow-up Time (Years)",
    ylab = "Survival Probability",
    title = "All-Cause Mortality by RIR Status (Statin Users with LDL-C <70 mg/dL)",
    subtitle = "NHANES 2005-2016 | HR = 1.80 (0.93-3.49), P=0.08",
    ggtheme = km_theme,
    break.time.by = 2
  )
  
  p_rir$plot <- p_rir$plot +
    labs(caption = "Statin users with LDL-C <70 mg/dL (N=343). RIR defined as hs-CRP >=2 mg/L.\nTrend toward higher mortality with RIR; underpowered for significance.")
  
  ggsave(file.path(output_dir, "Figure7_KM_RIR_Subgroup.png"), 
         print(p_rir), width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "Figure7_KM_RIR_Subgroup.pdf"), 
         print(p_rir), width = 8, height = 7, bg = "white")
  
  cat("\n✅ Saved: Figure7_KM_RIR_Subgroup.png/pdf\n")
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n--- Mortality Summary ---\n")
df_km %>%
  group_by(hscrp_elevated) %>%
  summarise(
    N = n(),
    Deaths = sum(event),
    Mortality_Rate = Deaths / N * 100,
    Mean_FU_years = mean(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

