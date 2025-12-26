# =============================================================================
# FIGURE: Clinical Characteristics by RIR Status
# Simplified visual counterpart to Table 1 showing key differences between
# RIR and Non-RIR groups among statin users with LDL-C <70 mg/dL
# =============================================================================

library(survey)
library(srvyr)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# Set survey option for lonely PSUs
options(survey.lonely.psu = "adjust")

# Define paths
data_path <- file.path("data", "processed", "rir_analytic_cohort.csv")
output_dir <- file.path("manuscript", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading data...\n")
df <- read_csv(data_path, show_col_types = FALSE)

# Filter to statin users with LDL-C <70 (primary cohort for RIR)
primary_cohort <- df %>%
  filter(statin_user == 1 & ldl_calc < 70) %>%
  mutate(
    rir_status = factor(ifelse(hscrp >= 2, "RIR", "No RIR"), 
                        levels = c("No RIR", "RIR"))
  )

cat("Primary cohort N:", nrow(primary_cohort), "\n")
cat("RIR cases:", sum(primary_cohort$rir_status == "RIR"), "\n")

# =============================================================================
# CALCULATE WEIGHTED MEANS AND PROPORTIONS BY RIR STATUS
# =============================================================================

# Continuous variables - weighted means
calc_weighted_mean <- function(data, var, group_var = "rir_status") {
  data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      mean = weighted.mean(.data[[var]], fasting_weight, na.rm = TRUE),
      sd = sqrt(sum(fasting_weight * (.data[[var]] - weighted.mean(.data[[var]], fasting_weight))^2) / sum(fasting_weight)),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      se = sd / sqrt(n),
      variable = var
    )
}

# Binary variables - weighted proportions
calc_weighted_prop <- function(data, var, group_var = "rir_status") {
  data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      prop = weighted.mean(.data[[var]], fasting_weight, na.rm = TRUE) * 100,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      se = sqrt(prop/100 * (1 - prop/100) / n) * 100,
      variable = var
    )
}

# Calculate for key variables
cat("Calculating weighted statistics...\n")

# Continuous variables
bmi_stats <- calc_weighted_mean(primary_cohort, "bmi")
hba1c_stats <- calc_weighted_mean(primary_cohort, "hba1c")
age_stats <- calc_weighted_mean(primary_cohort, "age")

# Binary variables
diabetes_stats <- calc_weighted_prop(primary_cohort, "diabetes")
hypertension_stats <- calc_weighted_prop(primary_cohort, "hypertension")
obese_stats <- calc_weighted_prop(primary_cohort, "obese")
smoker_stats <- calc_weighted_prop(primary_cohort, "current_smoker")

# Sex (female)
primary_cohort$female <- ifelse(primary_cohort$sex == 2, 1, 0)
female_stats <- calc_weighted_prop(primary_cohort, "female")

# =============================================================================
# PANEL A: KEY CONTINUOUS VARIABLES (BMI, HbA1c)
# =============================================================================

continuous_data <- bind_rows(
  bmi_stats %>% mutate(var_label = "BMI (kg/m²)"),
  hba1c_stats %>% mutate(var_label = "HbA1c (%)")
)

# Color scheme
rir_colors <- c("No RIR" = "#4A90A4", "RIR" = "#C44E52")

p_continuous <- ggplot(continuous_data, aes(x = var_label, y = mean, fill = rir_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se),
                position = position_dodge(width = 0.8), width = 0.2, color = "black") +
  geom_text(aes(label = sprintf("%.1f", mean), y = mean + 1.96*se + 1),
            position = position_dodge(width = 0.8), size = 3.5, vjust = 0) +
  scale_fill_manual(values = rir_colors, name = "") +
  labs(x = "", y = "Mean Value", 
       title = "A. Metabolic Characteristics") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_cartesian(ylim = c(0, 40))

# =============================================================================
# PANEL B: KEY COMORBIDITIES (Prevalence)
# =============================================================================

comorbidity_data <- bind_rows(
  diabetes_stats %>% mutate(var_label = "Diabetes"),
  hypertension_stats %>% mutate(var_label = "Hypertension"),
  obese_stats %>% mutate(var_label = "Obesity\n(BMI ≥30)"),
  smoker_stats %>% mutate(var_label = "Current\nSmoker")
)

p_comorbidities <- ggplot(comorbidity_data, aes(x = var_label, y = prop, fill = rir_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, prop - 1.96*se), ymax = pmin(100, prop + 1.96*se)),
                position = position_dodge(width = 0.8), width = 0.2, color = "black") +
  geom_text(aes(label = sprintf("%.0f%%", prop), y = pmin(100, prop + 1.96*se) + 3),
            position = position_dodge(width = 0.8), size = 3, vjust = 0) +
  scale_fill_manual(values = rir_colors, name = "") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "", y = "Prevalence (%)", 
       title = "B. Comorbidity Prevalence") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10)
  )

# =============================================================================
# COMBINE PANELS
# =============================================================================

cat("Combining panels...\n")

# Use patchwork if available, otherwise gridExtra
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  combined <- p_continuous + p_comorbidities +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Clinical Characteristics by Residual Inflammatory Risk Status",
      subtitle = "Statin users with LDL-C <70 mg/dL (N = 577)",
      caption = "Error bars represent 95% CI. RIR defined as hs-CRP ≥2 mg/L.",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 9, hjust = 0, color = "gray50")
      )
    )
} else {
  library(gridExtra)
  combined <- grid.arrange(p_continuous, p_comorbidities, ncol = 2,
                           top = "Clinical Characteristics by RIR Status")
}

# =============================================================================
# SAVE FIGURE
# =============================================================================

cat("Saving figure...\n")

ggsave(file.path(output_dir, "Figure4_Characteristics_by_RIR.png"), 
       plot = combined, width = 10, height = 5, dpi = 300)

ggsave(file.path(output_dir, "Figure4_Characteristics_by_RIR.pdf"), 
       plot = combined, width = 10, height = 5)

# Clean up old figure
if (file.exists(file.path(output_dir, "Figure_RIR_by_Characteristics.png"))) {
  file.remove(file.path(output_dir, "Figure_RIR_by_Characteristics.png"))
  file.remove(file.path(output_dir, "Figure_RIR_by_Characteristics.pdf"))
  cat("Removed old 10-panel figure\n")
}

cat("\n✅ Figure saved to:", file.path(output_dir, "Figure4_Characteristics_by_RIR.png"), "\n")

# =============================================================================
# PRINT SUMMARY STATISTICS
# =============================================================================

cat("\n--- Summary Statistics ---\n")
cat("\nContinuous Variables:\n")
print(continuous_data %>% select(rir_status, var_label, mean, se))

cat("\nComorbidity Prevalence:\n")
print(comorbidity_data %>% select(rir_status, var_label, prop, se))
