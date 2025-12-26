# =============================================================================
# RIR Study - Figures
# Forest plot and prevalence bar chart
# =============================================================================

library(tidyverse)
library(survey)
library(ggplot2)

# Set working directory
project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(project_root)

# Output directories
output_dir <- file.path(project_root, "output", "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================

data_path <- file.path(project_root, "data", "processed", "rir_cohort.csv")
df <- read_csv(data_path, show_col_types = FALSE)

# Filter to primary cohort
df_primary <- df %>%
  filter(primary_cohort == 1) %>%
  drop_na(fasting_weight, psu, strata)

# Create factors
df_primary <- df_primary %>%
  mutate(
    sex_f = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    age_group = cut(age, breaks = c(20, 45, 65, 100),
                    labels = c("20-44", "45-64", "≥65")),
    obesity_f = factor(obesity, levels = c(0, 1), labels = c("Non-obese", "Obese")),
    diabetes_f = factor(diabetes, levels = c(0, 1), labels = c("No Diabetes", "Diabetes")),
    smoking_f = factor(smoking_status, levels = c(0, 1, 2),
                       labels = c("Never", "Former", "Current")),
    cvd_f = factor(cvd_history, levels = c(0, 1), labels = c("No CVD", "CVD History"))
  )

# Survey design
n_cycles <- 5
design_primary <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~I(fasting_weight / n_cycles),
  data = df_primary,
  nest = TRUE
)

# =============================================================================
# FIGURE 1: BAR CHART OF RIR PREVALENCE BY SUBGROUPS
# =============================================================================

cat("\nGenerating Figure 1: RIR Prevalence by Subgroups\n")

# Calculate prevalence by subgroups
calc_prev <- function(design, var_name) {
  result <- svyby(~rir, as.formula(paste("~", var_name)), 
                  design, svymean, na.rm = TRUE)
  result %>%
    mutate(
      prevalence = rir * 100,
      se = se * 100,
      ci_low = prevalence - 1.96 * se,
      ci_high = prevalence + 1.96 * se,
      group_var = var_name
    ) %>%
    rename(group_level = 1) %>%
    select(group_var, group_level, prevalence, ci_low, ci_high)
}

# Calculate for each subgroup
prev_overall <- tibble(
  group_var = "Overall",
  group_level = "All",
  prevalence = 100 * coef(svymean(~rir, design_primary, na.rm = TRUE)),
  ci_low = prevalence - 1.96 * 100 * SE(svymean(~rir, design_primary, na.rm = TRUE)),
  ci_high = prevalence + 1.96 * 100 * SE(svymean(~rir, design_primary, na.rm = TRUE))
)

prev_data <- bind_rows(
  prev_overall,
  calc_prev(design_primary, "sex_f"),
  calc_prev(design_primary, "age_group"),
  calc_prev(design_primary, "obesity_f"),
  calc_prev(design_primary, "diabetes_f"),
  calc_prev(design_primary, "cvd_f")
)

# Clean up labels
prev_data <- prev_data %>%
  mutate(
    group_label = case_when(
      group_var == "Overall" ~ "Overall",
      group_var == "sex_f" ~ paste0("Sex: ", group_level),
      group_var == "age_group" ~ paste0("Age: ", group_level),
      group_var == "obesity_f" ~ group_level,
      group_var == "diabetes_f" ~ group_level,
      group_var == "cvd_f" ~ group_level
    ),
    group_category = case_when(
      group_var == "Overall" ~ "Overall",
      group_var == "sex_f" ~ "Sex",
      group_var == "age_group" ~ "Age Group",
      group_var %in% c("obesity_f", "diabetes_f", "cvd_f") ~ "Comorbidity"
    )
  )

# Order for plotting
prev_data <- prev_data %>%
  mutate(group_label = factor(group_label, levels = rev(unique(group_label))))

# Create bar chart
fig1 <- ggplot(prev_data, aes(x = group_label, y = prevalence)) +
  geom_col(aes(fill = group_category), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  geom_hline(yintercept = prev_overall$prevalence, linetype = "dashed", 
             color = "red", alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Overall" = "#2E86AB",
    "Sex" = "#A23B72",
    "Age Group" = "#F18F01",
    "Comorbidity" = "#C73E1D"
  )) +
  labs(
    title = "Prevalence of Residual Inflammatory Risk",
    subtitle = "Among statin-treated adults with LDL-C <70 mg/dL",
    x = NULL,
    y = "Prevalence (%)",
    fill = "Subgroup"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 60, 20))

ggsave(file.path(output_dir, "fig1_rir_prevalence_subgroups.pdf"), 
       fig1, width = 10, height = 7, dpi = 300)
ggsave(file.path(output_dir, "fig1_rir_prevalence_subgroups.png"), 
       fig1, width = 10, height = 7, dpi = 300)

cat("Figure 1 saved\n")

# =============================================================================
# FIGURE 2: FOREST PLOT OF ADJUSTED ORs
# =============================================================================

cat("\nGenerating Figure 2: Forest Plot of Adjusted ORs\n")

# Fit logistic regression
model <- svyglm(rir ~ age + sex_f + obesity_f + diabetes_f + 
                  trigly + hdl + cvd_f,
                design = design_primary,
                family = quasibinomial())

# Extract coefficients
forest_data <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = case_when(
      term == "age" ~ "Age (per year)",
      term == "sex_fFemale" ~ "Female vs Male",
      term == "obesity_fObese" ~ "Obese vs Non-obese",
      term == "diabetes_fDiabetes" ~ "Diabetes vs No Diabetes",
      term == "trigly" ~ "Triglycerides (per mg/dL)",
      term == "hdl" ~ "HDL-C (per mg/dL)",
      term == "cvd_fCVD History" ~ "CVD History vs None",
      TRUE ~ term
    ),
    significant = p.value < 0.05
  ) %>%
  arrange(desc(estimate))

# Forest plot
fig2 <- ggplot(forest_data, aes(x = estimate, y = reorder(term_label, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = significant), size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = significant), 
                 height = 0.2) +
  scale_color_manual(values = c("TRUE" = "#C73E1D", "FALSE" = "#2E86AB"),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  labs(
    title = "Predictors of Residual Inflammatory Risk",
    subtitle = "Adjusted Odds Ratios from Survey-Weighted Logistic Regression",
    x = "Odds Ratio (95% CI)",
    y = NULL,
    color = "Significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "fig2_forest_plot_ors.pdf"), 
       fig2, width = 9, height = 6, dpi = 300)
ggsave(file.path(output_dir, "fig2_forest_plot_ors.png"), 
       fig2, width = 9, height = 6, dpi = 300)

cat("Figure 2 saved\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("FIGURES GENERATED\n")
cat(strrep("=", 60), "\n")
cat("\nOutput files:\n")
cat("  - fig1_rir_prevalence_subgroups.pdf/png\n")
cat("  - fig2_forest_plot_ors.pdf/png\n")
cat("\nSaved to:", output_dir, "\n")

