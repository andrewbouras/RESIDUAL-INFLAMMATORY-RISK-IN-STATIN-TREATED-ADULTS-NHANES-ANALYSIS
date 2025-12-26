# =============================================================================
# RIR: DOSE-RESPONSE CURVE - hs-CRP Threshold Sensitivity
# =============================================================================
# Shows RIR prevalence across different hs-CRP thresholds
# Demonstrates robustness of the primary finding
# =============================================================================

library(survey)
library(tidyverse)

options(survey.lonely.psu = "adjust")

# Color palette (consistent with RIR manuscript)
rir_colors <- list(
  primary = "#1E4A6E",
  secondary = "#5B8FA8",
  light = "#A3C4D9",
  accent = "#C44E52"
)

# Load data
data_path <- file.path("data", "processed", "rir_analytic_cohort.csv")
output_dir <- file.path("manuscript", "figures")

cat("Loading data...\n")
df <- read_csv(data_path, show_col_types = FALSE)

# Filter to primary cohort: statin users with LDL-C <70
primary_cohort <- df %>% 
  filter(statin_user == 1 & ldl_calc < 70)

cat("Primary cohort N:", nrow(primary_cohort), "\n")

# Adjust weights
n_cycles <- length(unique(primary_cohort$cycle))
primary_cohort$adj_weight <- primary_cohort$fasting_weight / n_cycles

# Create survey design
design_primary <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~adj_weight,
  data = primary_cohort,
  nest = TRUE
)

# =============================================================================
# Calculate prevalence at each hs-CRP threshold
# =============================================================================

thresholds <- c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)
results <- data.frame()

cat("\nCalculating prevalence by hs-CRP threshold...\n")

for (thresh in thresholds) {
  # Create binary variable for this threshold
  primary_cohort_temp <- primary_cohort %>%
    mutate(above_thresh = ifelse(hscrp >= thresh, 1, 0))
  
  design_temp <- svydesign(
    id = ~psu,
    strata = ~strata,
    weights = ~adj_weight,
    data = primary_cohort_temp,
    nest = TRUE
  )
  
  prev <- svymean(~above_thresh, design = design_temp, na.rm = TRUE)
  ci <- confint(prev)
  
  results <- rbind(results, data.frame(
    Threshold = thresh,
    Prevalence = coef(prev) * 100,
    CI_low = ci[1] * 100,
    CI_high = ci[2] * 100,
    SE = sqrt(attr(prev, "var")) * 100
  ))
  
  cat(sprintf("  hs-CRP >= %.1f: %.1f%% (%.1f-%.1f)\n", 
              thresh, coef(prev)*100, ci[1]*100, ci[2]*100))
}

# Mark the primary threshold (2.0)
results$Primary <- ifelse(results$Threshold == 2.0, TRUE, FALSE)

# =============================================================================
# Create dose-response plot
# =============================================================================

cat("\nCreating dose-response figure...\n")

# Theme consistent with other RIR figures
theme_rir <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "gray40", hjust = 0),
      plot.caption = element_text(size = base_size - 2, color = "gray50", hjust = 0),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 20, 15, 15)
    )
}

fig <- ggplot(results, aes(x = Threshold, y = Prevalence)) +
  # Confidence band
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), 
              fill = rir_colors$light, alpha = 0.5) +
  # Main line
  geom_line(color = rir_colors$primary, linewidth = 1.2) +
  # Points
  geom_point(aes(shape = Primary, size = Primary), color = rir_colors$primary) +
  # Highlight primary threshold
  geom_point(data = filter(results, Primary), 
             aes(x = Threshold, y = Prevalence),
             color = rir_colors$accent, size = 5, shape = 18) +
  # Add prevalence labels
  geom_text(aes(label = sprintf("%.1f%%", Prevalence)), 
            vjust = -1.5, size = 3, color = "gray30") +
  # Annotation for primary threshold
  annotate("text", x = 2.0, y = results$Prevalence[results$Threshold == 2.0] + 8,
           label = "Primary\nthreshold", size = 3, color = rir_colors$accent,
           fontface = "bold", lineheight = 0.9) +
  # Scales
  scale_x_continuous(breaks = thresholds, 
                     labels = paste0(">=", thresholds),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10),
                     expand = expansion(mult = c(0, 0.08))) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5), guide = "none") +
  # Labels
  labs(
    title = "RIR Prevalence Across hs-CRP Thresholds",
    subtitle = "Statin users with LDL-C <70 mg/dL (N = 577) | NHANES 2005-2010, 2015-2020",
    x = "hs-CRP Threshold (mg/L)",
    y = "Prevalence (%)",
    caption = "Shaded region represents 95% CI. Primary threshold (>=2.0 mg/L) highlighted in red."
  ) +
  theme_rir()

# Save
ggsave(file.path(output_dir, "Figure5_Dose_Response.png"), fig, 
       width = 9, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure5_Dose_Response.pdf"), fig, 
       width = 9, height = 5.5, bg = "white")

cat("\n✅ Saved: Figure5_Dose_Response.png/pdf\n")

# Print summary table
cat("\n--- hs-CRP Threshold Sensitivity ---\n")
print(results %>% select(Threshold, Prevalence, CI_low, CI_high) %>%
        mutate(across(where(is.numeric), ~round(., 1))))

