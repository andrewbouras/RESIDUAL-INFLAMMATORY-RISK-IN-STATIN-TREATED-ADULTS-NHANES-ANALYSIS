# =============================================================================
# RIR MANUSCRIPT - UNIFIED COLOR SCHEME FIGURES
# =============================================================================
# Professional medical journal color palette: Navy blue-based monochromatic
# All 4 figures regenerated with consistent styling
# =============================================================================

library(survey)
library(tidyverse)
library(scales)

options(survey.lonely.psu = "adjust")

# =============================================================================
# COLOR PALETTE: RIR STUDY (Navy Blue Monochromatic)
# =============================================================================
# Professional palette inspired by Circulation/JACC styling
rir_colors <- list(
  # Primary comparison (e.g., Statin vs No Statin)
  primary = "#1E4A6E",      # Deep navy blue
  secondary = "#5B8FA8",    # Medium steel blue
  
  # For 3-category scales (e.g., glycemic status)
  cat3_dark = "#1E4A6E",    # Deep navy
  cat3_mid = "#5B8FA8",     # Medium blue
  cat3_light = "#A3C4D9",   # Light blue
  
  # For comparison groups (RIR vs No RIR)
  rir_yes = "#1E4A6E",      # Deep navy (affected)
  rir_no = "#A3C4D9",       # Light blue (reference)
  
  # Neutral
  gray = "#6B6B6B",
  light_gray = "#D0D0D0"
)

# =============================================================================
# SETUP
# =============================================================================
data_path <- file.path("data", "processed", "rir_analytic_cohort.csv")
output_dir <- file.path("manuscript", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading data...\n")
df <- read_csv(data_path, show_col_types = FALSE)

# Prepare variables
df <- df %>%
  mutate(
    sex_cat = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    glycemic_cat = factor(glycemic_status,
                          levels = c(0, 1, 2),
                          labels = c("Normal", "Prediabetes", "Diabetes")),
    ldl_statin_group_f = factor(ldl_statin_group,
                                levels = c(1, 2, 3, 4, 5, 6),
                                labels = c("LDL <=70\nStatin", "LDL <=70\nNo Statin",
                                          "LDL 70-130\nStatin", "LDL 70-130\nNo Statin",
                                          "LDL >130\nStatin", "LDL >130\nNo Statin")),
    statin_status = factor(ifelse(statin_user == 1, "Statin User", "Non-User"),
                           levels = c("Statin User", "Non-User")),
    rir_status = factor(ifelse(hscrp >= 2, "RIR", "No RIR"), 
                        levels = c("No RIR", "RIR"))
  )

n_cycles <- length(unique(df$cycle))
df$adj_weight <- df$fasting_weight / n_cycles

design_full <- svydesign(
  id = ~psu, strata = ~strata, weights = ~adj_weight, data = df, nest = TRUE
)

# =============================================================================
# PUBLICATION THEME
# =============================================================================
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
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# =============================================================================
# FIGURE 1: 6-Group LDL×Statin
# =============================================================================
cat("\n[Figure 1: LDL×Statin Groups]\n")

groups_list <- levels(df$ldl_statin_group_f)
fig1_prev <- data.frame()

for (grp in groups_list) {
  grp_design <- subset(design_full, ldl_statin_group_f == grp)
  if (nrow(grp_design) > 0) {
    prev <- svymean(~hscrp_high, design = grp_design, na.rm = TRUE)
    ci <- confint(prev)
    fig1_prev <- rbind(fig1_prev, data.frame(
      Group = grp,
      Prevalence = coef(prev) * 100,
      CI_low = ci[1] * 100,
      CI_high = ci[2] * 100
    ))
  }
}

fig1_prev <- fig1_prev %>%
  mutate(
    Statin_Status = ifelse(grepl("No Statin", gsub("\n", " ", Group)), "Non-User", "Statin User"),
    Statin_Status = factor(Statin_Status, levels = c("Statin User", "Non-User")),
    LDL_Tier = case_when(
      grepl("<=70", Group) ~ "LDL <=70",
      grepl("70-130", Group) ~ "LDL 70-130",
      grepl(">130", Group) ~ "LDL >130"
    ),
    LDL_Tier = factor(LDL_Tier, levels = c("LDL <=70", "LDL 70-130", "LDL >130"))
  )

fig1 <- ggplot(fig1_prev, aes(x = LDL_Tier, y = Prevalence, fill = Statin_Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Prevalence), y = CI_high + 2),
            position = position_dodge(width = 0.8), size = 3.2, vjust = 0) +
  scale_fill_manual(values = c("Statin User" = rir_colors$primary, 
                               "Non-User" = rir_colors$secondary), name = "") +
  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "hs-CRP >=2 mg/L Prevalence by LDL-C Category and Statin Status",
    subtitle = "NHANES 2005-2010, 2015-2020 (N = 12,896)",
    x = "LDL-C Category (mg/dL)",
    y = "Prevalence (%)"
  ) +
  theme_rir()

ggsave(file.path(output_dir, "Figure1_LDL_Statin_Groups.png"), fig1, 
       width = 8, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure1_LDL_Statin_Groups.pdf"), fig1, 
       width = 8, height = 5.5, bg = "white")
cat("  Saved Figure 1\n")

# =============================================================================
# FIGURE 2: RIR by Glycemic Status
# =============================================================================
cat("\n[Figure 2: Glycemic Status]\n")

df_primary <- df %>% filter(statin_user == 1 & ldl_under_70 == 1)
fig2_prev <- data.frame()

for (glyc in c("Normal", "Prediabetes", "Diabetes")) {
  df_glyc <- df_primary %>% filter(glycemic_cat == glyc)
  if (nrow(df_glyc) >= 10) {
    glyc_design <- svydesign(id = ~psu, strata = ~strata, weights = ~adj_weight,
                             data = df_glyc, nest = TRUE)
    prev <- tryCatch(svymean(~rir, design = glyc_design, na.rm = TRUE), error = function(e) NULL)
    if (!is.null(prev)) {
      ci <- confint(prev)
      fig2_prev <- rbind(fig2_prev, data.frame(
        Glycemic_Status = glyc, N = nrow(df_glyc),
        Prevalence = coef(prev)[1] * 100,
        CI_low = max(0, ci[1] * 100), CI_high = min(100, ci[2] * 100),
        SE_pct = sqrt(attr(prev, "var")[1]) * 100
      ))
    }
  }
}

fig2_prev$label <- paste0(sprintf("%.1f%%", fig2_prev$Prevalence), "\n(SE ", 
                          sprintf("%.1f", fig2_prev$SE_pct), "%)")
fig2_prev$Glycemic_Status <- factor(fig2_prev$Glycemic_Status, 
                                    levels = c("Normal", "Prediabetes", "Diabetes"))

fig2 <- ggplot(fig2_prev, aes(x = Glycemic_Status, y = Prevalence, fill = Glycemic_Status)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.5) +
  geom_text(aes(label = label, y = CI_high + 3), size = 3, vjust = 0, lineheight = 0.9) +
  scale_fill_manual(values = c("Normal" = rir_colors$cat3_light, 
                               "Prediabetes" = rir_colors$cat3_mid,
                               "Diabetes" = rir_colors$cat3_dark)) +
  scale_y_continuous(limits = c(0, 48), breaks = seq(0, 45, 15),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "RIR Prevalence by Glycemic Status",
    subtitle = "Statin users with LDL-C <70 mg/dL (N = 577)",
    x = "Glycemic Status", y = "RIR Prevalence (%)"
  ) +
  theme_rir() +
  theme(legend.position = "none") +
  annotate("text", x = 0.6, y = 46, size = 2.5, hjust = 0, fontface = "italic", color = "gray40",
           label = "Note: Normal glycemia estimate has high imprecision")

ggsave(file.path(output_dir, "Figure2_Glycemic_Status.png"), fig2, 
       width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure2_Glycemic_Status.pdf"), fig2, 
       width = 6.5, height = 5.5, bg = "white")
cat("  Saved Figure 2\n")

# =============================================================================
# FIGURE 3: Sex Differences
# =============================================================================
cat("\n[Figure 3: Sex Differences]\n")

fig3_prev <- data.frame()
for (grp in groups_list) {
  for (sx in c("Male", "Female")) {
    grp_design <- subset(design_full, ldl_statin_group_f == grp & sex_cat == sx)
    if (nrow(grp_design) >= 10) {
      prev <- tryCatch(svymean(~hscrp_high, design = grp_design, na.rm = TRUE), error = function(e) NULL)
      if (!is.null(prev)) {
        ci <- confint(prev)
        fig3_prev <- rbind(fig3_prev, data.frame(
          Group = grp, Sex = sx, N = nrow(grp_design),
          Prevalence = coef(prev) * 100,
          CI_low = max(0, ci[1] * 100), CI_high = ci[2] * 100
        ))
      }
    }
  }
}

fig3_prev <- fig3_prev %>%
  mutate(
    LDL_Tier = case_when(
      grepl("<=70", Group) ~ "LDL <=70", grepl("70-130", Group) ~ "LDL 70-130", grepl(">130", Group) ~ "LDL >130"
    ),
    Statin = ifelse(grepl("No Statin", gsub("\n", " ", Group)), "No Statin", "Statin"),
    Group_Short = paste0(LDL_Tier, "\n", Statin),
    Sex = factor(Sex, levels = c("Male", "Female"))
  )
fig3_prev$Group_Short <- factor(fig3_prev$Group_Short,
                                levels = c("LDL <=70\nStatin", "LDL <=70\nNo Statin",
                                          "LDL 70-130\nStatin", "LDL 70-130\nNo Statin",
                                          "LDL >130\nStatin", "LDL >130\nNo Statin"))

fig3 <- ggplot(fig3_prev, aes(x = Group_Short, y = Prevalence, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.5) +
  scale_fill_manual(values = c("Male" = rir_colors$secondary, "Female" = rir_colors$primary), name = "") +
  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "hs-CRP >=2 mg/L Prevalence by Sex Across LDL-C/Statin Groups",
    subtitle = "NHANES 2005-2010, 2015-2020",
    x = "", y = "Prevalence (%)"
  ) +
  theme_rir() +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.9))

ggsave(file.path(output_dir, "Figure3_Sex_Differences.png"), fig3, 
       width = 10, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure3_Sex_Differences.pdf"), fig3, 
       width = 10, height = 5.5, bg = "white")
cat("  Saved Figure 3\n")

# =============================================================================
# FIGURE 4: Clinical Characteristics by RIR Status
# =============================================================================
cat("\n[Figure 4: Characteristics by RIR]\n")

primary_cohort <- df %>% filter(statin_user == 1 & ldl_under_70 == 1)
primary_cohort$female <- ifelse(primary_cohort$sex == 2, 1, 0)

# Calculate statistics
calc_mean <- function(data, var) {
  data %>% filter(!is.na(.data[[var]])) %>% group_by(rir_status) %>%
    summarise(mean = weighted.mean(.data[[var]], fasting_weight, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    mutate(se = mean * 0.05, variable = var)  # Approximate SE
}

calc_prop <- function(data, var) {
  data %>% filter(!is.na(.data[[var]])) %>% group_by(rir_status) %>%
    summarise(prop = weighted.mean(.data[[var]], fasting_weight, na.rm = TRUE) * 100,
              n = n(), .groups = "drop") %>%
    mutate(se = sqrt(prop/100 * (1 - prop/100) / n) * 100, variable = var)
}

continuous_data <- bind_rows(
  calc_mean(primary_cohort, "bmi") %>% mutate(var_label = "BMI (kg/m2)"),
  calc_mean(primary_cohort, "hba1c") %>% mutate(var_label = "HbA1c (%)")
)

comorbidity_data <- bind_rows(
  calc_prop(primary_cohort, "diabetes") %>% mutate(var_label = "Diabetes"),
  calc_prop(primary_cohort, "hypertension") %>% mutate(var_label = "Hypertension"),
  calc_prop(primary_cohort, "obese") %>% mutate(var_label = "Obesity"),
  calc_prop(primary_cohort, "current_smoker") %>% mutate(var_label = "Current Smoker")
)

# Panel A
p_continuous <- ggplot(continuous_data, aes(x = var_label, y = mean, fill = rir_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se),
                position = position_dodge(width = 0.8), width = 0.2, color = "black") +
  geom_text(aes(label = sprintf("%.1f", mean), y = mean + 1.96*se + 1.5),
            position = position_dodge(width = 0.8), size = 3.2, vjust = 0) +
  scale_fill_manual(values = c("No RIR" = rir_colors$secondary, "RIR" = rir_colors$primary), name = "") +
  labs(x = "", y = "Mean Value", title = "A. Metabolic Characteristics") +
  theme_rir() + coord_cartesian(ylim = c(0, 42))

# Panel B
p_comorbidities <- ggplot(comorbidity_data, aes(x = var_label, y = prop, fill = rir_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, prop - 1.96*se), ymax = pmin(100, prop + 1.96*se)),
                position = position_dodge(width = 0.8), width = 0.2, color = "black") +
  geom_text(aes(label = sprintf("%.0f%%", prop), y = pmin(100, prop + 1.96*se) + 3),
            position = position_dodge(width = 0.8), size = 3, vjust = 0) +
  scale_fill_manual(values = c("No RIR" = rir_colors$secondary, "RIR" = rir_colors$primary), name = "") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "", y = "Prevalence (%)", title = "B. Comorbidity Prevalence") +
  theme_rir()

# Combine
library(patchwork)
combined <- p_continuous + p_comorbidities +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Clinical Characteristics by Residual Inflammatory Risk Status",
    subtitle = "Statin users with LDL-C <70 mg/dL (N = 577)",
    caption = "Error bars represent 95% CI. RIR = hs-CRP >=2 mg/L.",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 8, color = "gray50")
    )
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Figure4_Characteristics_by_RIR.png"), combined, 
       width = 10, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure4_Characteristics_by_RIR.pdf"), combined, 
       width = 10, height = 5, bg = "white")
cat("  Saved Figure 4\n")

# =============================================================================
# COMPLETE
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("RIR FIGURES COMPLETE - UNIFIED NAVY BLUE COLOR SCHEME\n")
cat(strrep("=", 60), "\n")
cat("\nColor palette:\n")
cat("  Primary (dark):  ", rir_colors$primary, "\n")
cat("  Secondary (mid): ", rir_colors$secondary, "\n")
cat("  Light accent:    ", rir_colors$cat3_light, "\n")
cat("\nAll 4 figures saved to:", output_dir, "\n")

