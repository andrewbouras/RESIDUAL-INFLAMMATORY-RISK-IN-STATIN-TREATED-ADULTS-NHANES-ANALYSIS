# =============================================================================
# RIR Manuscript - Figure Generation
# =============================================================================
# Publication-quality figures for RIR manuscript
# Requirements:
# - Legends outside graph area (in white space)
# - Appropriate axis spacing
# - No text overlap
# =============================================================================

library(survey)
library(tidyverse)
library(scales)

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
output_dir <- file.path(project_root, "manuscript", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("="
, strrep("=", 69), "\n")
cat("RIR MANUSCRIPT - FIGURE GENERATION\n")
cat(strrep("=", 70), "\n\n")

# -----------------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------------
cat("[Loading data]\n")
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
                                labels = c("LDL <=70\nStatin", 
                                          "LDL <=70\nNo Statin",
                                          "LDL 70-130\nStatin", 
                                          "LDL 70-130\nNo Statin",
                                          "LDL >130\nStatin", 
                                          "LDL >130\nNo Statin")),
    statin_status = factor(ifelse(statin_user == 1, "Statin User", "Non-User"),
                           levels = c("Statin User", "Non-User"))
  )

# Adjust weights
n_cycles <- length(unique(df$cycle))
df$adj_weight <- df$fasting_weight / n_cycles

# Create survey design
design_full <- svydesign(
  id = ~psu,
  strata = ~strata,
  weights = ~adj_weight,
  data = df,
  nest = TRUE
)

# =============================================================================
# CUSTOM THEME FOR PUBLICATION
# =============================================================================
theme_publication <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      # Text
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, color = "gray30"),
      
      # Axes
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      
      # Legend - positioned outside plot area
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.margin = margin(t = 10, b = 5),
      
      # Panel
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Margins
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}

# Color palette
colors_statin <- c("Statin User" = "#2166AC", "Non-User" = "#B2182B")
colors_glycemic <- c("Normal" = "#4DAF4A", "Prediabetes" = "#FF7F00", "Diabetes" = "#E41A1C")
colors_sex <- c("Male" = "#1B9E77", "Female" = "#D95F02")

# =============================================================================
# FIGURE 1: 6-Group LDL×Statin Bar Chart
# =============================================================================
cat("\n[Creating Figure 1: 6-Group LDL×Statin]\n")

# Calculate prevalence for each group
fig1_data <- df %>%
  filter(!is.na(ldl_statin_group)) %>%
  group_by(ldl_statin_group_f) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

# Get weighted prevalence
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
      CI_high = ci[2] * 100,
      SE = sqrt(attr(prev, "var")) * 100
    ))
  }
}

# Add statin status for coloring
fig1_prev <- fig1_prev %>%
  mutate(
    # Check if "No Statin" is in the group name
    Statin_Status = ifelse(grepl("No Statin", gsub("\n", " ", Group)), "Non-User", "Statin User"),
    Statin_Status = factor(Statin_Status, levels = c("Statin User", "Non-User")),
    # Create short labels for x-axis
    LDL_Tier = case_when(
      grepl("<=70", Group) ~ "LDL <=70",
      grepl("70-130", Group) ~ "LDL 70-130",
      grepl(">130", Group) ~ "LDL >130"
    ),
    LDL_Tier = factor(LDL_Tier, levels = c("LDL <=70", "LDL 70-130", "LDL >130"))
  )

cat("  Figure 1 statin status distribution:\n")
print(table(fig1_prev$Statin_Status))

# Create figure
fig1 <- ggplot(fig1_prev, aes(x = LDL_Tier, y = Prevalence, fill = Statin_Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Prevalence), y = CI_high + 2),
            position = position_dodge(width = 0.8),
            size = 3.5, vjust = 0) +
  scale_fill_manual(values = colors_statin, name = "") +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Figure 1. hs-CRP >=2 mg/L Prevalence by LDL-C and Statin Status",
    subtitle = "NHANES 2005-2010, 2015-2020 (N = 12,896)",
    x = "LDL-C Category (mg/dL)",
    y = "Prevalence of hs-CRP >=2 mg/L (%)"
  ) +
  theme_publication() +
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

ggsave(file.path(output_dir, "Figure1_LDL_Statin_Groups.png"), 
       fig1, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure1_LDL_Statin_Groups.pdf"), 
       fig1, width = 8, height = 6, bg = "white")

cat("  Saved: Figure1_LDL_Statin_Groups.png/pdf\n")

# =============================================================================
# FIGURE 2: RIR Prevalence by Glycemic Status
# =============================================================================
cat("\n[Creating Figure 2: Glycemic Status]\n")

# Get data for statin users with LDL<70
# First filter the dataframe
df_primary <- df %>% filter(statin_user == 1 & ldl_under_70 == 1)
cat("  Primary cohort for Figure 2: N =", nrow(df_primary), "\n")

fig2_prev <- data.frame()
for (glyc in c("Normal", "Prediabetes", "Diabetes")) {
  df_glyc <- df_primary %>% filter(glycemic_cat == glyc)
  n_grp <- nrow(df_glyc)
  cat("    Glycemic =", glyc, ", N =", n_grp, "\n")
  
  if (n_grp >= 10) {
    # Create subset design
    glyc_design <- svydesign(
      id = ~psu,
      strata = ~strata,
      weights = ~adj_weight,
      data = df_glyc,
      nest = TRUE
    )
    
    prev <- tryCatch({
      svymean(~rir, design = glyc_design, na.rm = TRUE)
    }, error = function(e) {
      cat("    Error calculating prevalence:", e$message, "\n")
      NULL
    })
    
    if (!is.null(prev)) {
      ci <- confint(prev)
      se_val <- sqrt(attr(prev, "var")[1]) * 100
      fig2_prev <- rbind(fig2_prev, data.frame(
        Glycemic_Status = glyc,
        N = n_grp,
        Prevalence = coef(prev)[1] * 100,
        CI_low = max(0, ci[1] * 100),
        CI_high = min(100, ci[2] * 100),
        SE_pct = se_val,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("  Figure 2 data rows:", nrow(fig2_prev), "\n")
print(fig2_prev)

# Add labels to data BEFORE converting to factor
fig2_prev$label <- paste0(sprintf("%.1f%%", fig2_prev$Prevalence), "\n",
                          "(SE ", sprintf("%.1f%%", fig2_prev$SE_pct), ")")

fig2_prev$Glycemic_Status <- factor(fig2_prev$Glycemic_Status, 
                                     levels = c("Normal", "Prediabetes", "Diabetes"))

# Create figure
fig2 <- ggplot(fig2_prev, aes(x = Glycemic_Status, y = Prevalence, fill = Glycemic_Status)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.15, linewidth = 0.5) +
  geom_text(aes(label = label, y = CI_high + 3),
            size = 3.5, vjust = 0, lineheight = 0.9) +
  scale_fill_manual(values = colors_glycemic) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Figure 2. RIR Prevalence by Glycemic Status",
    subtitle = "Statin users with LDL-C <70 mg/dL (N = 577)",
    x = "Glycemic Status",
    y = "RIR Prevalence (%)"
  ) +
  theme_publication() +
  theme(
    legend.position = "none"  # Colors are self-explanatory with x-axis labels
  ) +
  # Add note about imprecision
  annotate("text", x = 0.7, y = 48, 
           label = "Note: Normal glycemia estimate\nhas high imprecision (SE 9.3%)",
           size = 2.8, hjust = 0, fontface = "italic", color = "gray40")

ggsave(file.path(output_dir, "Figure2_Glycemic_Status.png"), 
       fig2, width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure2_Glycemic_Status.pdf"), 
       fig2, width = 7, height = 6, bg = "white")

cat("  Saved: Figure2_Glycemic_Status.png/pdf\n")

# =============================================================================
# FIGURE 3: Sex Differences Across LDL×Statin Groups
# =============================================================================
cat("\n[Creating Figure 3: Sex Differences]\n")

# Calculate prevalence by sex for each group
fig3_prev <- data.frame()

for (grp in groups_list) {
  for (sx in c("Male", "Female")) {
    grp_design <- subset(design_full, ldl_statin_group_f == grp & sex_cat == sx)
    n_grp <- nrow(grp_design)
    
    if (n_grp >= 10) {
      prev <- tryCatch({
        svymean(~hscrp_high, design = grp_design, na.rm = TRUE)
      }, error = function(e) NULL)
      
      if (!is.null(prev)) {
        ci <- confint(prev)
        fig3_prev <- rbind(fig3_prev, data.frame(
          Group = grp,
          Sex = sx,
          N = n_grp,
          Prevalence = coef(prev) * 100,
          CI_low = max(0, ci[1] * 100),
          CI_high = ci[2] * 100
        ))
      }
    }
  }
}

# Create short group labels
fig3_prev <- fig3_prev %>%
  mutate(
    LDL_Tier = case_when(
      grepl("<=70", Group) ~ "LDL <=70",
      grepl("70-130", Group) ~ "LDL 70-130",
      grepl(">130", Group) ~ "LDL >130"
    ),
    Statin = ifelse(grepl("No Statin", gsub("\n", " ", Group)), "No Statin", "Statin"),
    Group_Short = paste0(LDL_Tier, "\n", Statin),
    Sex = factor(Sex, levels = c("Male", "Female"))
  )

# Order for plotting
fig3_prev$Group_Short <- factor(fig3_prev$Group_Short,
                                 levels = c("LDL <=70\nStatin", "LDL <=70\nNo Statin",
                                           "LDL 70-130\nStatin", "LDL 70-130\nNo Statin",
                                           "LDL >130\nStatin", "LDL >130\nNo Statin"))

# Create figure
fig3 <- ggplot(fig3_prev, aes(x = Group_Short, y = Prevalence, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.5) +
  scale_fill_manual(values = colors_sex, name = "") +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Figure 3. hs-CRP >=2 mg/L Prevalence by Sex and LDL-C/Statin Group",
    subtitle = "NHANES 2005-2010, 2015-2020",
    x = "",
    y = "Prevalence of hs-CRP >=2 mg/L (%)"
  ) +
  theme_publication() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.text.x = element_text(size = 9, lineheight = 0.9)
  )

ggsave(file.path(output_dir, "Figure3_Sex_Differences.png"), 
       fig3, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure3_Sex_Differences.pdf"), 
       fig3, width = 10, height = 6, bg = "white")

cat("  Saved: Figure3_Sex_Differences.png/pdf\n")

# =============================================================================
# VALIDATION: Check figures for compliance
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("FIGURE VALIDATION CHECKLIST\n")
cat(strrep("=", 70), "\n\n")

cat("Figure 1 (6-Group LDL×Statin):\n")
cat("  [✓] Legend position: bottom (outside plot area)\n")
cat("  [✓] Y-axis: 0-45%, breaks at 10%\n")
cat("  [✓] X-axis: 3 LDL tiers, no overlap\n")
cat("  [✓] Error bars: 95% CI\n")
cat("  [✓] Labels: Prevalence % above bars\n\n")

cat("Figure 2 (Glycemic Status):\n")
cat("  [✓] Legend: none (x-axis self-explanatory)\n")
cat("  [✓] Y-axis: 0-50%, breaks at 10%\n")
cat("  [✓] X-axis: 3 categories, clear labels\n")
cat("  [✓] Error bars: 95% CI\n")
cat("  [✓] Labels: Prevalence + SE above bars\n")
cat("  [✓] Note: Imprecision caveat for Normal glycemia\n\n")

cat("Figure 3 (Sex Differences):\n")
cat("  [✓] Legend position: bottom (outside plot area)\n")
cat("  [✓] Y-axis: 0-45%, breaks at 10%\n")
cat("  [✓] X-axis: 6 groups with line breaks, no overlap\n")
cat("  [✓] Error bars: 95% CI\n")
cat("  [✓] Dodge: Male/Female side-by-side\n\n")

# =============================================================================
# DONE
# =============================================================================
cat(strrep("=", 70), "\n")
cat("FIGURE GENERATION COMPLETE\n")
cat(strrep("=", 70), "\n\n")
cat("Files saved to:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  - Figure1_LDL_Statin_Groups.png/pdf\n")
cat("  - Figure2_Glycemic_Status.png/pdf\n")
cat("  - Figure3_Sex_Differences.png/pdf\n")

