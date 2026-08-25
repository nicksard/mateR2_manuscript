# ==============================================================================
# Script 04: Generate mateR2 MCMC Diagnostics & Supplementary Fidelity
# Project: mateR2 Manuscript
# Output:
#   1. figures/Figure_03_MCMC_Diagnostics.png (8.5 x 11 in, no title/legends, tags a/b)
#   2. figures/Figure_S01_Target_Fidelity_Error.png (Supplementary Error Grid)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(patchwork)

figures_dir <- "figures"
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

summary_csv <- "data/outputs/mcmc_convergence_summary.csv"
if (!file.exists(summary_csv)) {
  stop(sprintf("Error: Summary file '%s' not found.", summary_csv))
}

# ------------------------------------------------------------------------------
# 1. Void Condition & Data Processing
# ------------------------------------------------------------------------------
is_void <- function(sr, mm) { mm < ((sr + 1) / 2) }

final_agg <- read_csv(summary_csv, show_col_types = FALSE) %>%
  rowwise() %>%
  mutate(is_structurally_void = is_void(target_SR, target_MM)) %>%
  ungroup() %>%
  mutate(
    target_Np_factor = factor(target_Np, levels = c(100, 200, 400, 800, 1600)),
    target_MM_factor = factor(target_MM, levels = c(1, 2, 4)),
    target_SR_factor = factor(target_SR, levels = c(1.0, 2.0, 4.0), labels = c("SR: 1.0", "SR: 2.0", "SR: 4.0")),
    
    # Text labels
    label_rhat_Np = ifelse(is_structurally_void | is.na(rhat_Np), "N/A", sprintf("%.2f", rhat_Np)),
    label_rhat_SR = ifelse(is_structurally_void | is.na(rhat_SR), "N/A", sprintf("%.2f", rhat_SR)),
    label_rhat_MM = ifelse(is_structurally_void | is.na(rhat_MM), "N/A", sprintf("%.2f", rhat_MM)),
    
    label_ess_Np = ifelse(is_structurally_void | is.na(ess_Np), "N/A", sprintf("%.0f", ess_Np)),
    label_ess_SR = ifelse(is_structurally_void | is.na(ess_SR), "N/A", sprintf("%.0f", ess_SR)),
    label_ess_MM = ifelse(is_structurally_void | is.na(ess_MM), "N/A", sprintf("%.0f", ess_MM)),
    
    label_err_Np = ifelse(is_structurally_void, "N/A", sprintf("%.1f%%", error_Np_pct)),
    label_err_SR = ifelse(is_structurally_void, "N/A", sprintf("%.1f%%", error_SR_pct)),
    label_err_MM = ifelse(is_structurally_void, "N/A", sprintf("%.1f%%", error_MM_pct)),
    
    # Categorical bins
    cat_rhat_Np = cut(rhat_Np, breaks = c(-Inf, 1.05, 1.10, Inf), labels = c("Optimal", "Warning", "Poor")),
    cat_rhat_SR = cut(rhat_SR, breaks = c(-Inf, 1.05, 1.10, Inf), labels = c("Optimal", "Warning", "Poor")),
    cat_rhat_MM = cut(rhat_MM, breaks = c(-Inf, 1.05, 1.10, Inf), labels = c("Optimal", "Warning", "Poor")),
    
    cat_ess_Np = cut(ess_Np, breaks = c(-Inf, 200, 1000, 10000, Inf), labels = c("Low", "Good", "High", "Very High")),
    cat_ess_SR = cut(ess_SR, breaks = c(-Inf, 200, 1000, 10000, Inf), labels = c("Low", "Good", "High", "Very High")),
    cat_ess_MM = cut(ess_MM, breaks = c(-Inf, 200, 1000, 10000, Inf), labels = c("Low", "Good", "High", "Very High"))
  )

# ------------------------------------------------------------------------------
# 2. Minimalist Theme (No Legends, Clean Lines)
# ------------------------------------------------------------------------------
diag_theme <- function(show_y = TRUE, show_x = TRUE) {
  t <- theme_minimal(base_size = 11) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
      strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5),
      strip.text       = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 6)),
      axis.title       = element_text(face = "bold", size = 11),
      axis.text        = element_text(color = "black", size = 9.5),
      legend.position  = "none" # Strictly no legends (described in caption)
    )
  if (!show_y) {
    t <- t + theme(
      axis.text.y  = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  if (!show_x) {
    t <- t + theme(
      axis.title.x = element_blank()
    )
  }
  return(t)
}

# ------------------------------------------------------------------------------
# 3. Plot Engines (Grayscale / Print-Friendly Palettes)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ESS Engine: Aligned so "Low" (Warning) is Dark Gray, "Very High" is White
# ------------------------------------------------------------------------------
gen_ess_plot <- function(data_subset, var_cat, var_label, title, show_y = TRUE, show_x = TRUE) {
  ggplot(data_subset, aes(x = target_Np_factor, y = target_MM_factor)) +
    geom_tile(aes(fill = {{var_cat}}), color = "black", linewidth = 0.35) +
    geom_tile(data = data_subset %>% filter(is_structurally_void), 
              fill = "gray95", color = "black", linewidth = 0.5) +
    # Dynamic text color: White on dark 'Low' cells, Black on all other light cells
    geom_text(aes(label = {{var_label}}, 
                  color = ifelse({{var_cat}} == "Low", "white", "black")), 
              size = 2.5, fontface = "bold") +
    scale_color_identity() +
    scale_fill_manual(
      values = c(
        "Very High" = "#ffffff", # Excellent (Clean White)
        "High"      = "#f0f0f0", # Strong (Very Light Gray)
        "Good"      = "#cccccc", # Acceptable (Light-Medium Gray)
        "Low"       = "#737373"  # Warning Flag (Dark Gray)
      ),
      drop = FALSE, na.value = "gray95"
    ) +
    diag_theme(show_y, show_x) +
    facet_grid(target_SR_factor ~ .) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = title,
      x = expression(bold("Target Parents (") * bolditalic(N[P]) * bold(")")),
      y = expression(bold("Target Mean Mates (") * bolditalic(MM) * bold(")"))
    )
}

# ------------------------------------------------------------------------------
# R-hat Engine: Aligned so "Poor" (Warning) is Dark Gray, "Optimal" is White
# ------------------------------------------------------------------------------
gen_rhat_plot <- function(data_subset, var_cat, var_label, title, show_y = TRUE, show_x = TRUE) {
  ggplot(data_subset, aes(x = target_Np_factor, y = target_MM_factor)) +
    geom_tile(aes(fill = {{var_cat}}), color = "black", linewidth = 0.35) +
    geom_tile(data = data_subset %>% filter(is_structurally_void), 
              fill = "gray95", color = "black", linewidth = 0.5) +
    # Dynamic text color: White on dark 'Poor' cells, Black on other cells
    geom_text(aes(label = {{var_label}}, 
                  color = ifelse({{var_cat}} == "Poor", "white", "black")), 
              size = 2.5, fontface = "bold") +
    scale_color_identity() +
    scale_fill_manual(
      values = c(
        "Optimal" = "#ffffff", # Fully converged (Clean White)
        "Warning" = "#cccccc", # Intermediate warning (Light-Medium Gray)
        "Poor"    = "#737373"  # Non-converged warning (Dark Gray)
      ),
      drop = FALSE, na.value = "gray95"
    ) +
    diag_theme(show_y, show_x) +
    facet_grid(target_SR_factor ~ .) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = title,
      x = expression(bold("Target Parents (") * bolditalic(N[P]) * bold(")")),
      y = expression(bold("Target Mean Mates (") * bolditalic(MM) * bold(")"))
    )
}

# ------------------------------------------------------------------------------
# 4. Build Figure 03 (Full-Page 8.5 x 11 in, Tagged a & b in Top-Left Corners)
# ------------------------------------------------------------------------------
# Top Row: ESS Subplots (No "a)" in titles)
p_ess_np <- gen_ess_plot(final_agg, cat_ess_Np, label_ess_Np, expression(bold("Parents (") * bolditalic(N[P]) * bold(")")), show_y = TRUE, show_x = FALSE)
p_ess_sr <- gen_ess_plot(final_agg, cat_ess_SR, label_ess_SR, expression(bold("Sex Ratio (") * bolditalic(SR) * bold(")")), show_y = FALSE, show_x = FALSE)
p_ess_mm <- gen_ess_plot(final_agg, cat_ess_MM, label_ess_MM, expression(bold("Mean Mates (") * bolditalic(MM) * bold(")")), show_y = FALSE, show_x = FALSE)

tier_a <- (p_ess_np | p_ess_sr | p_ess_mm)

# Bottom Row: R-hat Subplots (No "b)" in titles)
p_rhat_np <- gen_rhat_plot(final_agg, cat_rhat_Np, label_rhat_Np, expression(bold("Parents (") * bolditalic(N[P]) * bold(")")), show_y = TRUE, show_x = TRUE)
p_rhat_sr <- gen_rhat_plot(final_agg, cat_rhat_SR, label_rhat_SR, expression(bold("Sex Ratio (") * bolditalic(SR) * bold(")")), show_y = FALSE, show_x = TRUE)
p_rhat_mm <- gen_rhat_plot(final_agg, cat_rhat_MM, label_rhat_MM, expression(bold("Mean Mates (") * bolditalic(MM) * bold(")")), show_y = FALSE, show_x = TRUE)

tier_b <- (p_rhat_np | p_rhat_sr | p_rhat_mm)

# Stack vertically with independent 'a)' and 'b)' tags in the top-left outer corners
fig03_combined <- (tier_a / tier_b) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = "",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 16, margin = margin(r = 5, b = 5)),
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
    )
  )
fig03_combined

file_fig03 <- file.path(figures_dir, "Figure_03_MCMC_Diagnostics.png")
ggsave(file_fig03, fig03_combined, width = 8.5, height = 11.0, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Build Figure S01 (Supplementary Target Fidelity Percent Error Grid)
# ------------------------------------------------------------------------------
max_err_val <- max(c(final_agg$error_Np_pct, final_agg$error_SR_pct, final_agg$error_MM_pct), na.rm = TRUE)
err_limits  <- c(0, ceiling(max_err_val))

gen_fidelity_plot <- function(data_subset, var_fill, var_label, title, show_y = TRUE) {
  ggplot(data_subset, aes(x = target_Np_factor, y = target_MM_factor)) +
    geom_tile(aes(fill = {{var_fill}}), color = "black", linewidth = 0.4) +
    geom_tile(data = data_subset %>% filter(is_structurally_void), fill = "gray95", color = "black", linewidth = 0.8) +
    geom_text(aes(label = {{var_label}}), size = 3.6, fontface = "bold", color = "black") +
    scale_fill_gradient(low = "#ffffff", high = "#969696", na.value = "gray95", limits = err_limits) +
    diag_theme(show_y = show_y, show_x = TRUE) +
    facet_grid(target_SR_factor ~ .) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = title,
      x = expression(bold("Target Parents (") * bolditalic(N[P]) * bold(")")),
      y = expression(bold("Target Mean Mates (") * bolditalic(MM) * bold(")"))
    )
}

p_fid_np <- gen_fidelity_plot(final_agg, error_Np_pct, label_err_Np, expression(bold("Parents (") * bolditalic(N[P]) * bold(")")), show_y = TRUE)
p_fid_sr <- gen_fidelity_plot(final_agg, error_SR_pct, label_err_SR, expression(bold("Sex Ratio (") * bolditalic(SR) * bold(")")), show_y = FALSE)
p_fid_mm <- gen_fidelity_plot(final_agg, error_MM_pct, label_err_MM, expression(bold("Mean Mates (") * bolditalic(MM) * bold(")")), show_y = FALSE)

fig_s01_fidelity <- (p_fid_np | p_fid_sr | p_fid_mm) +
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = "",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.margin = margin(10, 10, 10, 10)
    )
  )
fig_s01_fidelity
file_figs01 <- file.path(figures_dir, "Figure_S01_Target_Fidelity_Error.png")
ggsave(file_figs01, fig_s01_fidelity, width = 14, height = 7.0, dpi = 300)

cat(sprintf("\n[SUCCESS] Production figures generated:\n -> %s\n -> %s\n", file_fig03, file_figs01))
