# ==============================================================================
# Script 05: Bipartite Network Topology & Micro/Macro Validation
# Evaluates Connected Components (Macro) and Co-Mating Pair Density (Micro)
# Across Geometric Parent Pool Scales (Np = 100, 200, 400, 800, 1600)
# ==============================================================================

library(mateR2)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

set.seed(42)
# Run from the repository root (open the .Rproj, or `Rscript scripts/05_...R`).
# No setwd(): a hard-coded absolute path makes the script unrunnable for anyone
# who is not the author.

# ------------------------------------------------------------------------------
# 1. Helper Function: Calculate Co-Mating Pair Density & Connected Components
# ------------------------------------------------------------------------------
calc_network_topology_stats <- function(mat_binary) {
  # Delegates to mateR2. The previous in-script implementation named the row
  # projection "females" and the column projection "males"; mateR2 matrices are
  # rows = males, columns = females. Because the reported statistic is the mean
  # of the two and numerator/denominator were transposed consistently, the
  # published Figure 4B was unaffected -- but the package functions keep the
  # labels straight and also return per-sex densities.
  tibble::tibble(
    Connected_Components  = mateR2::count_network_components(mat_binary),
    CoMating_Pair_Density = mateR2::comating_pair_density(mat_binary)$mean_density
  )
}

# 2. Configuration & Parameter Grid Linking
# ------------------------------------------------------------------------------
mixing_grid      <- rev(c(1, 0.75, 0.50, 0.25, 0.10, 0.05, 0.02, 0.01, 0))
n_swap_reps        <- 5     # Rewiring replicates per matrix state
# Posterior draws are now produced by script 02, pooled across chains, and
# Section 3.2 reports 10 per scenario. Read them rather than re-deriving.
data_dir          <- "data/outputs/chains"

# Read job/scenario metadata
job_grid <- read_csv("input/mcmc_scenarios_grid.csv", show_col_types = FALSE)

# Filter for the Main Text Anchor Scenarios: SR = 2.0 and MM = 2.0
target_jobs <- job_grid %>%
  filter(target_SR == 2.0, target_MM == 2.0)

cat(sprintf("Selected %d jobs matching Anchor Scenario (SR=2.0, MM=2.0) across Np scales.\n", nrow(target_jobs)))

# Locate saved job files. Script 02 writes one file per scenario x seed set,
# named scenario_<id>_seed_<set>.rds.
rds_files <- list.files(data_dir, pattern = "^scenario_\\d+_seed_\\d+\\.rds$",
                        full.names = TRUE)
if (!length(rds_files)) {
  stop("No chain outputs in ", data_dir, ". Run scripts/02_run_mcmc_chains.R first.")
}
# One row per scenario: the topology sweep does not vary by seed set.
target_jobs <- target_jobs[!duplicated(target_jobs$scenario_id), ]

results_list <- list()
list_idx     <- 1

# ------------------------------------------------------------------------------
# 3. Main Topology Sweep Across Targeted Jobs
# ------------------------------------------------------------------------------
cat("\n======================================================\n")
cat("Starting Topology Sweep Across Target Chains\n")
cat("======================================================\n")

for (i in 1:nrow(target_jobs)) {
  
  curr_job <- target_jobs[i, ]
  
  pattern_str <- sprintf("^scenario_%02d_seed_", curr_job$scenario_id)
  f_path      <- rds_files[grep(pattern_str, basename(rds_files))]
  
  if (length(f_path) == 0) {
    warning(sprintf("No output for scenario_id %d. Skipping.", curr_job$scenario_id))
    next
  }
  
  f_path <- f_path[1]   # seed sets are equivalent here; take the first
  cat(sprintf("Processing scenario %d (Np = %d) -> %s\n",
              curr_job$scenario_id, curr_job$target_Np, basename(f_path)))
  
  mcmc_data <- readRDS(f_path)
  
  # Ensure MAP table format
  map_df <- mcmc_data$map_table
  if (!"MAP_Count" %in% colnames(map_df) && "Count" %in% colnames(map_df)) {
    map_df$MAP_Count <- map_df$Count
  }
  
  mat_map <- mp_table_to_matrix(mp_table = map_df)
  
  # Posterior draws come straight from script 02: 10 states pooled across the
  # four chains. The old code took 5 draws from ONE chain's sample list, giving
  # 20 per scenario where Section 3.2 reports 10.
  samples_pool <- mcmc_data$posterior_draws
  
  ensemble_matrices <- list()
  if (length(samples_pool)) {
    for (d_i in seq_along(samples_pool)) {
      draw_df <- map_df
      draw_df$MAP_Count <- samples_pool[[d_i]]
      ensemble_matrices[[d_i]] <- tryCatch(mp_table_to_matrix(mp_table = draw_df),
                                           error = function(e) NULL)
    }
    ensemble_matrices <- Filter(Negate(is.null), ensemble_matrices)
  }
  if (!length(ensemble_matrices)) ensemble_matrices <- list(mat_map)
  
  # --- A. MAP Sweep ---
  for (I_val in mixing_grid) {
    for (rep in 1:n_swap_reps) {
      mat_shuffled <- randomize_mating_structure(
        binary_IIM = mat_map,
        I          = I_val
      )
      
      graph_stats <- calc_network_topology_stats(mat_shuffled)
      
      results_list[[list_idx]] <- tibble(
        Job_ID                = curr_job$job_id,
        Scenario_ID           = curr_job$scenario_id,
        Target_Np             = curr_job$target_Np,
        StateType             = "MAP",
        Intensity             = I_val,
        Replicate             = rep,
        Connected_Components  = graph_stats$Connected_Components,
        CoMating_Pair_Density = graph_stats$CoMating_Pair_Density
      )
      list_idx <- list_idx + 1
    }
  }
  
  # --- B. Ensemble Draws Sweep ---
  for (ens_i in seq_along(ensemble_matrices)) {
    mat_ens <- ensemble_matrices[[ens_i]]
    for (I_val in mixing_grid) {
      mat_shuffled <- randomize_mating_structure(
        binary_IIM = mat_ens,
        I          = I_val
      )
      
      graph_stats <- calc_network_topology_stats(mat_shuffled)
      
      results_list[[list_idx]] <- tibble(
        Job_ID                = curr_job$job_id,
        Scenario_ID           = curr_job$scenario_id,
        Target_Np             = curr_job$target_Np,
        StateType             = sprintf("Posterior_Draw_%d", ens_i),
        Intensity             = I_val,
        Replicate             = 1,
        Connected_Components  = graph_stats$Connected_Components,
        CoMating_Pair_Density = graph_stats$CoMating_Pair_Density
      )
      list_idx <- list_idx + 1
    }
  }
}

network_results_df <- bind_rows(results_list) %>%
  mutate(Np = factor(Target_Np, levels = c(100, 200, 400, 800, 1600)))
head(network_results_df)

# ------------------------------------------------------------------------------
# 4. Summary Statistics & Visualizations (Posterior Ensemble Only)
# ------------------------------------------------------------------------------
summary_topology_df <- network_results_df %>%
  filter(StateType != "MAP") %>%
  group_by(Np, Intensity) %>%
  summarize(
    Mean_Components = mean(Connected_Components),
    SE_Components   = sd(Connected_Components) / sqrt(n()),
    Mean_Density    = mean(CoMating_Pair_Density),
    SE_Density      = sd(CoMating_Pair_Density) / sqrt(n()),
    .groups = "drop"
  )

summary_topology_df

# --- Plot A: Connected Components Count (Macro-Topology) ---
p_components <- ggplot(summary_topology_df, aes(x = Intensity, y = Mean_Components, color = Np, group = Np)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "A) Connected Components Count (Macro-Topology)",
    x     = NULL,
    y     = "Mean Component Count",
    color = expression("Parent Pool (" * N[P] * ")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.title      = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

# --- Plot B: Co-Mating Pair Density (Micro-Topology) ---
p_density <- ggplot(summary_topology_df, aes(x = Intensity, y = Mean_Density, color = Np, group = Np)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "B) Co-Mating Pair Density (Micro-Topology)",
    x     = "Mixing Intensity (I)",
    y     = "Mean One-Mode Projection Density",
    color = expression("Parent Pool (" * N[P] * ")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

# --- Combined Patchwork ---
combined_macro_micro_plot <- p_components / p_density +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_macro_micro_plot)

# Export assets
dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_04_Macro_Micro_Plot.png", combined_macro_micro_plot, width = 11, height = 8.5, dpi = 300)
write_csv(summary_topology_df, "figures/topology_grid_summary.csv")
