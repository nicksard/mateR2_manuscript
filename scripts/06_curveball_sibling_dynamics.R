# ==============================================================================
# Master Simulation Pipeline: Sibling Dynamics & Parental Class Distributions
# (Posterior Ensemble Approach using mateR2 Internal Functions)
# ==============================================================================

library(mateR2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(foreach)
library(doParallel)

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Experimental Scenario Grid Setup
# ------------------------------------------------------------------------------
np_vals      <- c(100, 200, 400, 800, 1600)
sr_vals      <- c(2)
mm_vals      <- c(2)
sample_sizes <- c(100, 500, 1000)
n_reps       <- 10

# MCMC & Curveball Parameters
n_iter <- 250000; burn_in <- 25000; thin <- 10
w_np <- 50; w_sr <- 25; w_mm <- 50
curveball_I <- 0.25
n_ensemble_draws <- 5  # Number of thinned posterior draws per replicate
min_fertility <- 1000; max_fertility <- 5000; fertility_type <- "uniform"

# Build Full Factorial Grid
scenarios <- expand.grid(
  Np       = np_vals,
  SR       = sr_vals,
  MM       = mm_vals,
  N_sample = sample_sizes,
  rep      = 1:n_reps
)

cat(sprintf("Total Simulation Scenarios to Execute: %d\n", nrow(scenarios)))

# Create Output Directories
if (!dir.exists("data/outputs")) dir.create("data/outputs", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Parallel Cluster Execution (Posterior Sampling)
# ------------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 1)
cat(sprintf("Registering %d CPU cores...\n", n_cores))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("n_iter", "burn_in", "thin", "w_np", "w_sr", 
                              "w_mm", "curveball_I", "n_ensemble_draws", 
                              "min_fertility", "max_fertility", 
                              "fertility_type"), envir = environment())

start_time <- Sys.time()

df_raw <- foreach(
  i = seq_len(nrow(scenarios)),
  .combine  = rbind,
  .packages = c("mateR2", "dplyr", "tidyr", "tibble")
) %dopar% {
  
  curr  <- scenarios[i, ]
  max_m <- max(10, curr$MM * 2)
  
  # Stage 1: MCMC Execution with Posterior Ensemble Extraction
  mcmc_res <- tryCatch({
    generate_map_table(
      Np_target            = curr$Np,
      sr_target            = curr$SR,
      mean_mates_target    = curr$MM,
      max_males_per_female = max_m,
      max_females_per_male = max_m,
      n_iter               = n_iter,
      burn_in              = burn_in,
      thin                 = thin,
      np_weight            = w_np,
      sr_weight            = w_sr,
      mm_weight            = w_mm,
      sample_ensemble      = TRUE,
      n_ensemble           = n_ensemble_draws,
      max_error_pct        = Inf  # Unbiased posterior sampling across stationary distribution
    )
  }, error = function(e) NULL)
  
  if (!is.null(mcmc_res)) {
    
    # Extract thinned posterior draws (fallback to map_table if empty)
    samples_pool <- mcmc_res$ensemble_c_k %||% mcmc_res$mcmc_output$samples
    base_df      <- mcmc_res$map_table
    
    rep_draw_results <- list()
    
    if (!is.null(samples_pool) && length(samples_pool) > 0) {
      draw_indices <- round(seq(1, length(samples_pool), length.out = min(n_ensemble_draws, length(samples_pool))))
      
      for (d_idx in seq_along(draw_indices)) {
        counts_vec        <- samples_pool[[draw_indices[d_idx]]]
        draw_df           <- base_df
        draw_df$MAP_Count <- counts_vec
        
        mat_bin <- tryCatch(mp_table_to_matrix(mp_table = draw_df), error = function(e) NULL)
        if (is.null(mat_bin)) next
        
        # Stage 2: Matrix Expansion & Curveball Rewiring
        mat_mixed <- randomize_mating_structure(mat_bin, I = curveball_I)
        
        # Stage 3: Fecundity Assignment & Sub-Sampling
        mat_fecund <- brd.mat.fitness(mat_mixed, min.fert = min_fertility, max.fert = max_fertility, type = fertility_type)
        df_sub     <- mat.sub.sample(mat_fecund, num_offspring = curr$N_sample)
        
        # Reconstruct numeric matrix for package functions
        mat_sub <- xtabs(off1 ~ dads + moms, data = df_sub)
        mat_sub <- matrix(as.numeric(mat_sub), nrow = nrow(mat_sub), ncol = ncol(mat_sub))
        
        # Stage 4: Execute Package Functions
        s_out <- sib.stats(mat_sub)
        p_out <- parent.class.stats(mat_sub)
        
        # Flatten 3-row parent class output (maternal, paternal, overall) into 1 wide row
        p_out_wide <- p_out %>%
          tibble::rownames_to_column(var = "parent_type") %>%
          pivot_wider(
            names_from  = parent_type,
            values_from = c(detected, singletons, doubletons, tripletons, quadrupletons),
            names_glue  = "{parent_type}_{.value}"
          )
        
        if (!is.null(s_out) && !is.null(p_out_wide)) {
          rep_draw_results[[d_idx]] <- cbind(s_out, p_out_wide)
        }
      }
    }
    
    # Average across posterior draws for this scenario replicate
    if (length(rep_draw_results) > 0) {
      avg_stats <- bind_rows(rep_draw_results) %>%
        summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
      
      data.frame(
        Np       = curr$Np,
        SR       = curr$SR,
        MM       = curr$MM,
        N_sample = curr$N_sample,
        rep      = curr$rep,
        avg_stats
      )
    }
  }
}

stopCluster(cl)

end_time <- Sys.time()
cat(sprintf("\nExecution completed in %0.2f minutes.\n", as.numeric(difftime(end_time, start_time, units = "mins"))))

saveRDS(df_raw, "data/outputs/extended_sib_and_parent_stats.rds")

# ==============================================================================
# FIGURE 05: Sibling Group Sizes (Facets on the Right, No Titles)
# ==============================================================================
#df_raw <- readRDS("data/outputs/extended_sib_and_parent_stats.rds")
df_sib_summary <- df_raw %>%
  filter(MM == 2, SR == 2) %>%
  mutate(mean_FS_size = (FS_pairs * 2) / N_sample) %>%
  pivot_longer(
    cols      = c(mean_FS_size, mean_MHS, mean_PHS),
    names_to  = "sib_type",
    values_to = "mean_siblings"
  ) %>%
  mutate(
    Sib_Label = factor(
      sib_type,
      levels = c("mean_FS_size", "mean_MHS", "mean_PHS"),
      labels = c("Full-Siblings / Offspring", "Maternal Half-Sibs / Offspring", "Paternal Half-Sibs / Offspring")
    ),
    Sample_Color = factor(N_sample, levels = c(100, 500, 1000), labels = c("N = 100", "N = 500", "N = 1000"))
  ) %>%
  group_by(Np, Sib_Label, Sample_Color) %>%
  summarise(
    mean_val = mean(mean_siblings, na.rm = TRUE),
    sd_val   = sd(mean_siblings, na.rm = TRUE),
    se_val   = sd(mean_siblings, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

df_sib_summary %>% as.data.frame()

# Write to file
write_csv(df_sib_summary, "figures/Figure_05_Sibling_Group_Sizes_Summary.csv")


# Plot: facet_grid(Sib_Label ~ .) places stripSib_Label# Plot: facet_grid(Sib_Label ~ .) places strips vertically on the RIGHT side
p_sibs <- ggplot(
  df_sib_summary, 
  aes(x = Np, y = mean_val, color = Sample_Color, group = Sample_Color)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
    width = 30, linewidth = 0.5
  ) +
  scale_x_continuous(breaks = np_vals) +
  scale_color_manual(
    name   = "Sample Size (N)", 
    values = c("N = 100" = "#E69F00", "N = 500" = "#56B4E9", "N = 1000" = "#009E73")
  ) +
  facet_grid(Sib_Label ~ ., scales = "free_y") +  # Right-side strip placement
  labs(
    x = expression(bold("Total Parent Pool Scale (") * bolditalic(N[P]) * bold(")")),
    y = "Mean Sibling Count per Sampled Juvenile"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.6),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold")
  )

print(p_sibs)
ggsave("figures/Figure_05_Sibling_Group_Sizes.png", plot = p_sibs, width = 8.5, height = 7.5, dpi = 300)

# ==============================================================================
# FIGURE 06: Parent Yield Classes (No Titles, Right/Top Standard Grid)
# ==============================================================================

df_yield_grid <- df_raw %>%
  filter(MM == 2, SR == 2) %>%
  pivot_longer(
    cols      = c(maternal_singletons, maternal_doubletons, paternal_singletons, paternal_doubletons),
    names_to  = "class_raw",
    values_to = "count"
  ) %>%
  mutate(
    Parent_Sex = ifelse(grepl("^maternal_", class_raw), "Maternal Parents", "Paternal Parents"),
    Yield_Tier = factor(
      ifelse(grepl("singletons$", class_raw), "singletons", "doubletons"),
      levels = c("singletons", "doubletons"),
      labels = c("Singletons (1 Offspring)", "Doubletons (2 Offspring)")
    ),
    Sample_Color = factor(N_sample, levels = c(100, 500, 1000), labels = c("N = 100", "N = 500", "N = 1000"))
  ) %>%
  group_by(Np, Parent_Sex, Yield_Tier, Sample_Color) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    sd_count   = sd(count, na.rm = TRUE),
    se_count   = sd(count, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

p_yield_grid <- ggplot(
  df_yield_grid,
  aes(x = Np, y = mean_count, color = Sample_Color, group = Sample_Color)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = mean_count - se_count, ymax = mean_count + se_count),
    width = 30, linewidth = 0.5
  ) +
  scale_x_continuous(breaks = np_vals) +
  scale_color_manual(
    name   = "Sample Size (N)", 
    values = c("N = 100" = "#E69F00", "N = 500" = "#56B4E9", "N = 1000" = "#009E73")
  ) +
  facet_grid(Parent_Sex ~ Yield_Tier, scales = "free_y") +
  labs(
    x = expression(bold("Total Parent Pool Scale (") * bolditalic(N[P]) * bold(")")),
    y = "Mean Detected Parents"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.6),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold")
  )

print(p_yield_grid)
ggsave("figures/Figure_06_Parent_Yield_Classes.png", plot = p_yield_grid, width = 9.0, height = 7.0, dpi = 300)

# Write to file
write_csv(df_yield_grid, "figures/Figure_06_Parent_Yield_Grid_Summary.csv")
df_yield_grid %>% as.data.frame()
