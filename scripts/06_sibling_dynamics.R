# ==============================================================================
# Script 06: Sibling Dynamics & Parental Class Distributions
# Posterior ensemble -> edge-swap rewiring -> fecundity -> sub-sampling
#
# CHANGED from the original (see manuscript revision notes):
#   * Stage 1 is no longer re-run here. The original called generate_map_table()
#     per replicate with weights 50/25/50, unscaled, and the generate_map_table
#     default decay_constant of -0.5 rather than -0.05, and with
#     max_error_pct = Inf which disables the ensemble quality filter entirely.
#     The result was that Sections 4.3/4.4 were NOT run at SR = 2.0: ensemble
#     draws sat at a realized SR near 1.4-1.6, and the half-sibling asymmetry
#     reported in 4.3 is a function of that realized SR, not of the target.
#     This script now reads the posterior draws produced by script 02, so it
#     uses the same chains, the same baseline weight profile and the same decay
#     constant as Sections 3.1/3.2 -- by construction rather than by assertion.
#   * Adds readr and tibble (write_csv and rownames_to_column were called
#     without them) and drops the use of `%||%`, which is base R only from 4.4.0.
#   * Re-seeds between Stage 1 and Stage 2. The C++ sampler now draws from R's
#     RNG, so downstream sampling would otherwise inherit the chain's RNG path.
#   * "curveball" naming removed: the routine is a pairwise checkerboard edge
#     swap (Gotelli & Entsminger 2003; Miklos & Podani 2004), not Curveball.
# ==============================================================================

library(mateR2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(parallel)

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Experimental Scenario Grid
# ------------------------------------------------------------------------------
np_vals      <- c(100, 200, 400, 800, 1600)
sr_target    <- 2
mm_target    <- 2
sample_sizes <- c(100, 500, 1000)
n_reps       <- 10

mixing_I      <- 0.25          # Section 4.2: topology saturates by I = 0.25
min_fertility <- 1000; max_fertility <- 5000; fertility_type <- "uniform"
base_seed     <- 20260101

chain_dir <- "data/outputs/chains"
if (!dir.exists(chain_dir)) {
  stop("No chain outputs in ", chain_dir, ". Run scripts/02_run_mcmc_chains.R first.")
}

job_grid <- read_csv("input/mcmc_scenarios_grid.csv", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. Load the posterior ensemble for each Np at the anchor scenario
# ------------------------------------------------------------------------------
load_pool <- function(np) {
  sc <- job_grid %>%
    filter(target_Np == np, target_SR == sr_target, target_MM == mm_target) %>%
    pull(scenario_id) %>% unique()
  if (!length(sc)) stop("No scenario for Np = ", np, " at the anchor SR/MM.")
  fs <- list.files(chain_dir, full.names = TRUE,
                   pattern = sprintf("^scenario_%02d_seed_\\d+\\.rds$", sc[1]))
  if (!length(fs)) stop("No chain files for scenario ", sc[1], " (Np = ", np, ").")
  draws <- list(); base_df <- NULL
  for (f in fs) {
    x <- readRDS(f)
    if (is.null(base_df)) base_df <- x$map_table
    draws <- c(draws, x$posterior_draws)
  }
  list(base_df = base_df, draws = draws)
}

pools <- setNames(lapply(np_vals, load_pool), as.character(np_vals))
cat(sprintf("Loaded posterior pools: %s\n",
            paste(sprintf("Np=%d (%d draws)", np_vals,
                          vapply(pools, function(p) length(p$draws), integer(1))),
                  collapse = ", ")))

scenarios <- expand.grid(Np = np_vals, N_sample = sample_sizes, rep = seq_len(n_reps))
scenarios$SR <- sr_target; scenarios$MM <- mm_target
scenarios$task_id <- seq_len(nrow(scenarios))
cat(sprintf("Total simulation tasks: %d\n", nrow(scenarios)))

if (!dir.exists("data/outputs")) dir.create("data/outputs", recursive = TRUE)
if (!dir.exists("figures"))      dir.create("figures", recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. One task: draw -> matrix -> rewire -> fecundity -> sub-sample -> stats
# ------------------------------------------------------------------------------
run_task <- function(k) {
  curr <- scenarios[k, ]
  pool <- pools[[as.character(curr$Np)]]
  if (!length(pool$draws)) return(NULL)

  # Deterministic per-task seed. Stage 1 is already done, so this seeds only the
  # rewiring, fecundity and sub-sampling steps.
  set.seed(base_seed + curr$task_id)

  # Cycle deterministically through the pooled posterior draws.
  d_idx   <- ((curr$rep - 1L) %% length(pool$draws)) + 1L
  draw_df <- pool$base_df
  draw_df$MAP_Count <- pool$draws[[d_idx]]

  mat_bin <- tryCatch(mateR2::mp_table_to_matrix(mp_table = draw_df),
                      error = function(e) NULL)
  if (is.null(mat_bin)) return(NULL)

  mat_mixed  <- mateR2::randomize_mating_structure(mat_bin, I = mixing_I)
  mat_fecund <- mateR2::brd.mat.fitness(mat_mixed, min.fert = min_fertility,
                                        max.fert = max_fertility,
                                        type = fertility_type)
  df_sub <- mateR2::mat.sub.sample(mat_fecund, num_offspring = curr$N_sample)

  # Realized demography of the draw actually used -- so the SR the figures were
  # produced under is recorded rather than assumed.
  realized_Np <- nrow(mat_bin) + ncol(mat_bin)
  realized_SR <- nrow(mat_bin) / ncol(mat_bin)

  mat_sub <- xtabs(off1 ~ dads + moms, data = df_sub)
  mat_sub <- matrix(as.numeric(mat_sub), nrow = nrow(mat_sub), ncol = ncol(mat_sub))

  s_out <- mateR2::sib.stats(mat_sub)
  p_out <- mateR2::parent.class.stats(mat_sub)

  p_wide <- p_out %>%
    tibble::rownames_to_column(var = "parent_type") %>%
    tidyr::pivot_wider(
      names_from  = parent_type,
      values_from = c(detected, singletons, doubletons, tripletons, quadrupletons),
      names_glue  = "{parent_type}_{.value}"
    )
  if (is.null(s_out) || is.null(p_wide)) return(NULL)

  cbind(data.frame(Np = curr$Np, SR = curr$SR, MM = curr$MM,
                   N_sample = curr$N_sample, rep = curr$rep,
                   realized_Np = realized_Np, realized_SR = realized_SR,
                   draw_index = d_idx),
        s_out, p_wide)
}

# ------------------------------------------------------------------------------
# 4. Execute
# ------------------------------------------------------------------------------
n_cores <- as.integer(Sys.getenv("MATER2_WORKERS",
                                 unset = max(1, parallel::detectCores() - 1)))
cat(sprintf("Running %d tasks on %d workers...\n", nrow(scenarios), n_cores))
start_time <- Sys.time()

if (n_cores > 1) {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, { library(mateR2); library(dplyr); library(tidyr); library(tibble) })
  clusterExport(cl, c("scenarios", "pools", "mixing_I", "min_fertility",
                      "max_fertility", "fertility_type", "base_seed"),
                envir = environment())
  res <- parLapplyLB(cl, seq_len(nrow(scenarios)), run_task)
} else {
  res <- lapply(seq_len(nrow(scenarios)), run_task)
}

df_raw <- bind_rows(Filter(Negate(is.null), res))
cat(sprintf("\nCompleted %d/%d tasks in %0.2f minutes.\n",
            nrow(df_raw), nrow(scenarios),
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat(sprintf("Realized demography of draws used: mean Np = %.1f, mean SR = %.3f\n",
            mean(df_raw$realized_Np), mean(df_raw$realized_SR)))

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
