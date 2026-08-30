#!/usr/bin/env Rscript

# ==============================================================================
# Script 01: Setup Validation Grid (Unified Linear-Scaled Weights)
# Project: mateR2 Manuscript
# Goal: Generate biologically viable validation scenarios and assign 
#       scaled Goldilocks MCMC parameters across geometric population sizes.
# Output: input/mcmc_cluster_jobs_grid.csv
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)

# Ensure input directory exists
if (!dir.exists("input")) dir.create("input", recursive = TRUE)

# 1. Base Demographics (Geometric Expansion up to 1,600 parents)
total_parents        <- c(100, 200, 400, 800, 1600)
sex_ratios           <- c(1.0, 2.0, 4.0)
mean_mates           <- c(1.0, 2.0, 4.0)
max_biological_mates <- 10

# Independent seed sets. run_mcmc_chains() runs and disperses its own chains,
# so the grid is crossed with seed sets rather than with chain indices. The
# worst case across seed sets is what Figure 3 reports: single-seed R-hat
# excursions do not replicate and should not be shown as convergence failures.
seed_sets    <- c(100, 200, 300)

# `thin` governs the retained sample pool (MAP search + posterior ensemble);
# `trace_thin` governs only the diagnostic traces. Unrelated choices -- coupling
# them would move every reported ESS without the chain having changed.
# thin = 100 is ~37x faster than thin = 10 with bit-identical chains.
thin_samples <- 100
thin_traces  <- 10

# 2. Build and Filter the Factorial Grid
base_demographics <- expand.grid(
  target_Np = total_parents,
  target_SR = sex_ratios,
  target_MM = mean_mates
) %>% 
  mutate(
    # Bipartite Algebra: Calculate exact average mates required per sex
    required_female_mates = target_MM * (2 * target_SR) / (target_SR + 1),
    required_male_mates   = target_MM * 2 / (target_SR + 1)
  )

# THE BIOLOGICAL FILTER: Require valid floor (>= 1.0) and ceiling (<= 10)
base_demographics_filtered <- base_demographics %>% 
  filter(
    required_female_mates >= 1.0,
    required_male_mates   >= 1.0,
    required_female_mates <= max_biological_mates,
    required_male_mates   <= max_biological_mates
  ) %>% 
  arrange(target_Np, target_SR, target_MM) %>% 
  mutate(scenario_id = row_number()) %>% 
  select(scenario_id, target_Np, target_SR, target_MM)

cat(sprintf("--- Filtered down to %d Biologically Viable Scenarios ---\n", nrow(base_demographics_filtered)))

# 3. Unified MCMC Profile -- LINEAR weight scaling with N_P.
# D_theta is a log-ratio and is scale-free; what changes with N_P is the
# granularity of reachable states (one block move shifts N_P by ~3, so
# dD ~ 3/N_P). For w * dD to keep consistent grip across population sizes, w
# must scale linearly with N_P. Measured at N_P=1600/SR=1/MM=4: sqrt weights
# give 6.70% N_P error and ESS 17; linear gives 0.069% and ESS 19,521.
unified_profile <- base_demographics_filtered %>% 
  mutate(
    profile_name   = "Equal_Linear",
    length_name    = "1M_Iter",
    decay_constant = -0.05,
    max_mates_cap  = max_biological_mates,
    scale_factor   = target_Np / 100,          # LINEAR (was sqrt(target_Np/100))
    np_weight      = 50.0 * scale_factor,
    sr_weight      = 50.0 * scale_factor,
    mm_weight      = 50.0 * scale_factor,
    n_iter         = 1000000,
    burn_in        = 100000,
    thin           = thin_samples,
    trace_thin     = thin_traces,
    n_chains       = 4
  ) %>% 
  select(-scale_factor) %>% 
  mutate(grid_id = row_number())

# 4. Cross with independent seed sets (each runs its own n_chains chains)
seeds <- data.frame(seed_set = seed_sets)

job_grid <- unified_profile %>% 
  cross_join(seeds) %>% 
  # Smallest N_P first: run_mcmc_chains() holds n_chains full histories live, so
  # the largest scenarios are the memory-hungry ones. Ordering this way lets a
  # parallel backend clear the cheap rows before memory pressure peaks.
  arrange(target_Np, target_SR, target_MM, seed_set) %>% 
  mutate(job_id = row_number()) %>% 
  select(job_id, grid_id, scenario_id, profile_name, length_name, 
         target_Np, target_SR, target_MM, decay_constant, max_mates_cap,
         np_weight, sr_weight, mm_weight, n_iter, burn_in,
         thin, trace_thin, n_chains, seed_set)

# 5. Export Master Job Grid CSV
write_csv(job_grid, "input/mcmc_scenarios_grid.csv")
cat(sprintf("[Success] Exported %d jobs (%d scenarios x %d seed sets) to input/mcmc_scenarios_grid.csv\n", 
            nrow(job_grid), nrow(unified_profile), length(seed_sets)))
