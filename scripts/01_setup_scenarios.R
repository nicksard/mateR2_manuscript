#!/usr/bin/env Rscript

# ==============================================================================
# Script 01: Setup Validation Grid (Unified Equal-Scaled Weights)
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

# 3. Unified "Equal-Scaled" MCMC Profile
# Dynamic weight scaling with N_P via sqrt(N_P / 100)
unified_profile <- base_demographics_filtered %>% 
  mutate(
    profile_name   = "Equal_Scaled",
    length_name    = "1M_Iter",
    decay_constant = -0.05,
    scale_factor   = sqrt(target_Np / 100),
    np_weight      = 50.0 * scale_factor,
    sr_weight      = 50.0 * scale_factor,
    mm_weight      = 50.0 * scale_factor,
    n_iter         = 1000000,
    burn_in        = 100000
  ) %>% 
  select(-scale_factor) %>% 
  mutate(grid_id = row_number())

# 4. Cross with 4 Independent MCMC Chains per Scenario
chains <- data.frame(chain = 1:4)

job_grid <- unified_profile %>% 
  cross_join(chains) %>% 
  mutate(job_id = row_number()) %>% 
  select(job_id, grid_id, scenario_id, profile_name, length_name, 
         target_Np, target_SR, target_MM, decay_constant, 
         np_weight, sr_weight, mm_weight, n_iter, burn_in, chain)

# 5. Export Master Job Grid CSV
write_csv(job_grid, "input/mcmc_scenarios_grid.csv")
cat(sprintf("[Success] Exported %d MCMC jobs across %d scenarios to input/mcmc_scenarios_grid.csv\n", 
            nrow(job_grid), nrow(unified_profile)))
