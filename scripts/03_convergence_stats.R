# ==============================================================================
# Script 03: MCMC Convergence & Fidelity Diagnostics (Local Loop)
# Project: mateR2 Manuscript
# Goal: Load 4 chains per scenario from data/outputs/chains/, discard burn-in,
#       and calculate per-parameter convergence (R-hat, ESS) and fidelity error.
# ==============================================================================

library(readr)
library(dplyr)
library(coda)

# 1. Setup Directories & Load Grid
# ------------------------------------------------------------------------------
chains_dir  <- "data/outputs/chains"
summary_dir <- "data/outputs"
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

grid_file <- "input/mcmc_scenarios_grid.csv"
if (!file.exists(grid_file)) {
  stop(sprintf("Error: Grid file '%s' not found.", grid_file))
}

jobs <- read_csv(grid_file, show_col_types = FALSE)

# Group jobs by scenario to pair the 4 independent chains per scenario
scenarios <- jobs %>%
  group_by(scenario_id, target_Np, target_SR, target_MM, decay_constant, burn_in, n_iter) %>%
  summarize(
    job_ids = list(job_id),
    chains  = list(chain),
    .groups = "drop"
  )

num_scenarios <- nrow(scenarios)
cat(sprintf("\n[Diagnostics] Starting sequential analysis across %d scenarios...\n\n", num_scenarios))

# 2. Sequential Processing Loop
# ------------------------------------------------------------------------------
results_list <- vector("list", num_scenarios)
i <- 1
for (i in seq_len(num_scenarios)) {
  scen <- scenarios[i, ]
  j_ids <- unlist(scen$job_ids)
  c_ids <- unlist(scen$chains)
  
  trace_Np <- list()
  trace_SR <- list()
  trace_MM <- list()
  
  obs_Np_vec <- c()
  obs_SR_vec <- c()
  obs_MM_vec <- c()
  
  valid_chains_loaded <- 0
  c_idx <- 1
  # Loop through all 4 chains for this scenario
  for (c_idx in seq_along(c_ids)) {
    # Match output filename format from Script 02
    file_path <- sprintf("%s/job%04d_scen%02d_chain%d.rds", 
                         chains_dir, j_ids[c_idx], scen$scenario_id, c_ids[c_idx])

    if (!file.exists(file_path)) next
    
    res <- readRDS(file_path)

    chain_history <- res$mcmc_output$history

    if (is.null(chain_history) || nrow(chain_history) == 0) next
    
    valid_chains_loaded <- valid_chains_loaded + 1

    
    # Slicing burn-in: account for C++ internal thin step (default 10)
    thin_step_used <- ifelse(!is.null(res$mcmc_output$thin_step), res$mcmc_output$thin_step, 10)
    burn_in_rows   <- floor(scen$burn_in)
    
    if (burn_in_rows < nrow(chain_history)) {
      post_burn <- chain_history[(burn_in_rows + 1):nrow(chain_history), ]
    } else {
      post_burn <- chain_history
    }
    
    # Store un-thinned trace vectors into coda mcmc objects
    trace_Np[[valid_chains_loaded]] <- mcmc(post_burn$Np, start = min(post_burn$iteration), thin = 1)
    trace_SR[[valid_chains_loaded]] <- mcmc(post_burn$sr, start = min(post_burn$iteration), thin = 1)
    trace_MM[[valid_chains_loaded]] <- mcmc(post_burn$mean_mates, start = min(post_burn$iteration), thin = 1)
    
    # Capture MAP or mean realized values across post-burn-in sample
    obs_Np_vec <- c(obs_Np_vec, mean(post_burn$Np, na.rm = TRUE))
    obs_SR_vec <- c(obs_SR_vec, mean(post_burn$sr, na.rm = TRUE))
    obs_MM_vec <- c(obs_MM_vec, mean(post_burn$mean_mates, na.rm = TRUE))
  }
  
  # Skip scenario if no valid chain files were loaded
  if (valid_chains_loaded == 0) {
    cat(sprintf("  -> Scenario %d: No completed chain files found. Skipping.\n", scen$scenario_id))
    next
  }
  
  # Convert lists to coda mcmc.list objects
  mcmc_list_Np <- mcmc.list(trace_Np)
  mcmc_list_SR <- mcmc.list(trace_SR)
  mcmc_list_MM <- mcmc.list(trace_MM)
  
  # Compute Gelman-Rubin R-hat safely
  rhat_Np <- NA; rhat_SR <- NA; rhat_MM <- NA
  if (valid_chains_loaded > 1) {
    rhat_Np <- tryCatch(gelman.diag(mcmc_list_Np, autoburnin = FALSE, multivariate = FALSE)$psrf[1, 1], error = function(e) NA_real_)
    rhat_SR <- tryCatch(gelman.diag(mcmc_list_SR, autoburnin = FALSE, multivariate = FALSE)$psrf[1, 1], error = function(e) NA_real_)
    rhat_MM <- tryCatch(gelman.diag(mcmc_list_MM, autoburnin = FALSE, multivariate = FALSE)$psrf[1, 1], error = function(e) NA_real_)
  }
  
  # Compute Effective Sample Size (ESS) safely across all chains
  ess_Np <- tryCatch(sum(effectiveSize(mcmc_list_Np)), error = function(e) NA_real_)
  ess_SR <- tryCatch(sum(effectiveSize(mcmc_list_SR)), error = function(e) NA_real_)
  ess_MM <- tryCatch(sum(effectiveSize(mcmc_list_MM)), error = function(e) NA_real_)
  
  # Calculate mean realized targets across chains
  obs_Np <- mean(obs_Np_vec, na.rm = TRUE)
  obs_SR <- mean(obs_SR_vec, na.rm = TRUE)
  obs_MM <- mean(obs_MM_vec, na.rm = TRUE)
  
  # Construct diagnostic summary row
  results_list[[i]] <- data.frame(
    scenario_id  = scen$scenario_id,
    profile_name = "Equal_Scaled",
    target_Np    = scen$target_Np,
    target_SR    = scen$target_SR,
    target_MM    = scen$target_MM,
    
    # Realized MAP / Mean Demographic Values
    obs_Np       = obs_Np,
    obs_SR       = obs_SR,
    obs_MM       = obs_MM,
    
    # Percentage Fidelity Errors
    error_Np_pct = abs(obs_Np - scen$target_Np) / scen$target_Np * 100,
    error_SR_pct = abs(obs_SR - scen$target_SR) / scen$target_SR * 100,
    error_MM_pct = abs(obs_MM - scen$target_MM) / scen$target_MM * 100,
    
    # Explicit Convergence Diagnostics
    chains_used  = valid_chains_loaded,
    rhat_Np      = rhat_Np,  ess_Np = ess_Np,
    rhat_SR      = rhat_SR,  ess_SR = ess_SR,
    rhat_MM      = rhat_MM,  ess_MM = ess_MM,
    
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("Processed Scenario %02d/%02d (Chains: %d) | Rhat_MM: %.2f | ESS_MM: %.0f\n", 
              i, num_scenarios, valid_chains_loaded, rhat_MM, ess_MM))
}

# 3. Consolidate and Save Final Summary
# ------------------------------------------------------------------------------
final_diagnostics <- bind_rows(results_list)

if (nrow(final_diagnostics) > 0) {
  output_csv <- file.path(summary_dir, "mcmc_convergence_summary.csv")
  write_csv(final_diagnostics, output_csv)
  
  cat("\n====================================================\n")
  cat(sprintf("Diagnostics Complete! Summary saved to: %s\n", output_csv))
  cat("====================================================\n")
  print(head(final_diagnostics, 10))
} else {
  cat("\n[Error] No scenario outputs were processed. Please check chain RDS paths.\n")
}
