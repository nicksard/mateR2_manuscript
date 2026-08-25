# ==============================================================================
# Script 02: Run MCMC Chains (Local Parallel Runner for Windows)
# Project: mateR2 Manuscript
# Goal: Load input/mcmc_scenarios_grid.csv and run chains in parallel across
#       6 local CPU cores, saving .rds trace outputs to data/outputs/chains/
# ==============================================================================

library(readr)
library(dplyr)
library(mateR2)
library(future)
library(future.apply)

source(file = "scripts/coancestry3.R")
coancestry(pedigree_df = ,info = T)
?coances
# 1. Setup Output Directories
out_dir <- "data/outputs/chains"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

summary_dir <- "data/outputs"
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

# 2. Load Master Scenario Grid
grid_file <- "input/mcmc_scenarios_grid.csv"
if (!file.exists(grid_file)) {
  stop(sprintf("Error: Grid file '%s' not found.", grid_file))
}

jobs <- read_csv(grid_file, show_col_types = FALSE)
total_jobs <- nrow(jobs)
head(jobs)
cat(sprintf("\n[MCMC Pipeline] Loaded %d total jobs from %s\n", total_jobs, grid_file))

# 3. Configure Windows Parallel Cluster (6 Cores)
num_cores <- 6
plan(multisession, workers = num_cores)

cat(sprintf("[Parallel Setup] Initialized Windows multisession cluster using %d cores.\n\n", num_cores))

# 4. Worker Function to Execute Single Job
run_single_job <- function(i) {
  # Load library inside background worker session
  library(mateR2)
  
  job_row <- jobs[i, ]
  
  # Standard output filename matching local structure
  output_file <- sprintf("%s/job%04d_scen%02d_chain%d.rds", 
                         out_dir, job_row$job_id, job_row$scenario_id, job_row$chain)
  
  # Fast skip check for completed runs
  if (file.exists(output_file)) {
    return(data.frame(
      job_id = job_row$job_id, scenario_id = job_row$scenario_id, 
      chain = job_row$chain, status = "SKIPPED", run_time_sec = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  max_mates_cap <- 10
  
  # Deterministic integer seed generation (prevents POSIXct integer overflow)
  set.seed(12345 + job_row$job_id)
  run_seed <- sample.int(1e6, 1)
  
  tryCatch({
    start_time <- Sys.time()
    
    # Configuration and initial state setup
    config <- create_config_info(
      max_males_per_female = max_mates_cap, 
      max_females_per_male = max_mates_cap
    )
    
    initial_state <- create_initial_counts(
      config_info = config, 
      Np_target   = job_row$target_Np, 
      sr_target   = job_row$target_SR
    )
    
    # Silence C++/R console printing to prevent parallel I/O lockups on Windows
    capture.output({
      mcmc_results <- generate_map_table(
        Np_target            = job_row$target_Np,
        sr_target            = job_row$target_SR,
        mean_mates_target    = job_row$target_MM,
        max_males_per_female = max_mates_cap,
        max_females_per_male = max_mates_cap,
        decay_constant       = job_row$decay_constant,
        np_weight            = job_row$np_weight,
        sr_weight            = job_row$sr_weight,
        mm_weight            = job_row$mm_weight,
        n_iter               = job_row$n_iter,
        burn_in              = job_row$burn_in,
        thin                 = 10,
        initial_method       = initial_state,
        seed                 = run_seed
      )
    })
    
    end_time <- Sys.time()
    run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Attach execution metadata
    mcmc_results$performance_metrics <- data.frame(
      job_id       = job_row$job_id,
      scenario_id  = job_row$scenario_id,
      chain        = job_row$chain,
      run_time_sec = run_time
    )
    
    saveRDS(mcmc_results, output_file)
    
    data.frame(
      job_id = job_row$job_id, scenario_id = job_row$scenario_id, 
      chain = job_row$chain, status = "SUCCESS", run_time_sec = run_time,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    data.frame(
      job_id = job_row$job_id, scenario_id = job_row$scenario_id, 
      chain = job_row$chain, status = sprintf("FAILED: %s", e$message), run_time_sec = 0,
      stringsAsFactors = FALSE
    )
  })
}

# 5. Execute Parallel Runs with Progress Reporting
cat("Starting parallel MCMC sampling...\n")
start_total <- Sys.time()

results_list <- future_lapply(seq_len(total_jobs), run_single_job, future.seed = TRUE)

end_total <- Sys.time()
total_sec <- as.numeric(difftime(end_total, start_total, units = "secs"))

# Reset parallel plan back to sequential
plan(sequential)

# 6. Consolidate and Save Execution Log
summary_df <- bind_rows(results_list)
write_csv(summary_df, "data/outputs/mcmc_execution_log.csv")

cat("\n====================================================\n")
cat(sprintf("Parallel Execution Complete in %.1f minutes!\n", total_sec / 60))
cat(sprintf("Successful: %d | Skipped: %d | Failed: %d\n", 
            sum(summary_df$status == "SUCCESS"),
            sum(summary_df$status == "SKIPPED"),
            sum(startsWith(summary_df$status, "FAILED"))))
cat(sprintf("Execution log written to: data/outputs/mcmc_execution_log.csv\n"))
cat("====================================================\n")
