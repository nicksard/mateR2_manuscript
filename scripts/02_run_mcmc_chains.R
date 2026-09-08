#!/usr/bin/env Rscript

# ==============================================================================
# Script 02: Run MCMC Chains
# Project: mateR2 Manuscript
#
# Goal: For every (scenario x seed set) job in input/mcmc_scenarios_grid.csv,
#       run mateR2::run_mcmc_chains(), which executes n_chains dispersed chains
#       and returns coda traces plus per-parameter R-hat / ESS diagnostics.
#       One .rds per job is written to data/outputs/chains/.
#
# CHANGED from the original:
#   * Removed the source() of scripts/coancestry3.R -- that file is not in the
#     repository, so the script died on line 14 of a clean clone.
#   * The unit of work is the SCENARIO, not the chain. run_mcmc_chains() runs
#     and disperses its own chains, which implements Section 3.1's "dispersed
#     network states". The previous per-chain runner called
#     create_initial_counts(), which is deterministic, so all four chains
#     started from the *identical* state and R-hat was not measuring what
#     Section 5.1 claims it measures.
#   * Ten posterior draws are pooled across the n_chains chains and persisted.
#     Section 3.2 reports 10 draws per scenario; the previous script 05 took 5
#     draws PER CHAIN (20 per scenario) from a single chain's sample list.
#     Pooling also mixes better than drawing from one chain.
#   * Full per-iteration histories are NOT persisted. Four 1e6-row histories per
#     job is what made data/outputs/chains 720 MB and unarchivable. Only the
#     thinned diagnostic traces, the MAP table and the diagnostics survive.
#   * Runs from the repository root. No setwd().
# ==============================================================================

library(readr)
library(dplyr)
library(mateR2)

if (utils::packageVersion("mateR2") < "0.1.0") {
  stop("mateR2 >= 0.1.0 required (the C++ sampler must draw from R's RNG). ",
       "Install with pak::pak('nicksard/mateR2').")
}

# --- Parallel workers ---------------------------------------------------------
# run_mcmc_chains() holds n_chains full histories live before cutting traces,
# so each worker needs roughly 250 MB at N_P = 1600. Cap workers accordingly;
# override with MATER2_WORKERS.
n_workers <- as.integer(Sys.getenv("MATER2_WORKERS", unset = "4"))

out_dir <- "data/outputs/chains"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

grid_file <- "input/mcmc_scenarios_grid.csv"
if (!file.exists(grid_file)) {
  stop(sprintf("Grid file '%s' not found. Run scripts/01_setup_scenarios.R first.",
               grid_file))
}
jobs <- read_csv(grid_file, show_col_types = FALSE)
cat(sprintf("Loaded %d jobs across %d scenarios.\n",
            nrow(jobs), length(unique(jobs$scenario_id))))

# --- One job ------------------------------------------------------------------
run_job <- function(job_row) {
  out_file <- file.path("data/outputs/chains",
                        sprintf("scenario_%02d_seed_%d.rds",
                                job_row$scenario_id, job_row$seed_set))
  # Restartable: skip completed jobs.
  if (file.exists(out_file)) return(sprintf("skip  %s", basename(out_file)))

  t0 <- Sys.time()
  res <- tryCatch(
    mateR2::run_mcmc_chains(
      Np_target            = job_row$target_Np,
      sr_target            = job_row$target_SR,
      mean_mates_target    = job_row$target_MM,
      max_males_per_female = job_row$max_mates_cap,
      max_females_per_male = job_row$max_mates_cap,
      n_chains             = job_row$n_chains,
      seed                 = job_row$seed_set,
      n_iter               = job_row$n_iter,
      burn_in              = job_row$burn_in,
      thin                 = job_row$thin,
      trace_thin           = job_row$trace_thin,
      decay_constant       = job_row$decay_constant,
      np_weight            = job_row$np_weight,
      sr_weight            = job_row$sr_weight,
      mm_weight            = job_row$mm_weight,
      verbose              = FALSE,
      # Silence the C++ progress bar and acceptance-rate line. A bar written
      # from a parallel worker is noise in a log and can stall console I/O on
      # Windows -- which is what the original runner's capture.output() wrapper
      # was working around.
      show_progress        = FALSE
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(sprintf("FAIL  scenario %02d seed %d: %s",
                   job_row$scenario_id, job_row$seed_set, conditionMessage(res)))
  }

  # Pool the retained samples across chains and take n_posterior_draws evenly
  # spaced states. This is the ensemble Section 3.2 describes, and it is what
  # scripts 05 and 06 consume.
  n_posterior_draws <- 10
  pooled <- unlist(lapply(res$chains, function(ch) ch$mcmc_output$samples),
                   recursive = FALSE)
  draw_idx <- if (length(pooled) <= n_posterior_draws) seq_along(pooled) else
    unique(round(seq(1, length(pooled), length.out = n_posterior_draws)))

  # Persist diagnostics, traces, MAP tables and posterior draws -- NOT histories.
  slim <- list(
    job          = as.list(job_row),
    diagnostics  = res$diagnostics,
    traces       = res$traces,
    map_tables   = lapply(res$chains, `[[`, "map_table"),
    map_stats    = lapply(res$chains, `[[`, "map_stats"),
    # Best MAP across chains, for scripts that want a single point state.
    map_table    = res$chains[[which.max(vapply(res$chains,
                     function(ch) ch$map_stats$MaxLogProb, numeric(1)))]]$map_table,
    posterior_draws = pooled[draw_idx],
    n_pooled_samples = length(pooled),
    run_time_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    mateR2_version = as.character(utils::packageVersion("mateR2"))
  )
  saveRDS(slim, out_file)
  sprintf("done  scenario %02d seed %d (%.1fs)",
          job_row$scenario_id, job_row$seed_set, slim$run_time_sec)
}

# --- Execute ------------------------------------------------------------------
job_list <- split(jobs, seq_len(nrow(jobs)))

if (n_workers > 1 && requireNamespace("parallel", quietly = TRUE)) {
  cat(sprintf("Running on %d workers.\n", n_workers))
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(mateR2))
  msgs <- parallel::parLapplyLB(cl, job_list, run_job)
} else {
  msgs <- lapply(job_list, run_job)
}

invisible(lapply(msgs, function(m) cat(m, "\n")))
fails <- grep("^FAIL", unlist(msgs), value = TRUE)
cat(sprintf("\n[Complete] %d job files in %s. %d failures.\n",
            length(list.files(out_dir, "\\.rds$")), out_dir, length(fails)))
if (length(fails)) { cat(fails, sep = "\n"); quit(status = 1) }
