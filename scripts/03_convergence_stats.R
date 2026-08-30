#!/usr/bin/env Rscript

# ==============================================================================
# Script 03: Convergence and Fidelity Summary
# Project: mateR2 Manuscript
#
# Goal: Aggregate the per-job outputs of script 02 into the two summary tables
#       the manuscript and script 04 consume.
#
# Outputs:
#   data/outputs/mcmc_convergence_by_seed.csv  -- one row per scenario x seed set
#   data/outputs/mcmc_convergence_summary.csv  -- one row per scenario, WORST CASE
#                                                 across seed sets. Original
#                                                 18-column schema, so script 04
#                                                 runs unchanged.
#
# CHANGED from the original:
#   * R-hat and ESS are no longer computed here. mateR2::run_mcmc_chains()
#     computes them, so the diagnostic definition now lives in one place, in a
#     tested package function, rather than being duplicated in a script.
#   * The original loaded 1e6-row histories from disk and used them UN-thinned,
#     while Section 3.1 describes a trace recorded every 10th iteration. The
#     traces are now thinned once, in the package, at `trace_thin`.
#
# WORST CASE is a per-metric envelope: max error, max R-hat, min ESS taken
# independently across seed sets. It is deliberately conservative and is NOT a
# single realisable run -- R-hat and ESS in one cell may come from different
# seed sets. Say so in the figure caption.
# ==============================================================================

library(dplyr)
library(readr)

in_dir <- "data/outputs/chains"
files  <- list.files(in_dir, pattern = "^scenario_\\d+_seed_\\d+\\.rds$",
                     full.names = TRUE)
if (!length(files)) {
  stop("No job outputs in ", in_dir, ". Run scripts/02_run_mcmc_chains.R first.")
}
cat(sprintf("Found %d job files.\n", length(files)))

# --- Long table: one row per scenario x seed set ------------------------------
per_seed <- bind_rows(lapply(files, function(f) {
  x <- readRDS(f)
  d <- x$diagnostics
  g <- function(p, col) d[[col]][d$parameter == p]
  data.frame(
    scenario_id  = x$job$scenario_id,
    seed_set     = x$job$seed_set,
    profile_name = x$job$profile_name,
    target_Np    = x$job$target_Np,
    target_SR    = x$job$target_SR,
    target_MM    = x$job$target_MM,
    obs_Np = g("Np","observed"), obs_SR = g("sr","observed"),
    obs_MM = g("mean_mates","observed"),
    error_Np_pct = g("Np","error_pct"), error_SR_pct = g("sr","error_pct"),
    error_MM_pct = g("mean_mates","error_pct"),
    chains_used  = g("Np","n_chains"),
    rhat_Np = g("Np","rhat"), ess_Np = g("Np","ess"),
    rhat_SR = g("sr","rhat"), ess_SR = g("sr","ess"),
    rhat_MM = g("mean_mates","rhat"), ess_MM = g("mean_mates","ess"),
    run_time_sec = x$run_time_sec,
    stringsAsFactors = FALSE
  )
})) %>% arrange(scenario_id, seed_set)

write_csv(per_seed, "data/outputs/mcmc_convergence_by_seed.csv")

# --- Envelope: one row per scenario, worst case across seed sets ---------------
# obs_* is taken from the seed set that produced the worst error for that
# parameter, so each observed value stays paired with its own error.
worst_obs <- function(obs, err) obs[which.max(err)]

envelope <- per_seed %>%
  group_by(scenario_id, profile_name, target_Np, target_SR, target_MM) %>%
  summarise(
    obs_Np       = worst_obs(obs_Np, error_Np_pct),
    obs_SR       = worst_obs(obs_SR, error_SR_pct),
    obs_MM       = worst_obs(obs_MM, error_MM_pct),
    error_Np_pct = max(error_Np_pct),
    error_SR_pct = max(error_SR_pct),
    error_MM_pct = max(error_MM_pct),
    chains_used  = sum(chains_used),          # total chains behind the cell
    rhat_Np = max(rhat_Np), ess_Np = min(ess_Np),
    rhat_SR = max(rhat_SR), ess_SR = min(ess_SR),
    rhat_MM = max(rhat_MM), ess_MM = min(ess_MM),
    n_seed_sets  = n(),
    .groups = "drop"
  ) %>%
  arrange(scenario_id)

# Original 18-column schema, in the original order, so script 04 is untouched.
schema <- c("scenario_id","profile_name","target_Np","target_SR","target_MM",
            "obs_Np","obs_SR","obs_MM","error_Np_pct","error_SR_pct",
            "error_MM_pct","chains_used","rhat_Np","ess_Np","rhat_SR","ess_SR",
            "rhat_MM","ess_MM")
write_csv(envelope[, schema], "data/outputs/mcmc_convergence_summary.csv")

# --- Report -------------------------------------------------------------------
cat(sprintf("\n[Complete] %d scenarios x %d seed sets.\n",
            nrow(envelope), max(envelope$n_seed_sets)))
for (p in c("Np","SR","MM")) {
  e <- envelope[[paste0("error_", p, "_pct")]]
  r <- envelope[[paste0("rhat_", p)]]
  s <- envelope[[paste0("ess_",  p)]]
  cat(sprintf("  %-3s  median err %6.3f%%  max %7.3f%%  |  rhat>1.05: %2d/%d  |  ESS<400: %2d/%d  min ESS %8.0f\n",
              p, median(e), max(e), sum(r > 1.05), nrow(envelope),
              sum(s < 400), nrow(envelope), min(s)))
}
