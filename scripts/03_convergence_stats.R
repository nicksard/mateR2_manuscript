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
#   data/outputs/mcmc_convergence_summary.csv  -- one row per scenario, MEDIAN
#                                                 across seed sets, plus _env envelope
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

# --- One row per scenario: median across seed sets, with the envelope kept -----
# The manuscript reports the MEDIAN across seed sets and states the envelope as
# a bound. A worst-case envelope cannot distinguish a cell that fails to
# converge from one that drew badly once: for MM, R-hat > 1.05 is 8 of 30 by
# envelope but 4 by median, and N_P=1600/SR=2/MM=2 gives 1.363, 1.003, 1.036.
#
# Bare column names carry the MEDIAN, because those are the values the paper
# reports and the values script 04 plots. The envelope is retained under an
# _env suffix so both live in one object and nothing is transcribed by hand.
#
# obs_* is taken from the seed set that produced the median (or worst) error for
# that parameter, never a median of the observations themselves: error_* is
# computed from obs_*, so averaging them independently would describe a cell
# that no seed set actually produced.
worst_obs <- function(obs, err) obs[which.max(err)]
med_obs   <- function(obs, err) obs[order(err)[ceiling(length(err) / 2)]]

summary_tbl <- per_seed %>%
  group_by(scenario_id, profile_name, target_Np, target_SR, target_MM) %>%
  summarise(
    # --- median: what Section 4.1 reports -------------------------------------
    obs_Np       = med_obs(obs_Np, error_Np_pct),
    obs_SR       = med_obs(obs_SR, error_SR_pct),
    obs_MM       = med_obs(obs_MM, error_MM_pct),
    error_Np_pct = median(error_Np_pct),
    error_SR_pct = median(error_SR_pct),
    error_MM_pct = median(error_MM_pct),
    chains_used  = sum(chains_used),          # total chains behind the cell
    rhat_Np = median(rhat_Np), ess_Np = median(ess_Np),
    rhat_SR = median(rhat_SR), ess_SR = median(ess_SR),
    rhat_MM = median(rhat_MM), ess_MM = median(ess_MM),

    # --- envelope: the bound quoted in the caption ----------------------------
    obs_Np_env       = worst_obs(obs_Np, error_Np_pct),
    obs_SR_env       = worst_obs(obs_SR, error_SR_pct),
    obs_MM_env       = worst_obs(obs_MM, error_MM_pct),
    error_Np_pct_env = max(error_Np_pct),
    error_SR_pct_env = max(error_SR_pct),
    error_MM_pct_env = max(error_MM_pct),
    rhat_Np_env = max(rhat_Np), ess_Np_env = min(ess_Np),
    rhat_SR_env = max(rhat_SR), ess_SR_env = min(ess_SR),
    rhat_MM_env = max(rhat_MM), ess_MM_env = min(ess_MM),

    n_seed_sets  = n(),
    .groups = "drop"
  ) %>%
  arrange(scenario_id)

# The original 18 columns first, in the original order, so script 04 is
# untouched -- it now plots medians because the bare names hold them.
schema <- c("scenario_id","profile_name","target_Np","target_SR","target_MM",
            "obs_Np","obs_SR","obs_MM","error_Np_pct","error_SR_pct",
            "error_MM_pct","chains_used","rhat_Np","ess_Np","rhat_SR","ess_SR",
            "rhat_MM","ess_MM")
env_cols <- c("obs_Np_env","obs_SR_env","obs_MM_env",
              "error_Np_pct_env","error_SR_pct_env","error_MM_pct_env",
              "rhat_Np_env","ess_Np_env","rhat_SR_env","ess_SR_env",
              "rhat_MM_env","ess_MM_env","n_seed_sets")
write_csv(summary_tbl[, c(schema, env_cols)],
          "data/outputs/mcmc_convergence_summary.csv")

envelope <- summary_tbl        # name kept for the report block below

# --- Report -------------------------------------------------------------------
cat(sprintf("\n[Complete] %d scenarios x %d seed sets.\n",
            nrow(envelope), max(envelope$n_seed_sets)))
for (p in c("Np","SR","MM")) {
  e  <- envelope[[paste0("error_", p, "_pct")]]
  r  <- envelope[[paste0("rhat_",  p)]]
  s_ <- envelope[[paste0("ess_",   p)]]
  re <- envelope[[paste0("rhat_",  p, "_env")]]
  se <- envelope[[paste0("ess_",   p, "_env")]]
  ee <- envelope[[paste0("error_", p, "_pct_env")]]
  cat(sprintf(
    "  %-3s  median err %6.3f%%  (envelope max %7.3f%%)  |  rhat>1.05: %2d/%d median, %2d/%d envelope  |  ESS<400: %2d/%d median, %2d/%d envelope  |  min ESS %8.0f\n",
    p, median(e), max(ee),
    sum(r > 1.05), nrow(envelope), sum(re > 1.05), nrow(envelope),
    sum(s_ < 400), nrow(envelope), sum(se < 400), nrow(envelope), min(se)))
}
