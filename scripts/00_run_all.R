#!/usr/bin/env Rscript

# ==============================================================================
# Script 00: Run the full analysis pipeline
# Project: mateR2 Manuscript
#
# Reproduces every number and figure in the manuscript from a clean clone.
#
#   Rscript scripts/00_run_all.R              # everything
#   Rscript scripts/00_run_all.R 02 03 04     # only these steps
#
# Environment variables:
#   MATER2_WORKERS   parallel workers (default 4). run_mcmc_chains() holds
#                    n_chains full histories live, so each worker needs roughly
#                    250 MB at N_P = 1600. Lower this if memory is tight.
#
# Step 02 is restartable: it skips any job whose output file already exists.
# Delete data/outputs/chains to force a full re-run.
#
# Requirements: R >= 4.1 and mateR2 >= 0.1.0. No system software is needed --
# every figure is produced by R packages that renv can restore.
# ==============================================================================

steps <- list(
  list(id = "01", file = "01_setup_scenarios.R",
       what = "Build the scenario x seed-set job grid"),
  list(id = "02", file = "02_run_mcmc_chains.R",
       what = "Run MCMC chains (the long one; restartable)"),
  list(id = "03", file = "03_convergence_stats.R",
       what = "Aggregate R-hat / ESS / error; write summary CSVs"),
  list(id = "04", file = "04_generate_mateR2_figures.R",
       what = "Figures 3 and S1 (MCMC diagnostics, target fidelity)"),
  list(id = "05", file = "05_topology_validation.R",
       what = "Figure 4 (edge-swap topology sweep)"),
  list(id = "06", file = "06_sibling_dynamics.R",
       what = "Figures 5 and 6 (sibling groups, parent yield classes)"),
  list(id = "07", file = "07_figure_1_Bipartite_Network_visuals.R",
       what = "Figure 1 (conceptual bipartite network)"),
  list(id = "08", file = "08_figure_2_worked_Example.R",
       what = "Figure 2 (worked example)"),
  list(id = "09", file = "09_sibling_asymmetry_mechanism.R",
       what = "Sibling asymmetry mechanism (Sections 5.2 / 5.3)")
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  keep <- vapply(steps, function(s) s$id %in% args, logical(1))
  if (!any(keep)) {
    stop("No matching steps. Available: ", paste(vapply(steps, `[[`, "", "id"),
                                                 collapse = ", "))
  }
  steps <- steps[keep]
}

# --- Preflight ----------------------------------------------------------------
if (getRversion() < "4.1") stop("R >= 4.1 required; found ", getRversion())
if (!requireNamespace("mateR2", quietly = TRUE)) {
  stop("mateR2 is not installed. ",
       "devtools::install_github('nicksard/mateR2')")
}
if (utils::packageVersion("mateR2") < "0.1.0") {
  stop("mateR2 >= 0.1.0 required (the C++ sampler must draw from R's RNG); ",
       "found ", utils::packageVersion("mateR2"))
}
if (!dir.exists("scripts")) {
  stop("Run this from the repository root, e.g. `Rscript scripts/00_run_all.R`.")
}

cat(strrep("=", 78), "\n")
cat("mateR2 manuscript pipeline\n")
cat("  R         ", as.character(getRversion()), "\n")
cat("  mateR2    ", as.character(utils::packageVersion("mateR2")), "\n")
cat("  workers   ", Sys.getenv("MATER2_WORKERS", unset = "4"), "\n")
cat("  steps     ", paste(vapply(steps, `[[`, "", "id"), collapse = " "), "\n")
cat(strrep("=", 78), "\n\n")

# --- Run ----------------------------------------------------------------------
t_all <- Sys.time()
timings <- data.frame()

for (s in steps) {
  path <- file.path("scripts", s$file)
  if (!file.exists(path)) stop("Missing script: ", path)

  cat(sprintf("\n[%s] %s\n     %s\n", s$id, s$file, s$what))
  cat(strrep("-", 78), "\n")

  t0 <- Sys.time()
  status <- system2("Rscript", path)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (!identical(status, 0L)) {
    cat(sprintf("\n[FAIL] %s exited with status %d after %.1f min.\n",
                s$file, status, dt))
    stop("Pipeline halted at step ", s$id, ".")
  }
  cat(sprintf("[ok] %s finished in %.1f min\n", s$id, dt))
  timings <- rbind(timings, data.frame(step = s$id, script = s$file, minutes = dt))
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("Pipeline complete in %.1f minutes.\n\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
print(timings, row.names = FALSE, digits = 3)

outs <- c("data/outputs/mcmc_convergence_summary.csv",
          "data/outputs/mcmc_convergence_by_seed.csv")
figs <- if (dir.exists("figures")) list.files("figures", "\\.png$") else character()
cat("\nKey outputs:\n")
for (f in outs) cat(sprintf("  [%s] %s\n", ifelse(file.exists(f), "x", " "), f))
cat(sprintf("  %d figure(s) in figures/\n", length(figs)))
