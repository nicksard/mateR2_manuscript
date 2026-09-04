#!/usr/bin/env Rscript

# ==============================================================================
# Provenance audit
# Project: mateR2 Manuscript
#
# Confirms that every number and figure the manuscript cites traces to ONE
# coherent execution of the current code. Checks, in order:
#
#   1. Working tree is clean (uncommitted edits mean outputs may not match code)
#   2. Every output postdates the newest source file
#   3. Outputs postdate the installed mateR2 build
#   4. The mateR2 version recorded inside the chain files matches the installed one
#   5. All 90 chain files are present and share one mateR2 version
#   6. Derived outputs postdate the chains they are derived from
#
# Not part of the pipeline. Run before submission:
#     Rscript scripts/audit_provenance.R
# Exits non-zero if anything is stale.
# ==============================================================================

fail <- character()
note <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", ifelse(ok, "ok", "!!"), msg))
  if (!ok) fail <<- c(fail, msg)
}
hdr <- function(x) cat(sprintf("\n%s\n%s\n", x, strrep("-", nchar(x))))

cat(strrep("=", 78), "\nProvenance audit\n", strrep("=", 78), "\n", sep = "")

# --- 1. Clean tree ------------------------------------------------------------
hdr("1. Repository state")
st <- suppressWarnings(try(system2("git", c("status", "--porcelain"),
                                   stdout = TRUE, stderr = TRUE), silent = TRUE))
if (inherits(st, "try-error")) {
  cat("  [--] git unavailable; skipping\n")
} else {
  dirty <- st[!grepl("^\\?\\?", st)]           # ignore untracked
  note(length(dirty) == 0,
       sprintf("working tree clean (%d tracked file(s) modified)", length(dirty)))
  if (length(dirty)) cat(paste0("       ", dirty, collapse = "\n"), "\n")
  sha <- try(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), silent = TRUE)
  if (!inherits(sha, "try-error")) cat("       HEAD:", sha, "\n")
}

# --- 2. Newest source vs outputs ---------------------------------------------
hdr("2. Outputs vs source code")
src <- list.files("scripts", "\\.R$", full.names = TRUE)
src_t <- max(file.mtime(src))
cat(sprintf("  newest source file: %s (%s)\n",
            basename(src[which.max(file.mtime(src))]), format(src_t, "%Y-%m-%d %H:%M")))

outs <- c(list.files("data/outputs", full.names = TRUE, recursive = TRUE),
          list.files("figures", full.names = TRUE))
outs <- outs[!grepl("weight_sensitivity", outs)]   # step 10 is opt-in and cached
if (!length(outs)) {
  note(FALSE, "no outputs found -- run scripts/00_run_all.R")
} else {
  stale <- outs[file.mtime(outs) < src_t]
  note(length(stale) == 0,
       sprintf("all %d outputs postdate the newest source file (%d stale)",
               length(outs), length(stale)))
  if (length(stale)) cat(paste0("       ", head(stale, 12), collapse = "\n"), "\n")
}

# --- 3 & 4. Package build and recorded version --------------------------------
hdr("3. mateR2 build")
pkg_v <- as.character(utils::packageVersion("mateR2"))
pkg_t <- file.mtime(system.file("DESCRIPTION", package = "mateR2"))
cat(sprintf("  installed mateR2 %s, built %s\n", pkg_v, format(pkg_t, "%Y-%m-%d %H:%M")))
if (length(outs)) {
  older <- outs[file.mtime(outs) < pkg_t]
  note(length(older) == 0,
       sprintf("all outputs postdate the installed package build (%d predate it)",
               length(older)))
}

hdr("4. Version recorded inside the chain files")
cf <- list.files("data/outputs/chains", "^scenario_.*\\.rds$", full.names = TRUE)
note(length(cf) == 90, sprintf("chain files present: %d (expected 90)", length(cf)))
if (length(cf)) {
  vs <- unique(vapply(cf, function(f) {
    x <- readRDS(f); v <- x$mateR2_version
    if (is.null(v)) NA_character_ else as.character(v)
  }, character(1)))
  note(length(vs) == 1, sprintf("single mateR2 version across chains: %s",
                                paste(vs, collapse = ", ")))
  note(identical(vs[1], pkg_v),
       sprintf("chains built with the installed version (%s vs %s)", vs[1], pkg_v))
}

# --- 5. Derived outputs postdate the chains -----------------------------------
hdr("5. Derivation order")
if (length(cf)) {
  ch_t <- max(file.mtime(cf))
  derived <- c("data/outputs/mcmc_convergence_summary.csv",
               "data/outputs/mcmc_convergence_by_seed.csv",
               "data/outputs/extended_sib_and_parent_stats.rds",
               "data/outputs/sibling_asymmetry_cells.csv",
               "figures/Figure_03_MCMC_Diagnostics.png",
               "figures/Figure_S01_Target_Fidelity_Error.png")
  derived <- derived[file.exists(derived)]
  bad <- derived[file.mtime(derived) < ch_t]
  note(length(bad) == 0,
       sprintf("all %d derived outputs postdate the newest chain file (%d stale)",
               length(derived), length(bad)))
  if (length(bad)) cat(paste0("       ", bad, collapse = "\n"), "\n")
  missing <- setdiff(c("data/outputs/mcmc_convergence_summary.csv",
                       "data/outputs/sibling_asymmetry_cells.csv"), derived)
  if (length(missing)) note(FALSE, paste("missing expected output:", paste(missing, collapse = ", ")))
}

# --- Verdict ------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
if (!length(fail)) {
  cat("PASS - every output traces to the current code and package build.\n")
} else {
  cat(sprintf("FAIL - %d issue(s):\n", length(fail)))
  cat(paste0("  - ", fail, collapse = "\n"), "\n")
  cat("\nRe-run the full pipeline:  Rscript scripts/00_run_all.R\n")
  cat("(delete data/outputs/chains first to force regeneration)\n")
}
cat(strrep("=", 78), "\n")
# Non-zero exit so this can gate a release step or CI check. Under
# `source()` in RStudio quit() would close the session, so only exit when
# running non-interactively.
if (!interactive()) quit(status = if (length(fail)) 1 else 0, save = "no")
