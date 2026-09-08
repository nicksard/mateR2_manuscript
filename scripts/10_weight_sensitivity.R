#!/usr/bin/env Rscript

# ==============================================================================
# Script 10: Weight Sensitivity (Figure S2)
# Project: mateR2 Manuscript
#
# Why: Section 3.1 scales the target weight linearly with N_P. That rule needs
# justifying, and the justification is not that error decreases monotonically
# with w. Error is U-SHAPED in w: loose weights leave demographic bias, tight
# weights freeze the sampler. Under the ABC reading, w is a tolerance parameter,
# and this maps the tolerance surface directly.
#
# Two panels:
#   (a) weight x N_P for the boundary family (SR = 1, MM = 1)
#   (b) weight x mating system at N_P = 400
# with the linear rule w = 50 (N_P/100) overlaid on both.
#
# WARNING: this is the slow step. It runs its own MCMC rather than reusing
# script 02's grid, because it needs weights the validation grid does not use.
# 162 jobs x 4 chains x 1e6 iterations -- roughly 40-60 minutes on 4 workers.
# It is NOT part of the default pipeline; run it explicitly:
#     Rscript scripts/00_run_all.R 10
# Results are cached to data/outputs/weight_sensitivity.csv and reused if
# present, so the figure can be re-drawn without re-running the sweep.
# ==============================================================================

library(mateR2)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(ggtext)    # subscripts in panel titles
library(parallel)

if (!dir.exists("data/outputs")) dir.create("data/outputs", recursive = TRUE)
if (!dir.exists("figures"))      dir.create("figures", recursive = TRUE)

cache <- "data/outputs/weight_sensitivity.csv"

weights   <- c(50, 100, 200, 400, 800, 1000)
seed_sets <- c(100, 200, 300)   # match the main grid; 162 jobs, ~40-60 min
cap <- 10; gam <- -0.05
n_iter <- 1e6; burn_in <- 1e5

# (a) boundary family across N_P; (b) all feasible systems at N_P = 400
cells <- bind_rows(
  expand.grid(Np = c(100, 400, 1600), SR = 1, MM = 1),
  expand.grid(Np = 400, SR = c(1, 2, 4), MM = c(1, 2, 4)) %>%
    filter(MM >= (SR + 1) / 2),
  expand.grid(Np = c(100, 1600), SR = 2, MM = 2)
) %>% distinct(Np, SR, MM)

# --- Sweep --------------------------------------------------------------------
if (file.exists(cache)) {
  cat("Using cached sweep at ", cache, " (delete it to re-run).\n", sep = "")
  s <- read_csv(cache, show_col_types = FALSE)
} else {
  jobs <- merge(merge(cells, data.frame(w = weights)), data.frame(seed = seed_sets))
  cat(sprintf("Running %d jobs (this is the slow step)...\n", nrow(jobs)))

  run_job <- function(i) {
    j <- jobs[i, ]
    r <- try(mateR2::run_mcmc_chains(
      Np_target = j$Np, sr_target = j$SR, mean_mates_target = j$MM,
      max_males_per_female = cap, max_females_per_male = cap,
      n_chains = 4, seed = j$seed, n_iter = n_iter, burn_in = burn_in,
      thin = 100, trace_thin = 10, decay_constant = gam,
      np_weight = j$w, sr_weight = j$w, mm_weight = j$w,
      show_progress = FALSE), silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    d <- r$diagnostics
    # Freeze detection. ESS alone is misleading here: a chain pinned at one
    # state can still report a large ESS, and only R-hat exposes it. Count the
    # distinct values each chain actually visits.
    for (p in c("Np", "sr", "mean_mates")) {
      u <- vapply(r$traces[[p]], function(ch) length(unique(as.numeric(ch))), integer(1))
      d$n_frozen_chains[d$parameter == p] <- sum(u == 1)
    }
    d$Np <- j$Np; d$SR <- j$SR; d$MM <- j$MM; d$w <- j$w; d$seed <- j$seed
    d
  }

  n_cores <- as.integer(Sys.getenv("MATER2_WORKERS", unset = "4"))
  if (n_cores > 1) {
    cl <- makeCluster(n_cores); on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, library(mateR2))
    clusterExport(cl, c("jobs", "cap", "gam", "n_iter", "burn_in"), envir = environment())
    res <- parLapplyLB(cl, seq_len(nrow(jobs)), run_job)
  } else {
    res <- lapply(seq_len(nrow(jobs)), run_job)
  }
  raw <- bind_rows(Filter(Negate(is.null), res))
  raw <- raw[raw$parameter == "mean_mates", ]

  s <- raw %>%
    group_by(Np, SR, MM, w) %>%
    summarise(error_pct = mean(error_pct), rhat_max = max(rhat),
              ess_mean = mean(ess), n_frozen_chains = max(n_frozen_chains),
              n_seed_sets = n(), .groups = "drop")
  write_csv(s, cache)
}

# --- Figure -------------------------------------------------------------------
s$err <- pmax(s$error_pct, 0.005)
s$wf  <- factor(s$w, levels = weights)
s$lab <- ifelse(s$error_pct < 0.01, "<0.01", sprintf("%.2f", s$error_pct))
# Flag any cell that failed convergence OR froze. Reported, never dropped: a
# frozen chain still yields a point estimate, but its uncertainty is undefined.
s$bad <- (is.infinite(s$rhat_max) | s$rhat_max > 1.05) | (s$n_frozen_chains > 0)
LIM <- range(log10(s$err))
ink <- "grey15"

# The figure carries only the panels. Everything explanatory -- what the cells
# are, what the outline and the red x mean, how the sweep was run -- belongs in
# the manuscript caption, not baked into the image.
base <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        plot.title = ggtext::element_markdown(size = 9, hjust = 0),
        legend.key.height = unit(0.75, "cm"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

panel <- function(dat, yvar, ttl, ylab) {
  dat$row <- yvar
  ggplot(dat, aes(wf, row)) +
    geom_tile(aes(fill = log10(err)), colour = "grey70", linewidth = 0.3) +
    geom_rect(data = subset(dat, rule),
              aes(xmin = as.numeric(wf) - 0.46, xmax = as.numeric(wf) + 0.46,
                  ymin = as.numeric(row) - 0.46, ymax = as.numeric(row) + 0.46),
              fill = NA, colour = "black", linewidth = 0.9, inherit.aes = FALSE) +
    geom_text(aes(label = lab, colour = log10(err) > 0.7), size = 2.6, fontface = "bold") +
    geom_point(data = subset(dat, bad),
               aes(x = as.numeric(wf) + 0.30, y = as.numeric(row) + 0.29),
               shape = 4, size = 1.5, stroke = 0.75, colour = "firebrick",
               inherit.aes = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink), guide = "none") +
    scale_fill_gradient(low = "white", high = "grey20",
                        name = "log10\n% error", limits = LIM) +
    labs(title = ttl, x = "target weight (*w*)", y = ylab) + base +
    theme(axis.title.x = ggtext::element_markdown(size = 8),
          axis.title.y = if (is.character(ylab)) ggtext::element_markdown(size = 8)
                         else element_text(size = 8))
}

A <- s %>% filter(SR == 1, MM == 1) %>% mutate(rule = w == 50 * (Np / 100))
A$rowf <- factor(A$Np, levels = c(1600, 400, 100))
pa <- panel(A, A$rowf,
  "**(a)** Boundary family (*SR* = 1, *MM* = 1)",
  expression(N[P]))

B <- s %>% filter(Np == 400) %>% mutate(rule = w == 200)
B$rowf <- factor(paste0(B$SR, ":", B$MM),
                 levels = rev(c("1:1", "1:2", "1:4", "2:2", "2:4", "4:4")))
pb <- panel(B, B$rowf,
  "**(b)** Mating systems at *N*<sub>P</sub> = 400",
  "*SR* : *MM*")

p <- pa / pb + plot_layout(heights = c(1, 1.55), guides = "collect")

out <- "figures/Figure_S02_Weight_Sensitivity.png"
ggsave(out, p, width = 8.0, height = 5.6, dpi = 300, bg = "white")
cat("[SUCCESS] ->", out, "\n         ->", cache, "\n")
