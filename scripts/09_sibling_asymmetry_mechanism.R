#!/usr/bin/env Rscript

# ==============================================================================
# Script 09: Sibling Asymmetry Mechanism (Sections 5.2 / 5.3)
# Project: mateR2 Manuscript
#
# Question: why are maternal and paternal half-sibling group sizes NOT simply
# proportional to the operational sex ratio?
#
# Half-sibling group size is offspring-weighted, so it is a second-moment
# quantity: the offspring-weighted mean parental reproductive success is
# kbar * (1 + CV^2), not kbar. Female reproductive success is a single bounded
# fecundity draw and its CV is near-constant. Male reproductive success is a
# compound sum over a variable number of mates, so its CV is larger and varies
# with the mating system. The two effects partly cancel:
#
#     MHS / PHS  ~  SR * (1 + CV_f^2) / (1 + CV_m^2)
#
# Consequence: the sex with more mates does NOT automatically have the larger
# half-sibling groups, and any near-symmetry is a property of a particular
# (SR, MM) combination rather than of the sex ratio.
#
# Reads the chain outputs written by scripts/02_run_mcmc_chains.R. The five
# (SR, MM) combinations used here are already in the 30-cell validation grid,
# so no additional MCMC is run.
#
# Outputs:
#   data/outputs/sibling_asymmetry_replicates.rds   per-replicate raw values
#   data/outputs/sibling_asymmetry_cells.csv        per-cell summary + tests
#   data/outputs/sibling_asymmetry_by_seedset.csv   between-draw-set variation
# ==============================================================================

library(mateR2)
library(dplyr)
library(readr)
library(parallel)

# --- Design -------------------------------------------------------------------
# Two overlapping sets of cells:
#   (a) MECHANISM GRID -- SR x MM at two population sizes, to test the algebra.
#   (b) ANCHOR COLUMN  -- SR = 2, MM = 2 at every N_P used in Section 4.3, so
#       that section rests on uniform replicate effort across all five scales
#       rather than on whichever cells happened to be run.
mech_np   <- c(400, 800)
sr_vals   <- c(1, 2, 4)
mm_vals   <- c(2, 4)
anchor_np <- c(100, 200, 400, 800, 1600)
anchor_sr <- 2; anchor_mm <- 2
n_sample  <- 1000        # sub-sample depth for the half-sibling statistics
# EQUAL replicate effort in every cell. Must be a multiple of the pooled draw
# count (10 draws x 3 seed sets = 30) so that each posterior draw -- and so each
# seed set -- is used the same number of times. At n = 100 the first ten draws
# are used four times and the rest three, which over-weights one seed set; since
# between-seed-set variation is comparable to the effect being measured, that
# imbalance is not harmless.
n_reps    <- 120
mixing_I  <- 0.25
min_fert  <- 1000; max_fert <- 5000; fert_type <- "uniform"
base_seed <- 20260201

chain_dir <- "data/outputs/chains"
if (!dir.exists(chain_dir)) {
  stop("No chain outputs in ", chain_dir, ". Run scripts/02_run_mcmc_chains.R first.")
}
job_grid <- read_csv("input/mcmc_scenarios_grid.csv", show_col_types = FALSE)

# Cells, dropping any excluded by the graph-theoretic floor MM >= (SR + 1)/2.
mech_cells <- expand.grid(Np = mech_np, SR = sr_vals, MM = mm_vals) %>%
  filter(MM >= (SR + 1) / 2)
anchor_cells <- data.frame(Np = anchor_np, SR = anchor_sr, MM = anchor_mm)
cells <- bind_rows(mech_cells, anchor_cells) %>%
  distinct(Np, SR, MM) %>%
  arrange(Np, SR, MM) %>%
  mutate(role = ifelse(SR == anchor_sr & MM == anchor_mm,
                       ifelse(Np %in% mech_np, "anchor + mechanism", "anchor"),
                       "mechanism"))
cat(sprintf("Cells: %d (%d mechanism grid, %d anchor scales)\n",
            nrow(cells), nrow(mech_cells), length(anchor_np)))
cat(sprintf("Replicates per cell: %d (equal across all cells)\n", n_reps))

# --- Load the posterior draws for one cell ------------------------------------
# Draws are tagged with the seed set they came from so that between-draw-set
# variation can be reported: the effect below is small relative to it, and a
# single-seed study would not recover it reliably.
load_pool <- function(np, sr, mm) {
  sc <- job_grid %>%
    filter(target_Np == np, target_SR == sr, target_MM == mm) %>%
    pull(scenario_id) %>% unique()
  if (!length(sc)) stop("No scenario for Np=", np, " SR=", sr, " MM=", mm)
  fs <- list.files(chain_dir, full.names = TRUE,
                   pattern = sprintf("^scenario_%02d_seed_\\d+\\.rds$", sc[1]))
  if (!length(fs)) stop("No chain files for scenario ", sc[1])
  base_df <- NULL; draws <- list(); seedset <- integer(0)
  for (f in fs) {
    x <- readRDS(f)
    if (is.null(base_df)) base_df <- x$map_table
    draws   <- c(draws, x$posterior_draws)
    seedset <- c(seedset, rep(x$job$seed_set, length(x$posterior_draws)))
  }
  list(base_df = base_df, draws = draws, seedset = seedset)
}

pools <- setNames(
  lapply(seq_len(nrow(cells)), function(i)
    load_pool(cells$Np[i], cells$SR[i], cells$MM[i])),
  sprintf("%d_%d_%d", cells$Np, cells$SR, cells$MM))

tasks <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  data.frame(cells[i, c("Np","SR","MM")], rep = seq_len(n_reps), row.names = NULL)))
tasks$task_id <- seq_len(nrow(tasks))

# --- One replicate ------------------------------------------------------------
run_task <- function(k) {
  tk  <- tasks[k, ]
  key <- sprintf("%d_%d_%d", tk$Np, tk$SR, tk$MM)
  pool <- pools[[key]]
  if (!length(pool$draws)) return(NULL)

  set.seed(base_seed + tk$task_id)
  d_idx <- ((tk$rep - 1L) %% length(pool$draws)) + 1L

  df <- pool$base_df; df$MAP_Count <- pool$draws[[d_idx]]
  m  <- tryCatch(mateR2::mp_table_to_matrix(mp_table = df), error = function(e) NULL)
  if (is.null(m)) return(NULL)

  mixed <- mateR2::randomize_mating_structure(m, I = mixing_I)
  K     <- mateR2::brd.mat.fitness(mixed, min.fert = min_fert,
                                   max.fert = max_fert, type = fert_type)

  # Reproductive success per parent, from the FULL cohort (rows = males).
  rs_m <- rowSums(K); rs_f <- colSums(K)
  CVm  <- stats::sd(rs_m) / mean(rs_m)
  CVf  <- stats::sd(rs_f) / mean(rs_f)

  # Offspring-weighted mean parental RS = sum(rs^2)/sum(rs) = kbar (1 + CV^2).
  # This is the FULL-SIB-INCLUSIVE quantity the algebra above describes.
  ow_f <- sum(rs_f^2) / sum(rs_f)
  ow_m <- sum(rs_m^2) / sum(rs_m)

  # What sib.stats() reports EXCLUDES an offspring's own full-sib clutch, so it
  # is not the same quantity; the two are compared in the summary below.
  ss  <- mateR2::mat.sub.sample(K, num_offspring = n_sample)
  sub <- xtabs(off1 ~ dads + moms, data = ss)
  sub <- matrix(as.numeric(sub), nrow = nrow(sub), ncol = ncol(sub))
  s   <- mateR2::sib.stats(sub)

  data.frame(
    Np = tk$Np, SR = tk$SR, MM = tk$MM, rep = tk$rep,
    seed_set    = pool$seedset[d_idx],
    realized_SR = nrow(m) / ncol(m),
    CVm = CVm, CVf = CVf,
    ow_ratio_inclusive = ow_f / ow_m,
    MHS = s$mean_MHS, PHS = s$mean_PHS, diff = s$mean_MHS - s$mean_PHS
  )
}

# --- Execute ------------------------------------------------------------------
n_cores <- as.integer(Sys.getenv("MATER2_WORKERS", unset = "4"))
cat(sprintf("Running %d tasks on %d workers...\n", nrow(tasks), n_cores))
t0 <- Sys.time()

if (n_cores > 1) {
  cl <- makeCluster(n_cores); on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, library(mateR2))
  clusterExport(cl, c("tasks", "pools", "mixing_I", "min_fert", "max_fert",
                      "fert_type", "base_seed", "n_sample"), envir = environment())
  res <- parLapplyLB(cl, seq_len(nrow(tasks)), run_task)
} else {
  res <- lapply(seq_len(nrow(tasks)), run_task)
}
raw <- bind_rows(Filter(Negate(is.null), res))
cat(sprintf("Completed %d/%d replicates in %.2f minutes.\n",
            nrow(raw), nrow(tasks), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(raw, "data/outputs/sibling_asymmetry_replicates.rds")

# --- Per-cell summary and paired tests ----------------------------------------
summ <- raw %>%
  group_by(Np, SR, MM) %>%
  summarise(
    n_reps      = n(),
    realized_SR = mean(realized_SR),
    # NB: do not name any output column MHS/PHS/CVf/CVm -- summarise() evaluates
    # sequentially, so shadowing an input column breaks every later expression.
    CV_f = mean(CVf), CV_m = mean(CVm),
    mean_MHS = mean(MHS), mean_PHS = mean(PHS),
    ratio_excl  = mean(MHS) / mean(PHS),                 # what Section 4.3 reports
    ratio_incl  = mean(ow_ratio_inclusive),              # full-sib inclusive
    n_seed_sets = dplyr::n_distinct(seed_set),
    predicted   = mean(realized_SR) * (1 + mean(CVf)^2) / (1 + mean(CVm)^2),
    mean_diff   = mean(diff),
    ci_lo       = stats::t.test(MHS, PHS, paired = TRUE)$conf.int[1],
    ci_hi       = stats::t.test(MHS, PHS, paired = TRUE)$conf.int[2],
    p_value     = stats::t.test(MHS, PHS, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(incl_over_pred = ratio_incl / predicted,
         p_holm = p.adjust(p_value, "holm"))

write_csv(summ, "data/outputs/sibling_asymmetry_cells.csv")

# --- Table 1, exactly as the manuscript prints it -----------------------------
# Written from the same object the analysis produces, so the manuscript table
# and the pipeline output cannot drift apart. Five columns, N_P = 400,
# N_sample = 1000, rounded here rather than in the document.
table_01 <- summ %>%
  filter(Np == 400) %>%
  transmute(
    SR,
    MM,
    mean_maternal_HS = round(mean_MHS, 2),
    mean_paternal_HS = round(mean_PHS, 2),
    ratio            = round(ratio_excl, 3)
  ) %>%
  arrange(SR, MM)

write_csv(table_01, "data/outputs/table_01_sibling_asymmetry.csv")
cat("\n=== Table 1 (N_P = 400, N_sample =", n_sample, ") ===\n")
print(as.data.frame(table_01), row.names = FALSE)

cat("\n=== Per-cell summary (N_sample =", n_sample, ", equal n per cell) ===\n")
print(as.data.frame(summ %>% select(Np, SR, MM, CV_f, CV_m, mean_MHS, mean_PHS,
                                    ratio_excl, ratio_incl, predicted,
                                    incl_over_pred, p_value)),
      digits = 3, row.names = FALSE)

cat("\nThe algebra describes the FULL-SIB-INCLUSIVE ratio. Agreement with the\n")
cat("prediction (ratio_incl / predicted) should be near 1; ratio_excl departs\n")
cat("from it because half-sib counts exclude an offspring's own full-sib clutch.\n")

# --- Between-seed-set variation -----------------------------------------------
by_seed <- raw %>%
  group_by(Np, SR, MM, seed_set) %>%
  summarise(n = n(), ratio = mean(MHS) / mean(PHS), mean_diff = mean(diff),
            .groups = "drop")
write_csv(by_seed, "data/outputs/sibling_asymmetry_by_seedset.csv")

cat("\n=== Between-draw-set variation at the anchor scenario (SR = 2, MM = 2) ===\n")
print(as.data.frame(by_seed %>% filter(SR == anchor_sr, MM == anchor_mm)),
      digits = 3, row.names = FALSE)
cat("\nVariation between seed sets is comparable to the effect itself: a\n")
cat("single-seed study would not recover it reliably.\n")

# --- Pooled anchor analysis ---------------------------------------------------
cat("\n=== Anchor scenario (SR = 2, MM = 2) pooled across seed sets ===\n")
cat("    This is the Section 4.3 column: all five scales, equal effort.\n")
for (np in anchor_np) {
  x  <- raw %>% filter(Np == np, SR == anchor_sr, MM == anchor_mm)
  tt <- stats::t.test(x$MHS, x$PHS, paired = TRUE)
  cat(sprintf("  Np = %4d  n = %3d  MHS %6.3f  PHS %6.3f  ratio %5.3f  diff %7.3f  CI [%6.3f, %6.3f]  p = %.4f\n",
              np, nrow(x), mean(x$MHS), mean(x$PHS), mean(x$MHS) / mean(x$PHS),
              tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))
}
cat("\n[SUCCESS] -> data/outputs/sibling_asymmetry_{replicates.rds,cells.csv,by_seedset.csv}\n")
