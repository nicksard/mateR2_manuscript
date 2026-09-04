# ==============================================================================
# Script 07: Build Figure 1 (Bipartite Networks and Constraint Matrices)
# Project: mateR2 Manuscript
# Goal: Create a visual, pedagogical 7x7 toy matrix mapping diverse mating 
#       blocks (1:1, 2:1, 1:2, and 3:3 mass-spawning). This script generates 
#       a 4-panel figure explicitly linking the bipartite network topology 
#       to the underlying matrix algebra (Phase 1 & Phase 2) and saves it.
# ==============================================================================

# Load required libraries
library(mateR2)
library(patchwork)
library(ggplot2)

#setting workign directory
# Run from the repository root (open the .Rproj, or `Rscript scripts/07_...R`).
# No setwd(): a hard-coded absolute path makes the script unrunnable for anyone
# who is not the author, which is the first thing a reviewer hits.

# ==========================================
# 1. CONSTRUCT THE STRUCTURAL PARENT POOL
# ==========================================
# Explicitly design a 7x7 matrix to capture the requested diversity
# of mating blocks: 1:1, 2:1, 1:2, and a 3:3 cluster.

Nm <- 7 # Males (Rows)
Nf <- 7 # Females (Columns)
binary_matrix <- matrix(0, nrow = Nm, ncol = Nf)
rownames(binary_matrix) <- paste0("M", 1:Nm)
colnames(binary_matrix) <- paste0("F", 1:Nf)

# Block 1: Strict Monogamy (1:1)
binary_matrix[1, 1] <- 1

# Block 2: Polygyny (1 Male, 2 Females)
binary_matrix[2, 2:3] <- 1

# Block 3: Polyandry (2 Males, 1 Female)
binary_matrix[3:4, 4] <- 1

# Block 4: Complex Polygynandry / Mass Spawning (3:3)
binary_matrix[5:7, 5:7] <- 1

# ==========================================
# 2. ASSIGN FECUNDITY (Phase 2 Weights)
# ==========================================
# Generate random offspring counts (1 to 15) only for successful matings
set.seed(42) # For reproducible paper exports
fecundity_weights <- matrix(sample(1:15, Nm * Nf, replace = TRUE), nrow = Nm, ncol = Nf)

# Multiply the binary structure by the weights to get the final realized matrix
weighted_matrix <- binary_matrix * fecundity_weights

# ==========================================
# 3. GENERATE THE 4 INDIVIDUAL PLOTS 
# ==========================================

# Top Row: Bipartite Networks 
plot_net_binary <- plot_bipartite_network(
  binary_matrix, 
  title = "Phase 1: Binary Topology\n(Mating Success)", 
  is_weighted = FALSE
)

plot_net_weighted <- plot_bipartite_network(
  weighted_matrix, 
  title = "Phase 2: Fecundity Allocation\n(Reproductive Success)", 
  is_weighted = TRUE
)

# Bottom Row: Constraint Matrices 
plot_mat_binary <- plot_realized_matrix(
  binary_matrix, 
  title = "Structural Mating Constraints", 
  fill_label = "Mated"
)

plot_mat_weighted <- plot_realized_matrix(
  weighted_matrix, 
  title = "Realized Breeding Matrix", 
  fill_label = "Offspring"
)

# ==========================================
# 4. STITCH INTO 2x2 GRID & SAVE TO PROJECT
# ==========================================
# Using patchwork: (Top Row) / (Bottom Row)
fig_1_schematic <- (plot_net_binary | plot_net_weighted) / 
  (plot_mat_binary | plot_mat_weighted) +
  plot_annotation(tag_levels = 'a',tag_suffix = ")") 

# Render the final figure in the RStudio viewer
print(fig_1_schematic)

# Export the figure to the project's 'figures' directory
ggsave(
  filename = "figures/Figure_01_Bipartite_Networks.png",
  plot = fig_1_schematic,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)
