# ==============================================================================
# Script 08: Build Figure 2 (Stacked Panels 2a and 2b)
# Project: mateR2 Manuscript
# Goal: Build a single vertically stacked composite figure containing:
#       - (a) Mate-Pair Summary Table (Compressed CBS State Space c_b)
#       - (b) Expanded Bipartite Mating Matrix with Demographic Proofs
# ==============================================================================

library(dplyr)
library(tidyr)
library(gt)
library(magick)

# Ensure the figures directory exists
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# ==============================================================================
# 1. BUILD PANEL (a): MATE-PAIR SUMMARY TABLE
# ==============================================================================

mp_summary_df <- data.frame(
  Block_ID      = c("1", "2", "3", "4"),
  Block_Type    = c("1M : 1F", "1M : 2F", "2M : 1F", "3M : 3F"),
  Males_m       = c(1, 1, 2, 3),
  Females_f     = c(1, 2, 1, 3),
  Count_c_b     = c(1, 1, 1, 1),
  stringsAsFactors = FALSE
)

fig2a_table <- mp_summary_df %>%
  gt() %>%
  cols_label(
    Block_ID      = md("**Block ($b$)**"),
    Block_Type    = md("**Mating Motif**"),
    Males_m       = md("**Males ($m_b$)**"),
    Females_f     = md("**Females ($f_b$)**"),
    Count_c_b     = md("**Count ($c_b$)**")
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_header(
    title = md("**a) Compressed Mate-Pair Summary Table**"),
    subtitle = md("State vector $\\vec{c} = [c_1, c_2, c_3, c_4]$ parameterizing CBS space")
  ) %>%
  tab_style(
    style = cell_fill(color = "#f2f7fb"),
    locations = cells_body(columns = Count_c_b)
  ) %>%
  opt_table_lines(extent = "none") %>%
  tab_options(
    heading.align = "left",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.title.font.size = px(15),
    heading.subtitle.font.size = px(12),
    table.background.color = "white"
  )

# ==============================================================================
# 2. BUILD PANEL (b): EXPANDED MATRIX & MARGINAL DEMOGRAPHICS
# ==============================================================================

Nm <- 7 # Males (Rows)
Nf <- 7 # Females (Columns)

binary_matrix <- matrix(0, nrow = Nm, ncol = Nf)
rownames(binary_matrix) <- paste0("M", 1:Nm)
colnames(binary_matrix) <- paste0("F", 1:Nf)

# Define blocks matching Panel (a)
binary_matrix[1, 1] <- 1        # Motif 1: 1M:1F
binary_matrix[2, 2:3] <- 1      # Motif 2: 1M:2F
binary_matrix[3:4, 4] <- 1      # Motif 3: 2M:1F
binary_matrix[5:7, 5:7] <- 1    # Motif 4: 3M:3F

matrix_data <- as.data.frame(binary_matrix)
matrix_data$Male_ID <- rownames(matrix_data)
matrix_data <- matrix_data %>% select(Male_ID, everything())

# Row Marginals (d_i^m)
matrix_data <- matrix_data %>%
  rowwise() %>%
  mutate(`Male Marginals (d_i^m)` = sum(c_across(starts_with("F")))) %>%
  ungroup()

# Column Marginals (d_j^f)
female_marginals <- matrix_data %>%
  summarise(across(starts_with("F"), sum)) %>%
  mutate(
    Male_ID = "**Marginals** ($d_j^f$)",
    `Male Marginals (d_i^m)` = sum(c_across(starts_with("F"))) # Total Edges (E)
  ) %>%
  select(Male_ID, everything())

plot_df <- bind_rows(matrix_data, female_marginals)

# Math proofs
total_edges <- female_marginals$`Male Marginals (d_i^m)`
Np <- Nm + Nf
SR <- Nm / Nf
Mean_Mates <- ((total_edges / Nm) + (total_edges / Nf)) / 2

fig2b_table <- plot_df %>%
  gt(rowname_col = "Male_ID") %>%
  fmt_markdown(columns = stub()) %>%
  tab_spanner(
    label = md("**Female Partners ($j$)**"),
    columns = starts_with("F")
  ) %>%
  cols_label(
    `Male Marginals (d_i^m)` = md("**Marginals**<br>($d_i^m$)")
  ) %>%
  cols_align(
    align = "center",
    columns = c(starts_with("F"), `Male Marginals (d_i^m)`)
  ) %>%
  tab_header(
    title = md("**b) Expanded Bipartite Mating Matrix**"),
    subtitle = "Mapping individual node degrees to population demographic targets"
  ) %>%
  tab_stubhead(label = md("**Males ($i$)**")) %>%
  text_transform(
    locations = cells_body(columns = starts_with("F")),
    fn = function(x) {
      ifelse(x == "0", paste0("<span style='color: #D3D3D3;'>", x, "</span>"), 
             paste0("<b>", x, "</b>"))
    }
  ) %>%
  tab_style(
    style = cell_fill(color = "#e6f2ff"),
    locations = cells_body(columns = `Male Marginals (d_i^m)`)
  ) %>%
  tab_style(
    style = cell_fill(color = "#e6f2ff"),
    locations = cells_body(rows = Nm + 1)
  ) %>%
  tab_style(
    style = cell_fill(color = "#cce5ff"),
    locations = cells_body(columns = `Male Marginals (d_i^m)`, rows = Nm + 1)
  ) %>%
  tab_source_note(
    source_note = md(sprintf(
      "**Legend:** $d_i^m$ = Male Mating Success; $d_j^f$ = Female Mating Success<br>
      <br>
      **Demographic Proofs Derived from Network Marginals:**<br>
      • **Total Parents ($N_P$):** $N_M + N_F = \\sum c_b m_b + \\sum c_b f_b = %d + %d = %d$<br>
      • **Sex Ratio ($SR$):** $N_M / N_F = %d / %d = %.2f$<br>
      • **Total Mating Edges ($E$):** $\\sum d_i^m = \\sum d_j^f = \\sum c_b (m_b \\cdot f_b) = %d$<br>
      • **Mean Mates ($MM$):** $\\frac{1}{2}\\left(\\frac{E}{N_M} + \\frac{E}{N_F}\\right) = \\frac{E \\cdot N_P}{2 \\cdot N_M \\cdot N_F} = \\frac{%d \\cdot %d}{2 \\cdot %d \\cdot %d} = %.2f$",
      Nm, Nf, Np, Nm, Nf, SR, total_edges, total_edges, Np, Nm, Nf, Mean_Mates
    ))
  ) %>%
  opt_table_lines(extent = "none") %>%
  opt_row_striping(row_striping = FALSE) %>%
  tab_options(
    heading.align = "left",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.title.font.size = px(15),
    heading.subtitle.font.size = px(12),
    source_notes.font.size = px(12),
    table.background.color = "white"
  )

# ==============================================================================
# 3. RENDER & STACK VERTICALLY INTO A SINGLE COMPOSITE IMAGE
# ==============================================================================

temp_a <- tempfile(fileext = ".png")
temp_b <- tempfile(fileext = ".png")

# High-resolution rendering
gtsave(fig2a_table, filename = temp_a, zoom = 3)
gtsave(fig2b_table, filename = temp_b, zoom = 3)

img_a <- image_read(temp_a)
img_b <- image_read(temp_b)

# Standardize canvas widths: align Panel (a) to the left (West)
target_width <- max(image_info(img_a)$width, image_info(img_b)$width)

img_a_aligned <- image_extent(img_a, geometry = paste0(target_width, "x", image_info(img_a)$height), gravity = "West", color = "white")
img_b_aligned <- image_extent(img_b, geometry = paste0(target_width, "x", image_info(img_b)$height), gravity = "West", color = "white")

# Add vertical spacer padding between panels
img_a_padded <- image_border(img_a_aligned, "white", "0x20")
img_b_padded <- image_border(img_b_aligned, "white", "0x20")

# Stack vertically
combined_fig2 <- image_append(c(img_a_padded, img_b_padded), stack = TRUE)

# Save the unified vertical publication graphic
output_path <- "figures/Figure_02_Matrix_Algebra.png"
image_write(combined_fig2, path = output_path, format = "png")

cat(sprintf("\n[Success] Left-aligned Figure 2 saved to:\n -> %s\n", output_path))
