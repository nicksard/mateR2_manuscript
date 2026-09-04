# ==============================================================================
# Script 08: Build Figure 2 (Stacked Panels 2a and 2b)
# Project: mateR2 Manuscript
# Goal: A single vertically stacked composite figure containing:
#       - (a) Mate-Pair Summary Table (compressed CBS state space, c_b)
#       - (b) Expanded bipartite mating matrix with demographic proofs
#
# CHANGED: rebuilt in ggplot2. The previous version rendered gt tables to PNG
# with gtsave() and stacked them with magick, which required TWO pieces of
# system software outside renv's reach: the ImageMagick C++ library, and a
# headless Chrome install (gt renders HTML and screenshots it). Neither can be
# restored by renv, so the script died for anyone without both. Figures 3-6 are
# already plain R; this makes the whole figure pipeline plain R.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggtext)   # real subscripts in the proofs block; pure R, no system deps

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

ink        <- "grey15"
fill_count <- "#f2f7fb"   # panel (a) count column
fill_marg  <- "#e6f2ff"   # panel (b) marginals
fill_corn  <- "#cce5ff"   # panel (b) corner (total edges)
zero_grey  <- "grey80"

base_tbl <- theme_void(base_size = 9) +
  theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0),
        plot.subtitle = element_text(size = 8.5, colour = "grey30", hjust = 0),
        # element_markdown() renders the proofs with real subscripts; plain
        # element_text() prints "N_M" and "c_b" literally, which misstates
        # notation Section 2.1 defines carefully.
        plot.caption  = ggtext::element_markdown(size = 8, colour = ink,
                                                 hjust = 0, lineheight = 1.5),
        plot.margin   = margin(6, 10, 6, 10))

# ==============================================================================
# 1. PANEL (a): MATE-PAIR SUMMARY TABLE
# ==============================================================================
mp <- data.frame(
  Block_ID   = c("1", "2", "3", "4"),
  Block_Type = c("1M : 1F", "1M : 2F", "2M : 1F", "3M : 3F"),
  Males_m    = c(1, 1, 2, 3),
  Females_f  = c(1, 2, 1, 3),
  Count_c_b  = c(1, 1, 1, 1)
)

hdr_a <- c("Block~(italic(b))", "Mating~motif", "Males~(italic(m)[italic(b)])",
           "Females~(italic(f)[italic(b)])", "Count~(italic(c)[italic(b)])")
a_long <- mp %>%
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = "col", values_to = "val",
               values_transform = as.character) %>%
  mutate(colx = match(col, names(mp)),
         val  = ifelse(grepl("^[0-9]+$", val), val, val))

p_a <- ggplot(a_long, aes(colx, -row)) +
  # count column highlight
  geom_tile(data = subset(a_long, colx == 5),
            fill = fill_count, colour = NA, width = 0.96, height = 0.9) +
  geom_text(aes(label = val), size = 3.1, colour = ink) +
  annotate("text", x = seq_along(hdr_a), y = 0, label = hdr_a, parse = TRUE,
           fontface = "bold", size = 3.1, colour = ink) +
  annotate("segment", x = 0.45, xend = 5.55, y = -0.55, yend = -0.55, colour = ink) +
  annotate("segment", x = 0.45, xend = 5.55, y =  0.55, yend =  0.55, colour = ink) +
  annotate("segment", x = 0.45, xend = 5.55, y = -4.55, yend = -4.55, colour = ink) +
  scale_x_continuous(limits = c(0.35, 5.7)) +
  scale_y_continuous(limits = c(-4.9, 0.9)) +
  labs(title = "a) Compressed Mate-Pair Summary Table",
       subtitle = expression("State vector "*c==group("[",list(c[1],c[2],c[3],c[4]),"]")*
                             " parameterizing CBS space")) +
  base_tbl

# ==============================================================================
# 2. PANEL (b): EXPANDED MATRIX & MARGINAL DEMOGRAPHICS
# ==============================================================================
Nm <- 7; Nf <- 7
B <- matrix(0, Nm, Nf, dimnames = list(paste0("M", 1:Nm), paste0("F", 1:Nf)))
B[1, 1]     <- 1   # 1M:1F
B[2, 2:3]   <- 1   # 1M:2F
B[3:4, 4]   <- 1   # 2M:1F
B[5:7, 5:7] <- 1   # 3M:3F

row_marg <- rowSums(B); col_marg <- colSums(B); E <- sum(B)
Np <- Nm + Nf; SR <- Nm / Nf; MM <- ((E / Nm) + (E / Nf)) / 2

cells <- expand.grid(i = 1:Nm, j = 1:Nf) %>%
  mutate(v = as.vector(B[cbind(i, j)]))

p_b <- ggplot() +
  # matrix body
  geom_tile(data = cells, aes(j, -i), fill = "white", colour = "grey88", linewidth = 0.3) +
  geom_text(data = cells, aes(j, -i, label = v,
                              fontface = ifelse(v == 1, "bold", "plain"),
                              colour = ifelse(v == 1, ink, zero_grey)), size = 3) +
  # row marginals
  geom_tile(data = data.frame(i = 1:Nm, v = row_marg),
            aes(Nf + 1, -i), fill = fill_marg, colour = "grey88", linewidth = 0.3) +
  geom_text(data = data.frame(i = 1:Nm, v = row_marg),
            aes(Nf + 1, -i, label = v), size = 3, fontface = "bold", colour = ink) +
  # column marginals
  geom_tile(data = data.frame(j = 1:Nf, v = col_marg),
            aes(j, -(Nm + 1)), fill = fill_marg, colour = "grey88", linewidth = 0.3) +
  geom_text(data = data.frame(j = 1:Nf, v = col_marg),
            aes(j, -(Nm + 1), label = v), size = 3, fontface = "bold", colour = ink) +
  # corner: total edges
  annotate("tile", x = Nf + 1, y = -(Nm + 1), fill = fill_corn,
           colour = "grey70", linewidth = 0.3) +
  annotate("text", x = Nf + 1, y = -(Nm + 1), label = E, size = 3,
           fontface = "bold", colour = ink) +
  # labels
  annotate("text", x = 1:Nf, y = 0, label = colnames(B), size = 2.9, colour = ink) +
  annotate("text", x = Nf + 1, y = 0, label = "d[i]^m", parse = TRUE,
           size = 2.9, fontface = "bold", colour = ink) +
  annotate("text", x = 0, y = -(1:Nm), label = rownames(B), size = 2.9, colour = ink) +
  annotate("text", x = 0, y = -(Nm + 1), label = "d[j]^f", parse = TRUE,
           size = 2.9, fontface = "bold", colour = ink) +
  annotate("text", x = (Nf + 1) / 2, y = 0.85, label = "Female partners (j)",
           size = 3.1, fontface = "bold", colour = ink) +
  annotate("text", x = -0.85, y = -(Nm + 1) / 2, label = "Males (i)", angle = 90,
           size = 3.1, fontface = "bold", colour = ink) +
  scale_colour_identity() +
  coord_fixed(clip = "off") +
  scale_x_continuous(limits = c(-1.2, Nf + 1.8)) +
  scale_y_continuous(limits = c(-(Nm + 1.7), 1.3)) +
  labs(title = "b) Expanded Bipartite Mating Matrix",
       subtitle = "Mapping individual node degrees to population demographic targets",
       caption = sprintf(paste0(
         "Legend: row marginals = male mating success (*d*<sub>*i*</sub><sup>*m*</sup>); ",
         "column marginals = female mating success (*d*<sub>*j*</sub><sup>*f*</sup>)<br><br>",
         "**Demographic proofs derived from network marginals:**<br>",
         "&bull; Total parents: *N*<sub>P</sub> = *N*<sub>M</sub> + *N*<sub>F</sub> = ",
         "&Sigma; *c*<sub>*b*</sub>*m*<sub>*b*</sub> + &Sigma; *c*<sub>*b*</sub>*f*<sub>*b*</sub> ",
         "= %d + %d = **%d**<br>",
         "&bull; Sex ratio: *SR* = *N*<sub>M</sub> / *N*<sub>F</sub> = %d / %d = **%.2f**<br>",
         "&bull; Total mating edges: *E* = &Sigma; *d*<sub>*i*</sub><sup>*m*</sup> = ",
         "&Sigma; *d*<sub>*j*</sub><sup>*f*</sup> = ",
         "&Sigma; *c*<sub>*b*</sub>(*m*<sub>*b*</sub> &middot; *f*<sub>*b*</sub>) = **%d**<br>",
         "&bull; Mean mates: *MM* = &frac12;(*E*/*N*<sub>M</sub> + *E*/*N*<sub>F</sub>) = ",
         "(*E* &middot; *N*<sub>P</sub>) / (2 &middot; *N*<sub>M</sub> &middot; *N*<sub>F</sub>) = ",
         "(%d &middot; %d) / (2 &middot; %d &middot; %d) = **%.2f**"),
         Nm, Nf, Np, Nm, Nf, SR, E, E, Np, Nm, Nf, MM)) +
  base_tbl

# ==============================================================================
# 3. STACK
# ==============================================================================
fig2 <- p_a / p_b + plot_layout(heights = c(1, 3.0))
out <- "figures/Figure_02_Worked_Example.png"
ggsave(out, fig2, width = 7.2, height = 8.0, dpi = 400, bg = "white")
cat("[SUCCESS] ->", out, "\n")
