#!/usr/bin/env Rscript
##############################################################################
# Visualize LDSC genetic correlation (RA vs immune-mediated diseases) as a
# one-column heatmap of rg with a p-value column.
#   Usage: Rscript plot_ldsc_rg.R stats/rg_summary.csv fig/rg_heatmap.png
##############################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

args   <- commandArgs(trailingOnly = TRUE)
INFILE <- if (length(args) >= 1) args[1] else "stats/rg_summary.csv"
OUTFILE<- if (length(args) >= 2) args[2] else "fig/rg_heatmap.png"

dt <- fread(INFILE)
n_all <- nrow(dt)
bonf  <- 0.05 / n_all                    # Bonferroni over all traits tested
dt[, sig := ifelse(p < bonf, "*", "")]
dt <- dt[!is.na(rg) & !is.na(p)]         # drop traits without an rg estimate

# clean trait labels
dt[, trait_label := gsub("_", " ", trait)]
dt[, rg_cap := pmin(pmax(rg, -1), 1)]    # cap to [-1,1] for display
dt[, p_label := ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                        formatC(p, format = "f", digits = 4))]
dt[, trait_label := factor(trait_label, levels = trait_label[order(rg_cap)])]

p1 <- ggplot(dt, aes(x = "rg", y = trait_label, fill = rg_cap)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f%s", rg_cap, sig)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = expression(r[g])) +
  labs(x = "", y = "") + theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10, face = "bold"),
        panel.grid = element_blank())

p2 <- ggplot(dt, aes(x = "p", y = trait_label)) +
  geom_text(aes(label = p_label, fontface = ifelse(p < bonf, "bold", "plain")), size = 2.8) +
  labs(x = "", y = "") + theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold"),
        panel.grid = element_blank())

p <- p1 + p2 + plot_layout(widths = c(3, 1.5)) +
  plot_annotation(
    title = "Genetic correlation with RA (immune-mediated diseases)",
    subtitle = paste0("* Bonferroni (p<", formatC(bonf, format = "e", digits = 1), ")"),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")))

h <- max(4, min(nrow(dt) * 0.25 + 2, 48))
ggsave(OUTFILE, p, width = 7, height = h, dpi = 300)
ggsave(sub("\\.[^.]+$", ".pdf", OUTFILE), p, width = 7, height = h)
cat(sprintf("wrote %s (%d traits, bonf=%.2g)\n", OUTFILE, nrow(dt), bonf))
