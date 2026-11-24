# ======================================================================
# figS6_microbiome_by_cohort.R — Figure S6
# #
# Gut microbiome composition across cohorts.
# #
# Inputs : results/diversity_*.csv
# Outputs: figures/png, figures/pdf
# Paper  : Figure S6
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
rm(list = ls())
library(data.table); library(ggplot2); library(patchwork); library(ggforce)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)

COL <- c("Discovery (BGI)"="#2b6cb0","LC subset"="#d98c00",
         "Guangzhou cohort"="#117a65")
SHORT <- c("Discovery (BGI)"="Discovery\n(BGI)","LC subset"="LC\nsubset",
           "Guangzhou cohort"="Guangzhou\ncohort")
th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))

uni <- fread(file = file.path(RES,"diversity_species_universe.csv"))
sh  <- fread(file = file.path(RES,"diversity_shannon.csv"))
pc  <- fread(file = file.path(RES,"diversity_pcoa.csv"))
pmv <- fread(file = file.path(RES,"diversity_permanova.csv"))

REN <- c("LC subset" = "LC subset", "Guangzhou cohort" = "Guangzhou cohort")
for (x in list(sh, pc)) x[cohorte %in% names(REN), cohorte := unname(REN[cohorte])]

n_amb  <- uni[en_LC & en_China, .N]
n_lc   <- uni[en_LC & !en_China, .N]
n_ch   <- uni[!en_LC & en_China, .N]
n_solo <- uni[!en_LC & !en_China, .N]
circ <- data.table(x0 = c(0.40, 0.60), y0 = c(0.46, 0.46), r = c(0.255, 0.255),
                   coh = c("LC subset","Guangzhou cohort"))
pa <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.06, ymax = 0.86,
           fill = NA, colour = "grey55", linewidth = 0.3) +
  ggforce::geom_circle(data = circ, aes(x0 = x0, y0 = y0, r = r, fill = coh),
                       colour = NA, alpha = 0.30, show.legend = FALSE) +
  scale_fill_manual(values = COL) +
  annotate("text", x = 0.50, y = 0.46, label = n_amb,  size = 3.6, fontface = "bold") +
  annotate("text", x = 0.26, y = 0.46, label = n_lc,   size = 2.8) +
  annotate("text", x = 0.745, y = 0.46, label = n_ch,  size = 2.8) +
  annotate("text", x = 0.50, y = 0.145, size = 2.5, colour = "grey30",
           label = sprintf("%d in discovery only", n_solo)) +
  annotate("text", x = 0.30, y = 0.785, size = 2.5, fontface = "bold",
           colour = COL[["LC subset"]],
           label = sprintf("LC subset\n%d", uni[(en_LC), .N])) +
  annotate("text", x = 0.71, y = 0.785, size = 2.5, fontface = "bold",
           colour = COL[["Guangzhou cohort"]],
           label = sprintf("Guangzhou cohort\n%d", uni[(en_China), .N])) +
  annotate("text", x = 0.02, y = 0.95, hjust = 0, size = 2.5, colour = "grey35",
           label = sprintf("%d discovery species", nrow(uni))) +
  coord_fixed(xlim = c(0,1), ylim = c(0,1.02), expand = FALSE) +
  labs(title = "a   Species shared across cohorts") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))

sh[, coh := factor(cohorte, levels = names(COL))]
sh[, sx := substr(sexo, 1, 1)]
wil <- wilcox.test(shannon ~ sexo, data = sh[coh == "Discovery (BGI)"])$p.value
brk <- data.table(coh = factor("Discovery (BGI)", levels = names(COL)),
                  x = c(1,1,2,2), y = c(4.63,4.72,4.72,4.63))
pb <- ggplot(sh, aes(x = sx, y = shannon, colour = coh)) +
  geom_jitter(width = 0.18, size = 0.35, alpha = 0.22, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45,
               linewidth = 0.3, colour = "grey10", show.legend = FALSE) +
  geom_path(data = brk, aes(x, y), colour = "grey45", linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(data = data.table(coh = factor("Discovery (BGI)", levels = names(COL)),
                              x = 1.5, y = 4.80,
                              lab = sprintf("Wilcoxon P = %.2f", wil)),
            aes(x, y, label = lab), size = 2.3, colour = "grey30",
            vjust = 0, inherit.aes = FALSE) +
  facet_grid(~ coh, scales = "free_x", space = "free_x",
             labeller = labeller(coh = SHORT)) +
  scale_colour_manual(values = COL) +
  coord_cartesian(ylim = c(0.60, 5.30)) +
  labs(title = "b   Alpha diversity", x = NULL, y = "Shannon index") +
  th + theme(strip.background = element_blank(), strip.text = element_text(size = 7),
             panel.spacing = unit(0.5, "lines"))

pc[, coh := factor(cohorte, levels = names(COL))]
r2 <- pmv[analisis == "cohortes" & termino == "Model", R2]
pp <- pmv[analisis == "cohortes" & termino == "Model", p]
pcc <- ggplot(pc, aes(PCo1, PCo2, colour = coh)) +
  geom_point(size = 0.40, alpha = 0.38) +
  scale_colour_manual(values = COL, labels = gsub("\n"," ",SHORT[names(COL)]), name = NULL) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.4, size = 2.2, colour = "grey30",
           label = sprintf("PERMANOVA by cohort\nR\u00b2 = %.3f,  P = %.3f", r2, pp)) +
  labs(title = "c   Community structure",
       x = sprintf("PCo1  (%.1f%%)", pc$var_PCo1[1]),
       y = sprintf("PCo2  (%.1f%%)", pc$var_PCo2[1])) +
  th + theme(legend.position = "none")

etq <- data.table(
  x   = c(-0.315, -0.325, 0.505),
  y   = c(-0.370,  0.408, 0.300),
  lab = c("Discovery (BGI)","Guangzhou\ncohort","LC subset"),
  coh = factor(c("Discovery (BGI)","Guangzhou cohort","LC subset"),
               levels = names(COL)),
  hj  = c(0, 0, 1))
pcc <- pcc +
  geom_text(data = etq, aes(x, y, label = lab, colour = coh, hjust = hj),
            size = 2.5, fontface = "bold", vjust = 1, lineheight = 0.95,
            inherit.aes = FALSE, show.legend = FALSE) +
  coord_cartesian(xlim = c(-0.345, 0.525), ylim = c(-0.465, 0.415))

p <- (pa | pb | pcc) + plot_layout(widths = c(1, 1, 1.15)) +
  plot_annotation(
    title = "Microbiome composition differs markedly between cohorts",
    subtitle = paste("Cohorts were sequenced on different platforms, so between-cohort",
                     "differences in diversity are confounded with platform and pipeline."),
    theme = theme(plot.title = element_text(size = 9.5, hjust = 0),
                  plot.subtitle = element_text(size = 6.6, colour = "grey38", hjust = 0)))

ggsave(file.path(FIG,"png","figS3_microbiome_by_cohort.png"), p, width = 8.2, height = 3.4, dpi = 300)
ggsave(file.path(FIG,"pdf","figS3_microbiome_by_cohort.pdf"), p, width = 8.2, height = 3.4, device = cairo_pdf)
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
