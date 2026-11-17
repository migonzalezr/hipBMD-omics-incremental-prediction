# ======================================================================
# fig2_prs_ceiling.R — Figure 2
#
# PRS threshold selection, theoretical ceiling from SNP heritability, and permutation null.
#
# Inputs : results/prs_*_threshold_selection.csv; prs_permutations.csv
# Outputs: figures/png|pdf/fig2_prs_ceiling.*
# Paper  : Figure 2
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls())
library(data.table); library(ggplot2); library(patchwork)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)

CL <- "#2b6cb0"
CQ <- "#c0392b"
CC <- "#117a65"
H2   <- 0.319
R2C  <- 0.380193

ESC  <- 1 - R2C
N_LOS <- 4379; N_QIAN <- 28938

los <- fread(file = file.path(RES,"prs_los_threshold_selection.csv"))[, fuente := "LOS GWAS (N = 4,379)"]
qia <- fread(file = file.path(RES,"prs_external_threshold_selection.csv"))[, fuente := "Qian 2025 (N = 28,938)"]
cur <- rbind(los, qia)
cur[, x := as.numeric(sub("PT_","",PT))]
cur[, fuente := factor(fuente, levels = c("LOS GWAS (N = 4,379)","Qian 2025 (N = 28,938)"))]
COLF <- setNames(c(CL,CQ), levels(cur$fuente))
gan <- cur[, .SD[which.max(dR2)], by = fuente]
stopifnot(nrow(gan) == 2)

per <- fread(file = file.path(RES,"prs_permutations.csv"))
res <- fread(file = file.path(RES,"prs_permutations_summary.csv"))
per[, et := fifelse(score == "LOS_PT0.01", "LOS", "Qian")]
res[, et := fifelse(score == "LOS_PT0.01", "LOS", "Qian")]
per[, et := factor(et, levels = c("LOS","Qian"))]; res[, et := factor(et, levels = c("LOS","Qian"))]
COLE <- c(LOS = CL, Qian = CQ)

th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))

pa <- ggplot(cur, aes(x, dR2, colour = fuente)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_errorbar(aes(ymin = dR2 - sd, ymax = dR2 + sd), width = 0, linewidth = 0.35) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.2) +

  geom_point(data = gan, shape = 21, size = 3.4, fill = NA, stroke = 0.6) +
  geom_text(data = gan, aes(label = format(variantes, big.mark = ",")),

            vjust = c(-1.9, -2.7), hjust = c(0.9, -0.12), size = 2.1, show.legend = FALSE) +
  scale_x_log10(breaks = 10^c(-6,-4,-2,0),
                labels = expression(10^-6, 10^-4, 10^-2, 1)) +
  scale_colour_manual(values = COLF, name = NULL) +
  coord_cartesian(ylim = c(-0.0034, 0.0078)) +
  labs(title = "a   Threshold selection", x = "P-value threshold",
       y = expression(Delta*R^2*" added to the clinical model")) +

  th + theme(legend.position = c(0.015, 0.99), legend.justification = c(0, 1),
             legend.text = element_text(size = 6.2),
             legend.key.height = unit(0.5,"lines"),
             legend.background = element_blank())

techo <- function(N, M, h2 = H2) ESC * h2 / (1 + M/(N*h2))
NN <- 10^seq(3, 6.2, length.out = 300)
banda <- data.table(N = NN, lo = techo(NN, 1e5), hi = techo(NN, 5e4), mid = techo(NN, 6e4))
pts <- data.table(N = c(N_LOS, N_QIAN), R2 = c(gan[1, dR2], gan[2, dR2]),
                  et = factor(c("LOS","Qian"), levels = c("LOS","Qian")))

pb <- ggplot() +
  geom_ribbon(data = banda, aes(N, ymin = lo, ymax = hi), fill = "grey82", alpha = 0.55) +
  geom_line(data = banda, aes(N, mid), colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = H2 * ESC, colour = "grey45", linewidth = 0.35, linetype = "22") +
  geom_hline(yintercept = R2C, colour = CC, linewidth = 0.55) +

  annotate("segment", x = N_QIAN, xend = N_QIAN, y = gan[2, dR2], yend = techo(N_QIAN, 6e4),
           colour = CQ, linewidth = 0.3,
           arrow = arrow(length = unit(0.05,"inches"), ends = "both", type = "open")) +
  annotate("text", x = 3.4e4, y = 0.013, hjust = 0, size = 2.05, colour = CQ,
           lineheight = 0.95, label = "transfer cost\nEurope \u2192 Louisiana") +
  geom_point(data = pts, aes(N, R2, colour = et), size = 1.5) +
  scale_colour_manual(values = COLE, guide = "none") +

  annotate("text", x = 1.3e3, y = 0.50, hjust = 0, size = 2.2, colour = CC,
           label = "clinical model, R\u00b2 = 0.38") +
  annotate("text", x = 1.3e3, y = 0.145, hjust = 0, size = 2.2, colour = "grey35",
           label = "SNP heritability, total scale = 0.20  (asymptote)") +
  annotate("text", x = 4.6e3, y = 0.0036, hjust = 0, vjust = 1, size = 2.2, colour = CL, label = "LOS") +
  annotate("text", x = 3.4e4, y = 0.0042, hjust = 0, vjust = 0.5, size = 2.2, colour = CQ, label = "Qian") +
  annotate("text", x = 1.3e4, y = 0.00262, hjust = 0.5, size = 2.05, colour = "grey35", label = "observed \u0394R\u00b2 (cross-validated)") +
  annotate("text", x = 1.15e5, y = 0.0062, hjust = 0, size = 2.05, colour = "grey35",
           label = "theoretical ceiling") +
  annotate("text", x = 1.15e5, y = 0.0044, hjust = 0, size = 2.05, colour = "grey35",
           parse = TRUE, label = "group('(', M[eff]~'50,000-100,000', ')')") +
  scale_x_log10(breaks = 10^(3:6), labels = c("1,000","10,000","100,000","1,000,000")) +
  scale_y_log10(breaks = c(0.003, 0.01, 0.03, 0.1, 0.3),
                labels = c("0.003","0.01","0.03","0.1","0.3")) +

  coord_cartesian(xlim = c(1e3, 1.6e6), ylim = c(0.0025, 0.75)) +
  labs(title = "b   What a PRS can reach, at best",
       x = "GWAS sample size", y = "Variance in hip BMD explained") + th

YMAX <- 1.12 * max(per[, max(hist(dR2_nulo, breaks = seq(-0.0020, 0.0024, by = 0.0002),
                                  plot = FALSE)$counts), by = et]$V1)
stopifnot(per[, all(dR2_nulo > -0.0020 & dR2_nulo < 0.0024)])
pc <- ggplot(per, aes(dR2_nulo, fill = et)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +

  geom_histogram(breaks = seq(-0.0020, 0.0024, by = 0.0002), alpha = 0.40, colour = NA) +
  geom_segment(data = res, aes(x = dR2_observado, xend = dR2_observado, y = 16, yend = 1,
                               colour = et),
               arrow = arrow(length = unit(0.05,"inches")), linewidth = 0.45,
               inherit.aes = FALSE) +
  geom_text(data = res, aes(x = dR2_observado, y = 18, colour = et,
                            label = sprintf("observed\nz = %.1f", z)),
            size = 2.05, vjust = 0, lineheight = 0.95, inherit.aes = FALSE) +
  geom_text(data = res, aes(x = -0.0024, y = 0.96 * YMAX, colour = et, label = et),
            size = 2.3, hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ et, ncol = 1) +
  scale_fill_manual(values = COLE, guide = "none") +
  scale_colour_manual(values = COLE, guide = "none") +
  scale_x_continuous(breaks = c(-0.002, 0, 0.002, 0.004)) +
  scale_y_continuous(breaks = c(0, 20, 40)) +

  coord_cartesian(xlim = c(-0.0026, 0.0053), ylim = c(0, YMAX)) +
  labs(title = "c   Score permutation",
       x = expression(Delta*R^2*" under the null"), y = "permutations") +
  th + theme(strip.text = element_blank(), strip.background = element_blank())

p <- (pa | pb | pc) + plot_layout(widths = c(1, 1.1, 0.85)) +
  plot_annotation(
    title = paste("A hip-BMD polygenic score cannot reach what age, sex, weight",
                  "and height already explain"),

    subtitle = paste0("Modeling discovery set (BGI), n = 1,430; 5 \u00d7 5-fold cross-validation. ",
                      "PRS weights come from the independent GWAS sample (N = 4,379),\n",
                      "separated by a kinship firewall. The association is real - 0 of 150 ",
                      "permutations of the score reached the observed value - but it adds 0.43% of variance."),
    theme = theme(plot.title = element_text(size = 9.3, hjust = 0),
                  plot.subtitle = element_text(size = 6.4, colour = "grey38",
                                               hjust = 0, lineheight = 1.4)))

ggsave(file.path(FIG,"png","fig2_prs_ceiling.png"), p, width = 8.7, height = 3.5, dpi = 300)
ggsave(file.path(FIG,"pdf","fig2_prs_ceiling.pdf"), p, width = 8.7, height = 3.5, device = cairo_pdf)
cat(sprintf("techo LOS  (M=60k): %.5f   observado %.5f  (%.0f%% del techo)\n",
            techo(N_LOS,6e4), gan[1,dR2], 100*gan[1,dR2]/techo(N_LOS,6e4)))
cat(sprintf("techo Qian (M=60k): %.5f   observado %.5f  (%.0f%% del techo)\n",
            techo(N_QIAN,6e4), gan[2,dR2], 100*gan[2,dR2]/techo(N_QIAN,6e4)))
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
