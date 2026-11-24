# ======================================================================
# fig2_prs_thresholds.R — Figure 2
# #
# Delta R2 over the clinical base model across nine P-value thresholds (a) and score-permutation reference (b).
# #
# Inputs : results/prs_*_threshold_selection.csv, results/prs_permutations.csv
# Outputs: figures/png, figures/pdf
# Paper  : Figure 2
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
rm(list = ls())
library(data.table); library(ggplot2); library(patchwork)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- Sys.getenv("FIG_OUT", unset = file.path(BASE, "figures"))
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)

CL <- "#2b6cb0"
CQ <- "#c0392b"
CC <- "#117a65"
H2   <- 0.319
R2C  <- fread(file.path(RES, "prs_los_threshold_selection.csv"))$R2_base[1]
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

sens <- cur[fuente == levels(fuente)[2] & PT %in% c("PT_1e-6","PT_1e-4")]
stopifnot(nrow(sens) == 2)
RX <- 0.21; RY <- 0.00042; tt <- seq(0, 2*pi, length.out = 120)
sens_circ <- sens[, .(cx = 10^(log10(x) + RX*cos(tt)), cy = dR2 + RY*sin(tt)), by = PT]

pa <- ggplot(cur, aes(x, dR2, colour = fuente)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_errorbar(aes(ymin = dR2 - sd, ymax = dR2 + sd), width = 0, linewidth = 0.35) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.2) +
  geom_point(data = gan, shape = 21, size = 3.4, fill = NA, stroke = 0.6) +
  geom_path(data = sens_circ, aes(cx, cy, group = PT), colour = CQ, linewidth = 0.45,
            linetype = "22", inherit.aes = FALSE) +
  geom_text(data = sens, aes(x, dR2 + sd, label = format(variantes, big.mark = ",")),
            colour = CQ, size = 2.1, vjust = -0.7, inherit.aes = FALSE) +
  geom_text(data = gan, aes(label = format(variantes, big.mark = ",")),
            vjust = c(-1.9, -2.7), hjust = c(0.9, -0.12), size = 2.1, show.legend = FALSE) +
  scale_x_log10(breaks = 10^c(-6,-4,-2,0),
                labels = expression(10^-6, 10^-4, 10^-2, 1)) +
  scale_colour_manual(values = COLF, name = NULL) +
  coord_cartesian(ylim = c(-0.0034, 0.0078)) +
  labs(tag = "a", x = "P-value threshold",
       y = expression(Delta*R^2*" added to the clinical model")) +
  th + theme(legend.position = c(0.015, 0.99), legend.justification = c(0, 1),
             legend.text = element_text(size = 6.2),
             legend.key.height = unit(0.5,"lines"),
             legend.background = element_blank())


bin_width <- 0.0002
bin_limits <- range(per$dR2_nulo)
bin_breaks <- seq(floor(bin_limits[1] / bin_width) * bin_width,
                  ceiling(bin_limits[2] / bin_width) * bin_width, by = bin_width)
YMAX <- 1.12 * max(per[, max(hist(dR2_nulo, breaks = bin_breaks, plot = FALSE)$counts), by = et]$V1)
left_limit <- min(bin_breaks) - 0.0004
right_limit <- max(max(bin_breaks), res$dR2_observado) + 0.0006
pc <- ggplot(per, aes(dR2_nulo, fill = et)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_histogram(breaks = bin_breaks, alpha = 0.40, colour = NA) +
  geom_segment(data = res, aes(x = dR2_observado, xend = dR2_observado, y = 16, yend = 1,
                               colour = et),
               arrow = arrow(length = unit(0.05,"inches")), linewidth = 0.45,
               inherit.aes = FALSE) +
  geom_text(data = res, aes(x = dR2_observado, y = 18, colour = et,
                            label = sprintf("observed\nz = %.1f", z)),
            size = 2.05, vjust = 0, lineheight = 0.95, inherit.aes = FALSE) +
  geom_text(data = res, aes(x = left_limit + 0.0002, y = 0.96 * YMAX, colour = et, label = et),
            size = 2.3, hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ et, ncol = 1) +
  scale_fill_manual(values = COLE, guide = "none") +
  scale_colour_manual(values = COLE, guide = "none") +
  scale_x_continuous(breaks = c(-0.002, 0, 0.002, 0.004)) +
  scale_y_continuous(breaks = c(0, 20, 40)) +
  coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0, YMAX)) +
  labs(tag = "b",
       x = expression(Delta*R^2*" under the null"), y = "permutations") +
  th + theme(strip.text = element_blank(), strip.background = element_blank())

p <- (pa | pc) + plot_layout(widths = c(1.15, 1))

ggsave(file.path(FIG,"png","fig2_prs_ceiling.png"), p, width = 6.6, height = 3.5, dpi = 300)
pdf_device <- if (Sys.info()[["sysname"]] == "Darwin") {
  function(filename, ...) quartz(type = "pdf", file = filename, ...)
} else {
  cairo_pdf
}
ggsave(file.path(FIG,"pdf","fig2_prs_ceiling.pdf"), p, width = 6.6, height = 3.5, device = pdf_device)
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
