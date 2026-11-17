# ======================================================================
# fig1_study_design.R — Figure 1
#
# Participant flow, kinship firewall and the data each analytic step may use.
#
# Inputs : derived/analysis_sets/*; results/guangzhou_predictions.csv
# Outputs: figures/png|pdf/fig1_design.*
# Paper  : Figure 1
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls())
library(data.table); library(ggplot2); library(patchwork)

source("config/paths.R")
AR   <- file.path(BASE, "derived/analysis_sets")
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)

cz <- readRDS(file.path(AR, "cohort_final.rds"))
tm <- fread(file = file.path(AR, "sample_sizes.csv"))
gz <- fread(file = file.path(RES, "guangzhou_predictions.csv"))
N <- list(univ = 6732, bmd = 6713, cand = 4424, kin = 4383, gwas = 4379,
          eval = 2289, bgi = 1794, bgi_f = 1781, disc = 1430, hold = 351, lc = 495,
          cross = 131, drop = 41, drop_lc = 13, gz0 = 499, gz = 492, qian = 28938)
obs <- list(univ = tm[nivel == 0 & grepl("los_qc", bloque), n],
            bmd  = tm[nivel == 1, n],
            cand = length(cz$gwas0), kin = length(cz$gwas),
            gwas = tm[grepl("analizados en el GWAS", bloque), n],
            eval = length(cz$evalset), bgi = length(cz$BGI),
            bgi_f = length(cz$discovery) + length(cz$holdout),
            disc = length(cz$discovery), hold = length(cz$holdout), lc = length(cz$LC),
            cross = nrow(cz$cross), drop = length(cz$drop), drop_lc = length(cz$drop_BGI_LC),
            gz0 = tm[grepl("Guangzhou", bloque), n], gz = nrow(gz), qian = 28938)
mal <- names(N)[unlist(N) != unlist(obs)]
if (length(mal)) stop("cifras que no coinciden con los archivos: ", paste(mal, collapse = ", "))

a1 <- cz$cross$A; a2 <- cz$cross$B
stopifnot(sum((a1 %in% cz$gwas & a2 %in% cz$evalset) | (a2 %in% cz$gwas & a1 %in% cz$evalset)) == 0,
          all(cz$drop %in% cz$gwas0), sum(a1 %in% cz$drop | a2 %in% cz$drop) == nrow(cz$cross))
k1 <- cz$kin$A; k2 <- cz$kin$B; bf <- c(cz$discovery, cz$holdout)
stopifnot(sum((k1 %in% bf & k2 %in% cz$LC) | (k2 %in% bf & k1 %in% cz$LC)) == 0,
          sum((k1 %in% cz$discovery & k2 %in% cz$holdout) |
              (k2 %in% cz$discovery & k1 %in% cz$holdout)) == 0)
cat("control OK: todas las cifras y las tres paredes coinciden con los archivos\n")
f <- function(x) format(x, big.mark = ",")

AZ <- "#2b6cb0"; NA_ <- "#d98c00"; VE <- "#117a65"; RO <- "#c0392b"; GR <- "grey40"

cajas <- rbindlist(list(
  list(30, 70, 90, 99, sprintf("Louisiana Osteoporosis Study (LOS)\ngenotyped after QC   n = %s", f(N$univ)), "white", AZ, "grey10"),
  list(30, 70, 79, 87, sprintf("with total hip BMD   n = %s", f(N$bmd)), "white", AZ, "grey10"),
  list(3, 41, 62, 71, sprintf("Without metagenomics: GWAS candidates\nn = %s", f(N$cand)), "white", AZ, "grey10"),
  list(3, 41, 47, 55, sprintf("after kinship pruning   n = %s", f(N$kin)), "white", AZ, "grey10"),
  list(3, 41, 29, 39, sprintf("GWAS sample   n = %s\nLOS PRS weights", f(N$gwas)), AZ, AZ, "white"),
  list(51, 97, 62, 71, sprintf("With metagenomics: evaluation set\nn = %s", f(N$eval)), "white", AZ, "grey10"),
  list(51, 77, 47, 55, sprintf("BGI platform   n = %s", f(N$bgi)), "white", AZ, "grey10"),
  list(80, 97, 47, 55, "LC Sciences\nplatform", "white", NA_, "grey10"),
  list(51, 63.5, 29, 39, sprintf("Modeling\ndiscovery\nn = %s", f(N$disc)), "#dbe6f3", AZ, "grey10"),
  list(64.5, 77, 29, 39, sprintf("Internal\nholdout\nn = %s", f(N$hold)), "white", AZ, "grey10"),
  list(80, 97, 29, 39, sprintf("LC subset\nn = %s", f(N$lc)), "#fbeed9", NA_, "grey10"),
  list(3, 41, 3, 16, sprintf("External GWAS: Qian et al. 2025\nUK Biobank Europeans, N = %s\nexternal PRS weights", f(N$qian)), "grey95", GR, "grey10"),
  list(51, 97, 3, 16, sprintf("Guangzhou cohort: peri/post-menopausal women\nn = %s, %s analysed\n(6 without PRS genotype, 1 missing alcohol)\nseparate genotype QC; never used for fitting",
                              f(N$gz0), f(N$gz)), "#dcefe9", VE, "grey10")))
setnames(cajas, c("x0","x1","y0","y1","txt","fill","col","tc"))

fl <- rbindlist(list(
  list(50, 90, 50, 87), list(50, 79, 50, 75), list(22, 75, 74, 75),
  list(22, 75, 22, 71), list(74, 75, 74, 71),
  list(22, 62, 22, 55), list(22, 47, 22, 39),
  list(64, 62, 64, 55), list(88.5, 62, 88.5, 55),
  list(64, 47, 64, 43), list(57.25, 43, 70.75, 43),
  list(57.25, 43, 57.25, 39), list(70.75, 43, 70.75, 39), list(88.5, 47, 88.5, 39)))
setnames(fl, c("x","y","xend","yend"))
fl[, punta := !(y == 75 & yend == 75) & !(y == 43 & yend == 43) & !(x == 50 & y == 79)]

notas <- data.table(
  x = c(52, 23, 23, 65, 65.5),
  y = c(77, 58.5, 43, 58.5, 45.2),
  txt = c(sprintf("\u2212%d without BMD", N$univ - N$bmd),
          sprintf("\u2212%d related to the\nevaluation set", N$drop),
          sprintf("\u2212%d anthropometric outliers", N$kin - N$gwas),
          sprintf("\u2212%d related to LC", N$drop_lc),
          "family-aware split"))

pa <- ggplot() +

  annotate("rect", xmin = 43.6, xmax = 48.4, ymin = 28, ymax = 72, fill = "#f6dcd8",
           colour = RO, linewidth = 0.4) +
  annotate("text", x = 46, y = 50, label = "kinship firewall", angle = 90, colour = RO,
           size = 2.6, fontface = "bold") +

  annotate("segment", x = 41.3, xend = 50.6, y = 31, yend = 31, colour = GR, linewidth = 0.45,
           arrow = arrow(length = unit(0.06, "inches"), type = "closed")) +
  annotate("text", x = 49.6, y = 23.6, hjust = 1, size = 2.05, colour = GR, lineheight = 0.95,
           label = "only effect sizes cross;\nno individual does") +
  annotate("segment", x = 22, xend = 22, y = 16, yend = 18.5, colour = GR, linewidth = 0.4) +
  annotate("segment", x = 22, xend = 50.6, y = 18.5, yend = 18.5, colour = GR, linewidth = 0.4) +
  annotate("segment", x = 50.6, xend = 50.6, y = 18.5, yend = 28.4, colour = GR, linewidth = 0.4,
           arrow = arrow(length = unit(0.05, "inches"), type = "closed")) +
  geom_rect(data = cajas, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = cajas$fill, colour = cajas$col, linewidth = 0.45) +
  geom_text(data = cajas, aes(x = (x0 + x1)/2, y = (y0 + y1)/2, label = txt),
            colour = cajas$tc, size = 2.35, lineheight = 0.95) +
  geom_segment(data = fl[(punta)], aes(x, y, xend = xend, yend = yend), colour = AZ, linewidth = 0.45,
               arrow = arrow(length = unit(0.05, "inches"), type = "closed")) +
  geom_segment(data = fl[!(punta)], aes(x, y, xend = xend, yend = yend), colour = AZ, linewidth = 0.45) +
  geom_text(data = notas, aes(x, y, label = txt), hjust = 0, size = 2.1, colour = "grey30",
            fontface = "italic") +

  annotate("text", x = 98, y = 22.0, hjust = 1, size = 2.05, colour = RO, lineheight = 1.0,
           label = sprintf(paste0("related pairs (KING, \u2265 3rd degree) after pruning:\n",
                                  "GWAS\u2013evaluation %d \u2192 0  |  BGI\u2013LC 0  |  discovery\u2013holdout 0"),
                           N$cross)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(2, 100), expand = FALSE) +
  labs(title = "a   Participant flow and kinship firewall") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(b = 3)))

COLS <- c(sprintf("GWAS\nsample\n%s", f(N$gwas)), sprintf("Modeling\ndiscovery\n%s", f(N$disc)),
          sprintf("Internal\nholdout\n%s", f(N$hold)), sprintf("LC\nsubset\n%s", f(N$lc)),
          sprintf("Guangzhou\ncohort\n%s", f(N$gz)), sprintf("External\nGWAS (Qian)\n%s", f(N$qian)))

FILAS <- list(
  "LOS genotype QC, PCA, kinship (no phenotype)"      = c(1,1,1,1,0,0),
  "GWAS effect sizes (PRS weights)"                    = c(1,0,0,0,0,1),
  "PRS p-value threshold selection"                    = c(0,1,0,0,0,0),
  "MRS construction and balance selection"             = c(0,1,0,0,0,0),
  "Model fitting (Elastic Net, XGBoost)"               = c(0,1,0,0,0,0),
  "dR2" = c(0,1,0,0,0,0),
  "Microbiome-covariate associations"                  = c(0,1,0,1,1,0),
  "Confirmation of frozen models"                      = c(0,0,2,0,0,0),
  "Transfer of frozen models"                          = c(0,0,0,1,1,0),
  "Interpretation of frozen models (permutation, SHAP)" = c(0,0,2,0,0,0))
mx <- rbindlist(lapply(seq_along(FILAS), function(i)
  data.table(fila = names(FILAS)[i], col = COLS, v = FILAS[[i]])))
mx[, fila := factor(fila, levels = rev(names(FILAS)))]
mx[, col := factor(col, levels = COLS)]
mx[, lab := fcase(v == 1, "\u2713", v == 2, "\u2713\u2020", default = "\u2014")]

pb <- ggplot(mx, aes(col, fila)) +
  geom_tile(aes(fill = factor(v)), colour = "grey80", linewidth = 0.35) +
  geom_text(aes(label = lab, colour = factor(v)), size = 2.6) +
  scale_fill_manual(values = c(`0` = "white", `1` = AZ, `2` = AZ), guide = "none") +
  scale_colour_manual(values = c(`0` = "grey65", `1` = "white", `2` = "white"), guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +

  scale_y_discrete(expand = c(0, 0),
                   labels = function(x) lapply(x, function(v) if (v == "dR2")
                     quote(Delta*R^2~"estimation (5" %*% "5 CV) and permutation test") else v)) +
  labs(title = "b   Which data each step may see",
       caption = paste0("\u2020 The holdout was opened three times for estimation (free alpha; alpha fixed at 0.25; ",
                        "seventh specification added),\nand once to interpret the frozen models. ",
                        "No step that sets weights, thresholds, balances or model parameters uses it.")) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(b = 4)),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.caption = element_text(size = 6.4, colour = "grey38", hjust = 0, lineheight = 1.2),
        axis.title = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 6.6, lineheight = 0.95, colour = "grey15"),
        axis.text.y = element_text(size = 7.0, colour = "grey15"))

p <- wrap_elements(full = pa) / pb + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(
    title = "Weights, thresholds and models are fixed on one side of a kinship firewall and evaluated on the other",
    theme = theme(plot.title = element_text(size = 9.5, hjust = 0)))

ggsave(file.path(FIG,"png","fig1_design.png"), p, width = 8.4, height = 8.6, dpi = 300)
ggsave(file.path(FIG,"pdf","fig1_design.pdf"), p, width = 8.4, height = 8.6, device = cairo_pdf)
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
