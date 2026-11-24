# ======================================================================
# figS7_host_variables.R — Figure S7
# #
# Host variables associated with community composition and Shannon diversity.
# #
# Inputs : results/host_variables_*.csv
# Outputs: figures/png, figures/pdf
# Paper  : Figure S7
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
rm(list = ls())
library(data.table); library(ggplot2); library(patchwork)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)

COL <- c("Discovery (BGI)"="#2b6cb0","LC subset"="#d98c00",
         "Guangzhou cohort"="#117a65")
KEEP <- c("Age","Sex","Weight","Height","Hip BMD")
AJ   <- "Hip BMD (adjusted)"
FOCO <- c("Hip BMD", AJ)

pm <- fread(file = file.path(RES,"host_variables_permanova.csv"))[modelo=="marginal" & variable %in% KEEP]
sh <- fread(file = file.path(RES,"host_variables_shannon.csv"))[modelo=="marginal" & variable %in% KEEP]
sen <- fread(file = file.path(RES,"adjustment_sensitivity.csv"))[conjunto=="B"]
stopifnot(nrow(sen) == 3)

REN <- c("LC subset" = "LC subset", "Guangzhou cohort" = "Guangzhou cohort")
for (x in list(pm, sh, sen)) x[cohorte %in% names(REN), cohorte := unname(REN[cohorte])]

pm <- rbind(pm, sen[, .(cohorte, variable=AJ, modelo="ajustado", n,
                        R2=permanova_R2, p=permanova_p, estimable=TRUE)], fill=TRUE)
sh <- rbind(sh, sen[, .(cohorte, variable=AJ, modelo="ajustado", n,
                        beta=shannon_beta, lo=shannon_lo, hi=shannon_hi,
                        p=shannon_p, estimable=TRUE)], fill=TRUE)

ord <- pm[cohorte == "Discovery (BGI)"][order(-R2), variable]
prep <- function(d){
  d <- copy(d)
  d[, variable := factor(variable, levels = rev(ord))]
  d[, coh := factor(cohorte, levels = names(COL))]
  d[, sig := !is.na(p) & p < 0.05]
  d[(estimable)]
}
pm <- prep(pm); sh <- prep(sh)

bandas <- data.table(variable = factor(FOCO, levels = rev(ord)))
DG <- position_dodge(width = 0.62)

EJE <- rev(ord)
ETQ_EJE <- lapply(EJE, function(v) {
  if (v == AJ)        quote(atop(bold("Hip BMD"), bold("(adjusted)")))
  else if (v %in% FOCO) bquote(bold(.(v)))
  else                 v
})
GRP <- rev(names(COL))

th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.title = element_blank())

pa <- ggplot(pm, aes(R2, variable, colour = coh, group = factor(cohorte, levels = GRP))) +
  geom_tile(data = bandas, aes(x = 0.01, y = variable), width = 1, height = 0.9,
            fill = "#f2e7df", colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_linerange(aes(xmin = 0, xmax = R2), position = DG, linewidth = 0.35, alpha = 0.55) +
  geom_point(aes(fill = ifelse(sig, as.character(coh), "ns")), position = DG,
             shape = 21, size = 1.9, stroke = 0.6) +
  scale_fill_manual(values = c(COL, ns = "white"), guide = "none") +
  scale_colour_manual(values = COL) +
  scale_y_discrete(limits = EJE, labels = ETQ_EJE) +
  coord_cartesian(xlim = c(0, 0.0215), expand = FALSE, ylim = c(0.45, length(EJE) + 0.55)) +
  annotate("text", x = 0.0205, y = length(EJE) - 0.30, label = "single-sex cohorts",
           hjust = 1, size = 2.2, colour = "grey45", fontface = "italic") +
  labs(title = "a   Whole community composition",
       x = expression("PERMANOVA "*R^2*"  (variable alone)"), y = NULL) +
  th + theme(legend.position = c(0.97, 0.05), legend.justification = c(1,0),
             legend.text = element_text(size = 6.4),
             legend.key.height = unit(0.52,"lines"),
             legend.background = element_blank(),
             legend.key = element_blank())

pb <- ggplot(sh, aes(beta, variable, colour = coh, group = factor(cohorte, levels = GRP))) +
  geom_tile(data = bandas, aes(x = 0, y = variable), width = 10, height = 0.9,
            fill = "#f2e7df", colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = DG, linewidth = 0.45, alpha = 0.8) +
  geom_point(aes(fill = ifelse(sig, as.character(coh), "ns")), position = DG,
             shape = 21, size = 1.9, stroke = 0.6) +
  scale_fill_manual(values = c(COL, ns = "white"), guide = "none") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_y_discrete(limits = EJE) +
  coord_cartesian(xlim = c(-0.21, 0.235), expand = FALSE, ylim = c(0.45, length(EJE) + 0.55)) +
  labs(title = "b   Shannon diversity",
       x = expression("standardized "*beta*"  (SD of Shannon per SD)"), y = NULL) +
  th + theme(axis.text.y = element_blank())

p <- (pa | pb) + plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    title = paste("Gut community composition tracks sex and height;",
                  "its association with BMD does not survive adjustment"),
    subtitle = paste0("Bray-Curtis over the 306 species shared by the three cohorts, 999 permutations. ",
                      "Each variable tested alone;\nthe bottom row is BMD with age, sex, weight ",
                      "and height in the model. Filled markers, p < 0.05."),
    theme = theme(plot.title = element_text(size = 9.2, hjust = 0),
                  plot.subtitle = element_text(size = 6.4, colour = "grey38",
                                               hjust = 0, lineheight = 1.4)))

ggsave(file.path(FIG,"png","figS2_host_variables_microbiome.png"), p, width = 8.4, height = 3.6, dpi = 300)
ggsave(file.path(FIG,"pdf","figS2_host_variables_microbiome.pdf"), p, width = 8.4, height = 3.6, device = cairo_pdf)
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
