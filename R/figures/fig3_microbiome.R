# ======================================================================
# fig3_microbiome.R — Figure 3
#
# Microbiome associations before and after adjustment at the balance, species and community levels.
#
# Inputs : results of scripts 14-15
# Outputs: figures/png|pdf/fig3_microbiome_main.*
# Paper  : Figure 3
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls())
library(data.table); library(ggplot2); library(patchwork)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)
AZ <- "#2b6cb0"; RO <- "#c0392b"
th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

dec <- fread(file = file.path(RES,"balance_decomposition_covariates.csv"))
cor <- fread(file = file.path(RES,"balance_decomposition_correlations.csv"))
esp <- fread(file = file.path(RES,"species_association.csv"))
pmv <- fread(file = file.path(RES,"diversity_permanova.csv"))

EN <- c("sin ajuste"="No adjustment","solo edad"="Age","solo sexo"="Sex","solo peso"="Weight",
        "solo talla"="Height","solo BMI"="BMI","edad + sexo"="Age + sex",
        "edad + peso"="Age + weight","sexo + peso"="Sex + weight",
        "peso + talla"="Weight + height","edad + sexo + peso (esq.B)"="Age + sex + weight")
dec[, lab := EN[ajuste]][, nlf := -log10(FDR_min)]
dec[, g := fifelse(lab=="No adjustment", 0L,
           fifelse(lab %in% c("BMI","Age","Height","Sex","Weight"), 1L, 2L))]
setorder(dec, g, -nlf)
dec[, y := seq_len(.N) + 0.85*(g - min(g))]
dec[, cero := n_FDR05 == 0]

pa <- ggplot(dec, aes(nlf, -y)) +
  geom_vline(xintercept = -log10(0.05), colour = "grey55", linewidth = 0.3, linetype = "22") +
  geom_segment(aes(x = 0, xend = nlf, yend = -y), colour = "grey75", linewidth = 0.35) +
  geom_point(aes(colour = cero), size = 2.3) +
  geom_text(aes(x = nlf + 0.4, label = n_FDR05, colour = cero), hjust = 0, size = 2.3,
            show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE`=AZ, `TRUE`=RO),
                      labels = c("at least one balance at FDR < 0.05","no balance at FDR < 0.05"), name = NULL) +
  scale_y_continuous(breaks = -dec$y, labels = dec$lab) +
  scale_x_continuous(limits = c(-0.35, 15.2), expand = c(0,0)) +
  annotate("text", x = -log10(0.05)+0.25, y = -max(dec$y), label = "FDR = 0.05",
           hjust = 0, size = 2.2, colour = "grey45") +
  labs(title = "a   No single covariate removes it",
       x = expression("Strongest balance,  " * -log[10] * "FDR"), y = NULL) +
  th + theme(legend.position = c(0.72, 0.45), legend.background = element_blank(),
             legend.key = element_blank(), legend.title = element_text(size = 6.8),
             legend.text = element_text(size = 6.8))

conv <- data.table(
  lab = factor(c("Compositional\nbalances","Whole community\n(PERMANOVA)","Individual\nspecies"),
               levels = c("Individual\nspecies","Whole community\n(PERMANOVA)","Compositional\nbalances")),
  p0 = c(dec[lab=="No adjustment", FDR_min],
         pmv[analisis=="BMD sin ajustar" & termino=="Model", p],
         min(esp$FDR_crudo)),
  p1 = c(dec[lab=="Age + sex + weight", FDR_min],
         pmv[analisis=="BMD ajustado" & termino=="y", p],
         min(esp$FDR_ajust)),
  ann = c("15 \u2192 0 balances", "R\u00b2 0.63% \u2192 0.09%", "76 \u2192 0 species"))
conv[, `:=`(x0 = -log10(p0), x1 = -log10(p1))]

pb <- ggplot(conv, aes(y = lab)) +
  geom_vline(xintercept = -log10(0.05), colour = "grey55", linewidth = 0.3, linetype = "22") +
  geom_segment(aes(x = x0, xend = x1, yend = lab), colour = "grey68", linewidth = 0.4,
               arrow = arrow(length = unit(0.055,"in"), type = "closed")) +
  geom_point(aes(x = x0), colour = AZ, size = 2.2) +
  geom_point(aes(x = x1), colour = RO, size = 2.2) +
  geom_text(aes(x = x0 + 0.6, label = ann), hjust = 0, size = 2.2, colour = "grey30") +
  scale_x_continuous(limits = c(-0.5, 17.5), expand = c(0,0)) +
  labs(title = "b   Three analyses, same collapse",
       x = expression("Strongest signal,  " * -log[10] * "(FDR or " * italic(P) * ")"), y = NULL) +
  th

EN2 <- c("gender_dummy"="Sex","Height_cm"="Height","Weight_kg"="Weight","Age"="Age","BMI"="BMI")
cor <- cor[variable != "y"][, lab := EN2[variable]]
setorder(cor, -r_con_balance_top)
cor[, y := seq_len(.N)][, es_bmi := lab == "BMI"]

pc <- ggplot(cor, aes(r_con_balance_top, -y)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = r_con_balance_top, yend = -y), colour = "grey75", linewidth = 0.35) +
  geom_point(aes(colour = es_bmi), size = 2.3) +
  geom_text(aes(x = r_con_balance_top + ifelse(r_con_balance_top >= 0, 0.022, -0.022),
                label = sprintf("%.2f", r_con_balance_top),
                hjust = ifelse(r_con_balance_top >= 0, 0, 1)), size = 2.3, colour = "grey30") +
  scale_colour_manual(values = c(`FALSE`=AZ, `TRUE`=RO), guide = "none") +
  scale_y_continuous(breaks = -cor$y, labels = cor$lab) +
  scale_x_continuous(limits = c(-0.30, 0.52), expand = c(0,0)) +
  labs(title = "c   What the balance tracks",
       x = "Correlation with strongest balance", y = NULL) + th

p <- pa | (pb / pc)
p <- p + plot_layout(widths = c(1.38, 1)) +
  plot_annotation(
    title = "The microbiome-BMD association is not independent of age, sex and weight",
    theme = theme(plot.title = element_text(size = 9.5, hjust = 0)))

ggsave(file.path(FIG,"png","fig3_microbiome_main.png"), p, width = 8.1, height = 5.2, dpi = 300)
ggsave(file.path(FIG,"pdf","fig3_microbiome_main.pdf"), p, width = 8.1, height = 5.2, device = cairo_pdf)
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
