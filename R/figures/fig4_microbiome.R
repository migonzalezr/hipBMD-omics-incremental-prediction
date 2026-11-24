# ======================================================================
# fig4_microbiome.R — Figure 4
# #
# Balances, species, whole-community composition and host-covariate correlations before and after adjustment.
# #
# Inputs : results (microbiome association files)
# Outputs: figures/png, figures/pdf
# Paper  : Figure 4
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table)
library(ggplot2)
library(patchwork)

source("config/paths.R")
RES <- file.path(BASE, "results")
FIG <- Sys.getenv("FIG3_OUTPUT_DIR", unset = file.path(BASE, "figures"))
for (subdir in c("png", "pdf")) {
  dir.create(file.path(FIG, subdir), recursive = TRUE, showWarnings = FALSE)
}

alpha <- 0.05
blue <- "#216091"
orange <- "#A4531A"
groups <- c("Unadjusted", "Age + sex + weight")

figure_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 7)),
    axis.text = element_text(colour = "grey20"),
    axis.title = element_text(size = 9),
    plot.margin = margin(6, 9, 7, 6)
  )

adjustments <- c(
  "sin ajuste" = "Unadjusted", "solo edad" = "Age", "solo sexo" = "Sex",
  "solo peso" = "Weight", "solo talla" = "Height", "solo BMI" = "BMI",
  "edad + sexo" = "Age + sex", "edad + peso" = "Age + weight",
  "sexo + peso" = "Sex + weight", "peso + talla" = "Weight + height",
  "edad + sexo + peso (esq.B)" = "Age + sex + weight"
)
balances <- fread(file = file.path(RES, "balance_decomposition_covariates.csv"))
stopifnot(setequal(balances$ajuste, names(adjustments)))
balances[, adjustment := factor(adjustments[ajuste], levels = rev(adjustments))]
pa <- ggplot(balances, aes(n_FDR05, adjustment)) +
  geom_col(fill = blue, width = 0.62) +
  geom_text(aes(label = n_FDR05), nudge_x = 0.35, hjust = 0, size = 3.1) +
  scale_x_continuous(breaks = seq(0, 15, 5), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "a    Compositional balances", x = "Balances with BH-adjusted p < 0.05", y = NULL) +
  figure_theme +
  theme(axis.ticks.y = element_blank())

species <- fread(file = file.path(RES, "species_association.csv"))
stopifnot(nrow(species) == 428L,
          isTRUE(all.equal(species$FDR_crudo, p.adjust(species$p_crudo, "BH"))),
          isTRUE(all.equal(species$FDR_ajust, p.adjust(species$p_ajust, "BH"))))
species_p <- melt(species, measure.vars = c("FDR_crudo", "FDR_ajust"),
                  variable.name = "model", value.name = "p_adj")
species_p[, model := factor(model, labels = groups)]
at_threshold <- species_p[, .(n_species = sum(p_adj < alpha)), by = model]
line_labels <- setNames(sprintf("%s (%d microbial species)", at_threshold$model, at_threshold$n_species),
                        at_threshold$model)
ticks <- c(1e-9, 1e-6, 1e-3, alpha, 1)
pb <- ggplot(species_p, aes(p_adj, colour = model, linetype = model)) +
  geom_vline(xintercept = alpha, colour = "grey55", linewidth = 0.35, linetype = "dotted") +
  stat_ecdf(aes(y = after_stat(ecdf) * nrow(species)), pad = FALSE, linewidth = 0.75) +
  geom_point(data = at_threshold, aes(x = alpha, y = n_species), size = 2) +
  scale_x_log10(breaks = ticks, labels = c("1e-9", "1e-6", "0.001", "0.05", "1")) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400)) +
  scale_colour_manual(values = c(blue, orange), labels = line_labels, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = line_labels, name = NULL) +
  labs(title = "b    Individual species", x = "BH-adjusted p value (log scale)",
       y = "Cumulative number of species") +
  figure_theme +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification.inside = c(0, 1), legend.key.width = unit(13, "pt"),
        legend.background = element_blank(), legend.text = element_text(size = 8))

permanova <- fread(file = file.path(RES, "diversity_permanova.csv"))
community <- rbind(
  permanova[analisis == "BMD sin ajustar" & termino == "Model"],
  permanova[analisis == "BMD ajustado" & termino == "y"]
)
stopifnot(nrow(community) == 2L)
community[, model := factor(groups, levels = rev(groups))]
community[, percent := 100 * R2]
pc <- ggplot(community, aes(percent, model)) +
  geom_col(fill = blue, width = 0.45) +
  geom_text(aes(label = sprintf("%.2f%%  (p = %.3f)", percent, p)),
            nudge_x = 0.025, hjust = 0, size = 3.1) +
  scale_x_continuous(limits = c(0, 1.10), breaks = seq(0, 1.0, 0.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "c    Whole-community composition",
       x = "Bray–Curtis variation explained by hip BMD (%)", y = NULL) +
  figure_theme +
  theme(axis.ticks.y = element_blank())

host_names <- c("gender_dummy" = "Sex", "Height_cm" = "Height",
                "Weight_kg" = "Weight", "Age" = "Age", "BMI" = "BMI")
host <- fread(file = file.path(RES, "balance_decomposition_correlations.csv"))
stopifnot(setequal(host[variable != "y", variable], names(host_names)))
top_balance <- balances[ajuste == "sin ajuste", top]
host_long <- melt(host[variable != "y"], measure.vars = c("r_con_balance_top", "r_con_BMD"),
                  variable.name = "measure", value.name = "r")
host_long[, measure := factor(measure, labels = c(sprintf("Top balance (%s)", top_balance), "Hip BMD"))]
host_long[, variable := factor(host_names[variable], levels = rev(host_names))]
pd <- ggplot(host_long, aes(measure, variable, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.3) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(low = "#B65F2E", mid = "#F8F8F8", high = "#4284BC",
                       midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1), name = "Pearson r") +
  guides(fill = guide_colourbar(direction = "vertical", barwidth = unit(8, "pt"),
                                barheight = unit(81, "pt"), title.position = "top")) +
  labs(title = "d    Correlations with host covariates", x = NULL, y = NULL) +
  figure_theme +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.box.spacing = unit(12, "pt"),
        legend.margin = margin(0, 8, 0, 0),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8))

figure <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1.7, 1))
ggsave(file.path(FIG, "png", "fig3_microbiome_main.png"), figure,
       width = 8.3, height = 5.65, dpi = 300, bg = "white")
pdf_device <- if (Sys.info()[["sysname"]] == "Darwin") {
  function(filename, ...) quartz(type = "pdf", file = filename, ...)
} else {
  cairo_pdf
}
ggsave(file.path(FIG, "pdf", "fig3_microbiome_main.pdf"), figure,
       width = 8.3, height = 5.65, device = pdf_device, bg = "white")
