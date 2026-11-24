# ======================================================================
# figS9_model_interpretation.R — Figure S9
# #
# Elastic Net coefficients and XGBoost SHAP values in the internal holdout.
# #
# Inputs : results/interpretation_*.csv
# Outputs: figures/png, figures/pdf
# Paper  : Figure S9
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
rm(list = ls())
library(data.table); library(ggplot2); library(patchwork)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- Sys.getenv("FIGS6_OUTPUT_DIR", unset = file.path(BASE, "figures"))
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)
SP <- "S5q_base_PRSqian_MRS"
GRIS <- "#5b7fa6"; PRSC <- "#c0392b"; MRSC <- "#d98c00"

co <- fread(file = file.path(RES,"interpretation_coefficients.csv"))[especificacion == SP]
pm <- fread(file = file.path(RES,"interpretation_permutation.csv"))[especificacion == SP]
sh <- fread(file = file.path(RES,"interpretation_shap.csv"))[especificacion == SP]

NM <- c(Weight_kg="Weight", Age="Age", Height_cm="Height",
        Race_AfricanYes="African American", Race_CaucasianYes="White", Race_AsianYes="Asian",
        gender_dummyMale="Sex (male)", gPC1="Genetic PC1", gPC2="Genetic PC2",
        PRS_qian="PRS", MRS="MRS", Exercise_Regular_catYes="Exercise",
        Smoke_catYes="Smoking", Sun_Exposure_catYes="Sun exposure",
        Alcohol_Drinking_catYes="Alcohol", Cheese_Eat_catYes="Cheese",
        Milk_Drink_catYes="Milk", Yogurt_Eat_catYes="Yogurt")

d <- merge(merge(pm, co[, .(variable, coef_estandarizado)], by = "variable"),
           sh[, .(variable, shap_medio_abs)], by = "variable")
d[, lab := fifelse(is.na(NM[variable]), variable, NM[variable])]
setorder(d, -perdida_R2_EN)
d[, lab := factor(lab, levels = rev(lab))]
d[, col := fifelse(variable == "PRS_qian", PRSC, fifelse(variable == "MRS", MRSC, GRIS))]

th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())
sin_y <- theme(axis.text.y = element_blank())

pa <- ggplot(d, aes(coef_estandarizado, lab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = coef_estandarizado, yend = lab),
               colour = "grey78", linewidth = 0.35) +
  geom_point(aes(colour = I(col)), size = 1.9) +
  labs(tag = "a", x = expression("g/cm"^2*" per SD"), y = NULL) + th

pb <- ggplot(d, aes(perdida_R2_EN, lab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = perdida_R2_EN, yend = lab),
               colour = "grey78", linewidth = 0.35) +
  geom_point(aes(colour = I(col), shape = "Elastic Net"), size = 1.9) +
  geom_point(aes(x = perdida_R2_XGB, colour = I(col), shape = "XGBoost"),
             size = 1.6, fill = NA, stroke = 0.5) +
  scale_shape_manual(values = c(`Elastic Net` = 16, XGBoost = 22), name = NULL) +
  labs(tag = "b",
       x = expression("loss in holdout "*R^2*" when permuted"), y = NULL) +
  th + sin_y +
  theme(legend.position = c(0.97, 0.06), legend.justification = c(1,0),
        legend.text = element_text(size = 6.6), legend.key.height = unit(0.5,"lines"))

pc <- ggplot(d, aes(shap_medio_abs, lab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = shap_medio_abs, yend = lab),
               colour = "grey78", linewidth = 0.35) +
  geom_point(aes(colour = I(col)), size = 1.9) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.05))) +
  labs(tag = "c", x = "mean |SHAP|  (XGBoost)", y = NULL) + th + sin_y

p <- (pa | pb | pc) + plot_layout(widths = c(1, 1.25, 1))

ggsave(file.path(FIG,"png","figS6_variable_importance.png"), p, width = 8.3, height = 4.3, dpi = 300)
pdf_device <- if (Sys.info()[["sysname"]] == "Darwin") {
  function(filename, ...) quartz(type = "pdf", file = filename, ...)
} else {
  cairo_pdf
}
ggsave(file.path(FIG,"pdf","figS6_variable_importance.pdf"), p, width = 8.3, height = 4.3, device = pdf_device)
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
