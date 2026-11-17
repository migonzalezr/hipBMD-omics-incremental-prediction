# ======================================================================
# fig4_calibration.R — Figure 4
#
# Calibration of the base model in the internal holdout, the LC subset and the Guangzhou cohort.
#
# Inputs : results/*_predictions.csv; calibration_holdout_lc.csv; guangzhou_calibration.csv
# Outputs: figures/png|pdf/fig4_calibration.*
# Paper  : Figure 4
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls())
library(data.table); library(ggplot2)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")
dir.create(file.path(FIG,"png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIG,"pdf"), showWarnings = FALSE, recursive = TRUE)

MOD <- "S2_base__ElasticNet"
COH <- data.table(
  archivo = c("holdout_predictions.csv", "lc_predictions.csv", "guangzhou_predictions.csv"),
  clave   = c("holdout", "LC", "Guangzhou"),
  etiqueta= c("Internal holdout", "LC subset", "Guangzhou cohort"),
  col     = c("#2b6cb0", "#d98c00", "#117a65"))
COH[, etiqueta := factor(etiqueta, levels = etiqueta)]

dd <- rbindlist(lapply(seq_len(nrow(COH)), function(i) {
  d <- fread(file = file.path(RES, COH$archivo[i]))
  stopifnot(MOD %in% names(d))
  data.table(etiqueta = COH$etiqueta[i], obs = d$y, pred = d[[MOD]])
}))
stopifnot(dd[, !anyNA(obs) && !anyNA(pred)])

mt <- dd[, {
  f <- coef(lm(obs ~ pred))
  .(n = .N,
    R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    sesgo = mean(pred - obs),
    pendiente = unname(f[2]), b0 = unname(f[1]),
    r = cor(obs, pred))
}, by = etiqueta]
mt[, R2_recal := r^2]

rg <- dd[, .(x0 = min(pred), x1 = max(pred)), by = etiqueta]

ref <- rbind(
  fread(file = file.path(RES, "calibration_holdout_lc.csv"))[modelo == MOD,
        .(clave = cohorte, R2, sesgo, pendiente, r, R2_recal)],
  fread(file = file.path(RES, "guangzhou_calibration.csv"))[modelo == MOD,
        .(clave = "Guangzhou", R2, sesgo, pendiente, r, R2_recal)])
mt <- merge(mt, COH[, .(etiqueta, clave, col)], by = "etiqueta")
chk <- merge(mt, ref, by = "clave", suffixes = c("", "_ref"))
stopifnot(nrow(chk) == 3)
dif <- chk[, max(abs(c(R2 - R2_ref, sesgo - sesgo_ref, pendiente - pendiente_ref,
                       r - r_ref, R2_recal - R2_recal_ref)))]
if (dif > 1e-3) stop(sprintf("las metricas recalculadas difieren de las guardadas en %.4f", dif))
cat(sprintf("control OK: diferencia maxima con las tablas guardadas = %.2e\n", dif))
print(mt[, .(etiqueta, n, R2 = round(R2,4), sesgo = round(sesgo,4),
             pendiente = round(pendiente,3), r = round(r,3), R2_recal = round(R2_recal,4))],
      row.names = FALSE)

LIM <- c(0.55, 1.55)
lab <- mt[, .(etiqueta, txt = sprintf(
  "n = %d\nr = %.2f\nbias = %+.3f g/cm\u00b2\nslope = %.2f\nR\u00b2 = %.2f\nrecalibrated R\u00b2 = %.2f",
  n, r, sesgo, pendiente, R2, R2_recal))]

p <- ggplot(dd, aes(pred, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.4, linetype = "22") +
  geom_point(aes(colour = etiqueta), alpha = 0.40, size = 0.9, stroke = 0) +
  geom_segment(data = merge(mt, rg, by = "etiqueta"),
               aes(x = x0, xend = x1, y = b0 + pendiente * x0, yend = b0 + pendiente * x1),
               colour = "grey10", linewidth = 0.75, inherit.aes = FALSE) +
  geom_text(data = lab, aes(x = LIM[1] + 0.03, y = LIM[2] - 0.03, label = txt),
            hjust = 0, vjust = 1, size = 2.3, lineheight = 1.0, colour = "grey20",
            inherit.aes = FALSE) +
  facet_wrap(~ etiqueta, nrow = 1) +
  scale_colour_manual(values = setNames(COH$col, levels(COH$etiqueta)), guide = "none") +

  coord_fixed(xlim = LIM, ylim = LIM, expand = FALSE) +
  scale_x_continuous(breaks = seq(0.6, 1.5, 0.3)) +
  scale_y_continuous(breaks = seq(0.6, 1.5, 0.3)) +
  labs(title = "Transfer fails on calibration, not on ranking",
       subtitle = paste0(
         "Elastic Net base model trained in the modeling discovery set (BGI). Dashed: identity; ",
         "solid: calibration (observed ~ predicted), over the range of predictions.\n",
         "Guangzhou uses the base model reduced to the six predictors available there ",
         "(age, weight, height, BMI, alcohol, exercise).\n",
         "Bias = mean(predicted \u2212 observed). Recalibrated R\u00b2 = r\u00b2, the upper bound ",
         "reached by a linear recalibration fitted in the target cohort itself."),
       x = expression("Predicted total hip BMD (g/cm"^2*")"),
       y = expression("Observed total hip BMD (g/cm"^2*")")) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9.5, hjust = 0),
        plot.subtitle = element_text(size = 6.4, colour = "grey38", lineheight = 1.25),
        strip.background = element_blank(),
        strip.text = element_text(size = 8.5, face = "bold", hjust = 0),
        panel.spacing = unit(1.0, "lines"))

ggsave(file.path(FIG,"png","fig4_calibration.png"), p, width = 8.4, height = 3.9, dpi = 300)
ggsave(file.path(FIG,"pdf","fig4_calibration.pdf"), p, width = 8.4, height = 3.9, device = cairo_pdf)
cat("guardado en figures/png/ y figures/pdf/\n")
cat("===== FIN OK =====\n")
