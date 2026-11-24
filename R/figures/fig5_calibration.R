# ======================================================================
# fig5_calibration.R — Figure 5
# #
# Calibration and linear agreement of the clinical base model in the internal holdout, LC subset and Guangzhou cohort.
# #
# Inputs : results/*_predictions.csv, calibration files
# Outputs: figures/png, figures/pdf
# Paper  : Figure 5
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table)
library(ggplot2)
library(patchwork)

source("config/paths.R")
RES  <- file.path(BASE, "results")
FIG  <- Sys.getenv("FIG4_OUTPUT_DIR", unset = file.path(BASE, "figures"))
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
  fit <- lm(obs ~ pred)
  f <- coef(fit)
  .(n = .N,
    R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    sesgo = mean(pred - obs),
    pendiente = unname(f[2]), b0 = unname(f[1]),
    r = cor(obs, pred, method = "pearson"),
    R2_recal = summary(fit)$r.squared)
}, by = etiqueta]
stopifnot(isTRUE(all.equal(mt$R2_recal, mt$r^2, tolerance = 1e-12)))
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

lim <- c(0.55, 1.55)
stopifnot(all(dd$obs >= lim[1] & dd$obs <= lim[2]),
          all(dd$pred >= lim[1] & dd$pred <= lim[2]))
keys <- c("Identity (observed = predicted)", "Calibration fit (observed on predicted)")
line_types <- setNames(c("22", "solid"), keys)
line_colours <- setNames(c("grey45", "grey10"), keys)
identity_line <- data.table(line = keys[1], slope = 1, intercept = 0)
fit_segments <- merge(mt, rg, by = "etiqueta")[, line := keys[2]]
p <- ggplot(dd, aes(pred, obs)) +
  geom_abline(data = identity_line,
              aes(slope = slope, intercept = intercept, linetype = line, colour = line),
              linewidth = 0.58, key_glyph = draw_key_path) +
  geom_point(aes(fill = etiqueta), shape = 21, alpha = 0.50, size = 1.05, stroke = 0) +
  geom_segment(data = fit_segments,
               aes(x = x0, xend = x1, y = b0 + pendiente * x0, yend = b0 + pendiente * x1,
                   linetype = line, colour = line),
               linewidth = 0.75, inherit.aes = FALSE) +
  facet_wrap(~ etiqueta, nrow = 1) +
  scale_fill_manual(values = setNames(COH$col, levels(COH$etiqueta)), guide = "none") +
  scale_colour_manual(values = line_colours, breaks = keys, name = NULL) +
  scale_linetype_manual(values = line_types, breaks = keys, name = NULL) +
  coord_fixed(xlim = lim, ylim = lim, expand = FALSE) +
  scale_x_continuous(breaks = seq(0.6, 1.5, 0.3)) +
  scale_y_continuous(breaks = seq(0.6, 1.5, 0.3)) +
  labs(x = expression("Original predicted total hip BMD (g/cm"^2*")"),
       y = expression("Observed total hip BMD (g/cm"^2*")")) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold", hjust = 0),
        panel.spacing = unit(1.1, "lines"),
        legend.position = "bottom",
        legend.text = element_text(size = 8.5),
        legend.key.width = unit(24, "pt"))

metrics <- mt[, .(etiqueta,
  "Participants, n" = format(n),
  "Mean bias (g/cm²)" = sprintf("%+.3f", sesgo),
  "Calibration slope" = sprintf("%.2f", pendiente),
  "Raw model R²" = sprintf("%.3f", R2),
  "Pearson r" = sprintf("%.3f", r),
  "Recalibrated R²" = sprintf("%.3f", R2_recal)
)]
values <- melt(metrics, id.vars = "etiqueta", variable.name = "metric")
values[, metric := factor(metric, levels = rev(names(metrics)[-1]))]
values[, band := as.integer(metric) %% 2L]
summary_panel <- ggplot(values, aes(x = 1, y = metric)) +
  geom_tile(aes(fill = factor(band)), width = 1, height = 1) +
  geom_text(aes(label = value), size = 3.1) +
  scale_fill_manual(values = c("white", "#EFF3F6"), guide = "none") +
  scale_x_continuous(limits = c(0.5, 1.5), expand = expansion()) +
  scale_y_discrete(expand = expansion()) +
  facet_wrap(~ etiqueta, nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(axis.text.y = element_text(colour = "grey15", hjust = 1, margin = margin(r = 6)),
        strip.text = element_blank(),
        panel.spacing = unit(1.1, "lines"),
        plot.margin = margin(0, 5.5, 5.5, 5.5))
figure <- p / summary_panel + plot_layout(heights = c(3, 1.2))

ggsave(file.path(FIG, "png", "fig4_calibration.png"), figure,
       width = 8.4, height = 4.55, dpi = 300, bg = "white")
pdf_device <- if (Sys.info()[["sysname"]] == "Darwin") {
  function(filename, ...) quartz(type = "pdf", file = filename, ...)
} else {
  cairo_pdf
}
ggsave(file.path(FIG, "pdf", "fig4_calibration.pdf"), figure,
       width = 8.4, height = 4.55, device = pdf_device, bg = "white")
