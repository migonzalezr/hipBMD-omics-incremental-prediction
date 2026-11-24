# ======================================================================
# fig3_prediction_performance.R — Figure 3
# #
# Raw R2 of the clinical base model and score-augmented models in the four cohorts (Elastic Net and XGBoost), with bootstrap CIs.
# #
# Inputs : results/*_predictions.csv, results/table2_results.csv
# Outputs: figures/png, figures/pdf
# Paper  : Figure 3
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table)
library(ggplot2)
library(patchwork)
library(grid)

source("config/paths.R")

res <- file.path(BASE, "results")
cohorts <- c("Modeling discovery", "Internal holdout", "LC subset", "Guangzhou cohort")
table_cohorts <- c(cohorts[1:2], "LC subset", "Guangzhou cohort")
files <- c("discovery_predictions_oof.csv", "holdout_predictions.csv",
           "lc_predictions.csv", "guangzhou_predictions.csv")
algorithms <- c("ElasticNet", "XGBoost")
models <- c("S2_base", "S3_base_PRS", "S3q_base_PRSqian", "S4_base_MRS", "S5q_base_PRSqian_MRS")
keys <- c("base", "PRS", "PRS_qian", "MRS", "QMRS")
axis_labels <- c("Clinical base", "+ LOS PRS", "+ UK Biobank\nPRS", "+ MRS\nonly",
                 "+ UK Biobank\nPRS and MRS")
raw_r2 <- function(y, pred) 1 - sum((y - pred)^2) / sum((y - mean(y))^2)

d <- rbindlist(lapply(seq_along(files), function(i) {
  pred <- fread(file.path(res, files[i]))
  columns <- as.vector(outer(models, algorithms, paste, sep = "__"))
  fitted <- as.matrix(pred[, ..columns])
  stopifnot(!anyDuplicated(pred$IID), all(is.finite(pred$y)), all(is.finite(fitted)))
  evaluate <- function(j) apply(fitted[j, , drop = FALSE], 2, function(x) raw_r2(pred$y[j], x))
  set.seed(2026)
  bootstrap <- replicate(2000, evaluate(sample.int(nrow(pred), replace = TRUE)))
  ci <- apply(bootstrap, 1, quantile, probs = c(0.025, 0.975))
  data.table(cohort = cohorts[i], table_cohort = table_cohorts[i], n = nrow(pred),
             algorithm = rep(algorithms, each = length(models)), position = seq_along(models),
             score_key = keys, r2 = unname(evaluate(seq_len(nrow(pred)))),
             lower = ci[1, ], upper = ci[2, ])
}))
stopifnot(nrow(d) == 40L, all(d$lower <= d$r2), all(d$r2 <= d$upper))
d[, delta := r2 - r2[position == 1], by = .(cohort, algorithm)]
reference <- fread(file.path(res, "table2_results.csv"))[
  , .(table_cohort = cohorte, algorithm = algoritmo, score_key = clave, dR2_raw, table_n = n)]
check <- merge(d[position > 1], reference, by = c("table_cohort", "algorithm", "score_key"))
stopifnot(nrow(check) == 32L, all(check$n == check$table_n),
          max(abs(check$delta - check$dR2_raw)) < 0.00016)
print(d[, .(min = min(lower), max = max(upper), min_r2 = min(r2), max_r2 = max(r2)),
        by = .(algorithm, cohort)])
print(d[, as.list(setNames(round(r2, 4), keys)), by = .(algorithm, cohort)])

colors <- setNames(c("#1b1b1b", "#65b3e3", "#a3650a", "#117863"), cohorts)
text_colors <- setNames(c("#1b1b1b", "#2175a5", "#8a5208", "#116a59"), cohorts)
shapes <- setNames(c(16, 17, 15, 18), cohorts)
d[, cohort := factor(cohort, levels = cohorts)]
d[, x := position + fifelse(cohort == cohorts[1], -0.085,
                             fifelse(cohort == cohorts[2], 0.085, 0))]
d[, label_y := r2 + fifelse(cohort == cohorts[1], 0.015,
                            fifelse(cohort == cohorts[2], -0.017,
                                    fifelse(cohort == cohorts[3], 0.015, 0.024)))]
d[, label_x := x + fifelse(cohort == cohorts[1], -0.052, 0.065)]
d[, label_hjust := as.integer(cohort == cohorts[1])]
d[position == 5 & cohort %in% cohorts[3:4],
  `:=`(label_x = x - 0.06, label_hjust = 1L)]
d[position == 5 & cohort %in% cohorts[3:4], label_y := r2 + fifelse(cohort == cohorts[3], 0.023, 0.038)]
d[, label := gsub("-", "\u2212", sprintf("%.4f", r2), fixed = TRUE)]
right_labels <- d[position == 5]
right_labels[, name_y := r2 + c(0.022, 0.047, 0.031, 0.061)[as.integer(cohort)]]
left_edge <- 0.45

band <- function(algo, groups, limits, ticks, x_axis = FALSE) {
  dt <- d[algorithm == algo & cohort %in% groups]
  end <- right_labels[algorithm == algo & cohort %in% groups]
  stopifnot(all(dt$lower > limits[1]), all(dt$upper < limits[2]))
  p <- ggplot(dt, aes(x, r2, colour = cohort, shape = cohort, group = cohort)) +
    geom_vline(xintercept = 1:5, colour = "#e9edf1", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 0.52, width = 0.095, alpha = 0.58) +
    geom_line(linewidth = 0.83) +
    geom_point(size = 2.35) +
    geom_text(aes(x = label_x, y = label_y, label = label, hjust = label_hjust),
              colour = text_colors[as.character(dt$cohort)], size = 2.8, show.legend = FALSE) +
    geom_segment(data = end, aes(x = 5.31, xend = 5.52, y = r2, yend = name_y),
                 linewidth = 0.48, alpha = 0.85) +
    geom_text(data = end, aes(x = 5.57, y = name_y, label = cohort),
              colour = text_colors[as.character(end$cohort)], hjust = 0, size = 3.2, fontface = "bold") +
    scale_x_continuous(breaks = 1:5, labels = axis_labels, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = ticks, expand = expansion(mult = 0),
                       labels = function(x) gsub("-", "\u2212", sprintf("%.2f", x), fixed = TRUE)) +
    scale_colour_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    coord_cartesian(xlim = c(left_edge, 5.29), ylim = limits, clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(colour = "#e7ebef", linewidth = 0.41),
          axis.text = element_text(colour = "#263746"),
          axis.text.x = element_text(size = 9, margin = margin(t = 5)),
          axis.text.y = element_text(size = 9),
          axis.line = element_line(colour = "#303c48", linewidth = 0.43),
          plot.margin = margin(6, 114, 7, 8))
  if (!x_axis) {
    p <- p + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                   axis.text.x = element_blank()) +
      annotate("segment", x = left_edge - 0.033, xend = left_edge + 0.033,
               y = limits[1] - c(0.011, 0.021), yend = limits[1] - c(0.003, 0.013),
               linewidth = 0.50, colour = "#303c48")
  }
  p
}

heading <- function(letter, name) {
  wrap_elements(full = grobTree(
    textGrob(letter, x = 0.022, hjust = 0, gp = gpar(fontsize = 16, fontface = "bold")),
    textGrob(name, x = 0.106, hjust = 0, gp = gpar(fontsize = 13, fontface = "bold"))
  ))
}
panel <- function(algo, letter, title) {
  top <- band(algo, cohorts[1:2], c(0.28, 0.475), seq(0.28, 0.46, by = 0.03))
  mid <- band(algo, cohorts[3], c(0.105, 0.275), seq(0.12, 0.27, by = 0.03))
  bottom <- band(algo, cohorts[4], c(-0.405, -0.015), seq(-0.40, -0.05, by = 0.05),
                 x_axis = TRUE)
  heading(letter, title) / top / mid / bottom +
    plot_layout(heights = c(25, 155, 128, 179))
}
axis <- wrap_elements(full = textGrob("Raw R²", rot = 90, x = 0.64,
                                     gp = gpar(fontsize = 11, col = "#263746")))
a <- panel(algorithms[1], "A", "Elastic Net")
b <- panel(algorithms[2], "B", "XGBoost")
plots <- (axis | a | b) + plot_layout(widths = c(0.036, 1, 1))
note <- wrap_elements(full = textGrob("Whiskers show 95% bootstrap confidence intervals.",
                                      x = 0.055, hjust = 0,
                                      gp = gpar(fontsize = 9.5, col = "#465767")))
figure <- (plots / note) + plot_layout(heights = c(1, 0.047))
name <- "fig3_prediction_performance"
ggsave(file.path(BASE, "figures", "png", paste0(name, ".png")), figure,
       width = 13.5, height = 8, dpi = 320, bg = "white")
pdf_device <- function(filename, ...) quartz(type = "pdf", file = filename, ...)
ggsave(file.path(BASE, "figures", "pdf", paste0(name, ".pdf")), figure,
       width = 13.5, height = 8, device = pdf_device, bg = "white")
