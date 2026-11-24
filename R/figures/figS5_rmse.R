# ======================================================================
# figS5_rmse.R — Figure S5
#
# Root mean squared error (RMSE, g/cm2) of the clinical base model and score-augmented models in the
# four evaluation sets (Elastic Net and XGBoost), with 2,000-resample participant-level bootstrap
# percentile CIs (predictions held fixed). Discovery uses out-of-fold predictions.
#
# Inputs : results/discovery_predictions_oof.csv, results/holdout_predictions.csv,
#          results/lc_predictions.csv, results/guangzhou_predictions.csv
# Outputs: figures/png/figS5_rmse.png, figures/pdf/figS5_rmse.pdf, results/figS5_rmse_values.csv
# Paper  : Figure S5
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table)
library(ggplot2)
library(patchwork)
library(grid)

source("config/paths.R")

res <- file.path(BASE, "results")
cohorts <- c("Modeling discovery", "Internal holdout", "LC subset", "Guangzhou cohort")
files <- c("discovery_predictions_oof.csv", "holdout_predictions.csv",
           "lc_predictions.csv", "guangzhou_predictions.csv")
algorithms <- c("ElasticNet", "XGBoost")
models <- c("S2_base", "S3_base_PRS", "S3q_base_PRSqian", "S4_base_MRS", "S5q_base_PRSqian_MRS")
axis_labels <- c("Clinical\nbase model", "+ LOS PRS", "+ UK Biobank\nPRS", "+ MRS\nonly",
                 "+ UK Biobank\nPRS and MRS")
rmse <- function(y, pred) sqrt(mean((y - pred)^2))

d <- rbindlist(lapply(seq_along(files), function(i) {
  pred <- fread(file.path(res, files[i]))
  columns <- as.vector(outer(models, algorithms, paste, sep = "__"))
  fitted <- as.matrix(pred[, ..columns])
  stopifnot(!anyDuplicated(pred$IID), all(is.finite(pred$y)), all(is.finite(fitted)))
  evaluate <- function(j) apply(fitted[j, , drop = FALSE], 2, function(x) rmse(pred$y[j], x))
  set.seed(2026)
  bootstrap <- replicate(2000, evaluate(sample.int(nrow(pred), replace = TRUE)))
  ci <- apply(bootstrap, 1, quantile, probs = c(0.025, 0.975))
  data.table(cohort = cohorts[i], n = nrow(pred),
             algorithm = rep(algorithms, each = length(models)), model = rep(models, 2),
             position = rep(seq_along(models), 2),
             rmse = unname(evaluate(seq_len(nrow(pred)))), lower = ci[1, ], upper = ci[2, ])
}))
stopifnot(nrow(d) == 40L, all(d$lower <= d$rmse), all(d$rmse <= d$upper))
fwrite(d, file.path(res, "figS5_rmse_values.csv"))
print(d[, as.list(setNames(round(rmse, 4), models)), by = .(algorithm, cohort, n)])

colors <- setNames(c("#1b1b1b", "#65b3e3", "#a3650a", "#117863"), cohorts)
text_colors <- setNames(c("#1b1b1b", "#2175a5", "#8a5208", "#116a59"), cohorts)
shapes <- setNames(c(16, 17, 15, 18), cohorts)
d[, cohort := factor(cohort, levels = cohorts)]
d[, x := position + c(-0.12, -0.04, 0.04, 0.12)[as.integer(cohort)]]
limits <- c(0.098, 0.1475)

panel <- function(algo, letter, title) {
  dt <- d[algorithm == algo]
  stopifnot(all(dt$lower > limits[1]), all(dt$upper < limits[2]))
  # right-hand cohort names, ordered by the last RMSE and spaced to avoid overlap
  end <- dt[position == 5][order(rmse)]
  end[, name_y := rmse]
  gap <- 0.0055
  for (k in 2:nrow(end)) end$name_y[k] <- max(end$name_y[k], end$name_y[k - 1] + gap)
  ggplot(dt, aes(x, rmse, colour = cohort, shape = cohort, group = cohort)) +
    geom_vline(xintercept = 1:5, colour = "#e9edf1", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 0.5, width = 0.06, alpha = 0.6) +
    geom_line(linewidth = 0.83) +
    geom_point(size = 2.3) +
    geom_segment(data = end, aes(x = x + 0.05, xend = 5.42, y = rmse, yend = rmse),
                 linewidth = 0.45, alpha = 0.85) +
    geom_segment(data = end, aes(x = 5.42, xend = 5.55, y = rmse, yend = name_y),
                 linewidth = 0.45, alpha = 0.85) +
    geom_text(data = end, aes(x = 5.6, y = name_y, label = cohort),
              colour = text_colors[as.character(end$cohort)], hjust = 0, size = 3.1, fontface = "bold") +
    scale_x_continuous(breaks = 1:5, labels = axis_labels, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq(0.100, 0.145, by = 0.005), expand = expansion(mult = 0),
                       labels = function(x) sprintf("%.3f", x)) +
    scale_colour_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    coord_cartesian(xlim = c(0.45, 5.35), ylim = limits, clip = "off") +
    labs(x = NULL, y = NULL, title = title, tag = letter) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 13),
          plot.tag = element_text(face = "bold", size = 16),
          panel.grid.major.y = element_line(colour = "#e7ebef", linewidth = 0.41),
          axis.text = element_text(colour = "#263746"),
          axis.text.x = element_text(size = 9, margin = margin(t = 5)),
          axis.text.y = element_text(size = 9),
          axis.line = element_line(colour = "#303c48", linewidth = 0.43),
          plot.margin = margin(6, 125, 7, 8))
}
axis <- wrap_elements(full = textGrob("RMSE (g/cm\u00b2)", rot = 90, x = 0.6,
                                     gp = gpar(fontsize = 11, col = "#263746")))
figure <- (axis | panel(algorithms[1], "a", "Elastic Net") | panel(algorithms[2], "b", "XGBoost")) +
  plot_layout(widths = c(0.036, 1, 1))
name <- "figS5_rmse"
ggsave(file.path(BASE, "figures", "png", paste0(name, ".png")), figure,
       width = 14.4, height = 6.45, dpi = 300, bg = "white")
ggsave(file.path(BASE, "figures", "pdf", paste0(name, ".pdf")), figure,
       width = 14.4, height = 6.45, device = pdf, bg = "white")
