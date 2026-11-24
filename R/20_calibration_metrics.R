# ======================================================================
# 20_calibration_metrics.R — Calibration and recalibrated agreement
# #
# Calibration slope, bias, raw and recalibrated R2 in the internal holdout and LC subset; change in squared Pearson correlation after adding each score (2,000 bootstrap resamples).
# #
# Inputs : results/holdout_predictions.csv, results/lc_predictions.csv
# Outputs: results/calibration_holdout_lc.csv, results/recalibration_holdout_lc.csv
# Paper  : Figure 5, Table S11
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table)

source("config/paths.R")

out <- file.path(BASE, "results")
cohorts <- c(holdout = "holdout_predictions.csv", LC = "lc_predictions.csv")
B <- 2000L
seed <- 2026L
raw_r2 <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)

calibrate <- function(cohort, file) {
  d <- fread(file.path(out, file))
  y <- d$y
  columns <- grep("__(ElasticNet|XGBoost)$", names(d), value = TRUE)
  stopifnot(length(columns) > 0, !anyNA(d[, c("y", columns), with = FALSE]))
  set.seed(seed)
  indices <- replicate(B, sample.int(nrow(d), replace = TRUE), simplify = FALSE)
  calibrated <- lapply(columns, function(column) {
    p <- d[[column]]
    fit <- lm(y ~ p)
    pearson_r <- cor(y, p)
    fit_r2 <- summary(fit)$r.squared
    stopifnot(isTRUE(all.equal(fit_r2, pearson_r^2, tolerance = 1e-12)))
    data.table(cohorte = cohort, modelo = column, R2 = round(raw_r2(y, p), 4),
               RMSE = round(sqrt(mean((y - p)^2)), 4), sesgo = round(mean(p - y), 4),
               pendiente = round(unname(coef(fit)[2]), 3), r = round(pearson_r, 4),
               R2_recal = round(fit_r2, 4))
  })
  incremental <- lapply(columns[!startsWith(columns, "S2_base__")], function(column) {
    parts <- strsplit(column, "__", fixed = TRUE)[[1]]
    base <- d[[paste0("S2_base__", parts[2])]]
    p <- d[[column]]
    estimate <- cor(y, p)^2 - cor(y, base)^2
    reps <- vapply(indices, function(i) cor(y[i], p[i])^2 - cor(y[i], base[i])^2,
                   numeric(1))
    ci <- quantile(reps, c(0.025, 0.975))
    data.table(cohorte = cohort, algoritmo = parts[2], especificacion = parts[1],
               dR2_recal = round(estimate, 4), IC_inf = round(unname(ci[1]), 4),
               IC_sup = round(unname(ci[2]), 4), cruza_cero = ci[1] <= 0 && ci[2] >= 0)
  })
  list(calibration = rbindlist(calibrated), incremental = rbindlist(incremental))
}

result <- Map(calibrate, names(cohorts), cohorts)
fwrite(rbindlist(lapply(result, "[[", "calibration")),
       file.path(out, "calibration_holdout_lc.csv"))
fwrite(rbindlist(lapply(result, "[[", "incremental")),
       file.path(out, "recalibration_holdout_lc.csv"))
