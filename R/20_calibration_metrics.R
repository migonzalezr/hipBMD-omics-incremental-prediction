# ======================================================================
# 20_calibration_metrics.R — Calibration of the internal holdout and the LC subset
#
# R2, RMSE, bias = mean(predicted - observed), calibration slope
# (observed ~ predicted), correlation r and recalibrated R2 = r^2,
# for every model specification.
#
# Inputs : results/holdout_predictions.csv; results/lc_predictions.csv
# Outputs: results/calibration_holdout_lc.csv
# Paper  : Figure 4; Results (transfer across platforms and populations)
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

library(data.table)
source("config/paths.R")
RES <- file.path(BASE, "results")

calib <- function(y, p) {
  r <- cor(y, p)
  c(R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2),
    RMSE = sqrt(mean((y - p)^2)),
    sesgo = mean(p - y),
    pendiente = unname(coef(lm(y ~ p))[2]),
    r = r, R2_recal = r^2)
}

files <- c(LC = "lc_predictions.csv", holdout = "holdout_predictions.csv")
out <- rbindlist(lapply(files, function(f) {
  d <- fread(file.path(RES, f))
  rbindlist(lapply(setdiff(names(d), c("IID", "y")), function(m)
    as.data.table(as.list(c(modelo = m, calib(d$y, d[[m]]))))))
}), idcol = "cohorte")
num <- c("R2", "RMSE", "sesgo", "pendiente", "r", "R2_recal")
out[, (num) := lapply(.SD, as.numeric), .SDcols = num]
out[, (setdiff(num, "pendiente")) := lapply(.SD, round, 4), .SDcols = setdiff(num, "pendiente")]
out[, pendiente := round(pendiente, 3)]
fwrite(out, file.path(RES, "calibration_holdout_lc.csv"))
