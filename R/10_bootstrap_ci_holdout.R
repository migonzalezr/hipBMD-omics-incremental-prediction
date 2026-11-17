# ======================================================================
# 10_bootstrap_ci_holdout.R — Bootstrap confidence intervals
#
# Percentile 95% CIs for R2 and Delta R2 from 2,000 bootstrap resamples with the fitted models held fixed.
#
# Inputs : results/holdout_predictions.csv
# Outputs: results/holdout_ci.csv; holdout_ci_conditional.csv
# Paper  : Table 2
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table)

source("config/paths.R")
OUT  <- file.path(BASE,"results")
B <- 2000; SEED <- 2026

PRED <- fread(file = file.path(OUT,"holdout_predictions.csv"))
y <- PRED$y; n <- length(y)
cols <- setdiff(names(PRED), c("IID","y"))
cat(sprintf("holdout: %d personas | %d columnas de prediccion\n", n, length(cols)))

r2 <- function(obs,pred) 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)

set.seed(SEED)
idx <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)

R2B <- sapply(cols, function(cn) {
  p <- PRED[[cn]]
  vapply(idx, function(i) r2(y[i], p[i]), numeric(1))
})

partes <- tstrsplit(cols, "__", fixed = TRUE)
info <- data.table(col = cols, especificacion = partes[[1]], algoritmo = partes[[2]])

res <- rbindlist(lapply(unique(info$algoritmo), function(a) {
  cn_base <- info[algoritmo==a & especificacion=="S2_base", col]
  rbindlist(lapply(info[algoritmo==a, col], function(cn) {
    d <- R2B[, cn] - R2B[, cn_base]
    ic <- quantile(d, c(0.025, 0.975))
    data.table(algoritmo = a,
               especificacion = info[col==cn, especificacion],
               R2_holdout = round(r2(y, PRED[[cn]]), 4),
               dR2_vs_S2  = round(r2(y, PRED[[cn]]) - r2(y, PRED[[cn_base]]), 4),
               IC95_inf   = round(ic[1], 4),
               IC95_sup   = round(ic[2], 4),
               cruza_cero = (ic[1] <= 0 & ic[2] >= 0))
  }))
}))
setorder(res, algoritmo, especificacion)

cat("\n=== HOLDOUT CON INTERVALOS DE CONFIANZA (bootstrap, B = 2000) ===\n")
print(res, row.names = FALSE)
fwrite(res, file.path(OUT,"holdout_ci.csv"))

CONTR <- list(
  "dR2(MRS | PRS LOS)"  = c("S5_base_PRS_MRS",      "S3_base_PRS"),
  "dR2(MRS | PRS Qian)" = c("S5q_base_PRSqian_MRS", "S3q_base_PRSqian"),
  "dR2(PRS Qian | MRS)" = c("S5q_base_PRSqian_MRS", "S4_base_MRS")
)
cond <- rbindlist(lapply(unique(info$algoritmo), function(a) {
  rbindlist(lapply(names(CONTR), function(nm) {
    cn1 <- info[algoritmo==a & especificacion==CONTR[[nm]][1], col]
    cn0 <- info[algoritmo==a & especificacion==CONTR[[nm]][2], col]
    if (!length(cn1) || !length(cn0)) return(NULL)
    d  <- R2B[, cn1] - R2B[, cn0]
    ic <- quantile(d, c(0.025, 0.975))
    data.table(algoritmo=a, contraste=nm,
               dR2 = round(r2(y, PRED[[cn1]]) - r2(y, PRED[[cn0]]), 4),
               IC95_inf = round(ic[1],4), IC95_sup = round(ic[2],4),
               cruza_cero = (ic[1] <= 0 & ic[2] >= 0))
  }))
}))
cat("
=== CONTRASTES CONDICIONALES (responden a R1.1, complementariedad) ===
")
print(cond, row.names = FALSE)
fwrite(cond, file.path(OUT,"holdout_ci_conditional.csv"))

anchos <- res[especificacion!="S2_base", IC95_sup - IC95_inf]
cat(sprintf("\nancho tipico del IC95 del dR2: %.4f\n", median(anchos)))
cat(sprintf("para comparar, el dR2 del PRS en discovery fue 0.0043\n"))
cat(sprintf("-> con n = %d, el holdout NO puede resolver efectos de ese tamano\n", n))
cat("\nguardado en results/holdout_ci.csv\n")
cat("===== FIN OK =====\n")
