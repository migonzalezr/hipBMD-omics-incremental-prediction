# ======================================================================
# 22_table2_results.R — Table 2
#
# Association and incremental contribution of each score across cohorts (raw and recalibrated Delta R2).
#
# Inputs : results of scripts 09-13
# Outputs: results/table2_results.csv; TABLE2.md
# Paper  : Table 2
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls())
library(data.table)

source("config/paths.R")
P <- list(res  = file.path(BASE, "results"),
          docs = file.path(BASE, "results"))
SEED <- 2026; B <- 2000
rd <- function(f) { p <- file.path(P$res, f); if (!file.exists(p)) stop("falta ", f,
  " — correr antes el script que lo produce (16 para LC, 18 para Guangzhou)"); fread(file = p) }

SC <- data.table(
  clave = c("PRS", "PRS_qian", "MRS", "QMRS"),
  espec = c("S3_base_PRS", "S3q_base_PRSqian", "S4_base_MRS", "S5q_base_PRSqian_MRS"),
  etiqueta = c("PRS (own GWAS, LOS)", "PRS (external GWAS, Qian)", "MRS", "PRS (Qian) + MRS"))
COH <- c("Modeling discovery", "Internal holdout", "Validation 1: LC", "Validation 2: Guangzhou")

oof <- rd("discovery_predictions_oof.csv"); y <- oof$y; n <- length(y)
r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
set.seed(SEED); idx <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)
contr <- function(c1, c0) {
  p1 <- oof[[c1]]; p0 <- oof[[c0]]
  d <- vapply(idx, function(i) r2(y[i], p1[i]) - r2(y[i], p0[i]), numeric(1))
  ic <- quantile(d, c(0.025, 0.975))
  c(dR2 = round(r2(y, p1) - r2(y, p0), 4), lo = round(unname(ic[1]), 4), hi = round(unname(ic[2]), 4))
}
disc <- rbindlist(lapply(c("ElasticNet", "XGBoost"), function(a) rbindlist(lapply(seq_len(nrow(SC)), function(j) {
  v <- contr(paste0(SC$espec[j], "__", a), paste0("S2_base__", a))
  data.table(algoritmo = a, clave = SC$clave[j], dR2 = v[["dR2"]], lo = v[["lo"]], hi = v[["hi"]])
}))))

ref <- rd("discovery_conditional_contrast.csv")
chk <- merge(disc[clave %in% c("PRS_qian", "MRS")][, contraste := fifelse(clave == "MRS", "dR2(MRS)", "dR2(PRS Qian)")],
             ref, by = c("algoritmo", "contraste"), suffixes = c("", "_ref"))
stopifnot(nrow(chk) == 4)
dif <- chk[, max(abs(c(dR2 - dR2_ref, lo - IC95_inf, hi - IC95_sup)))]
if (dif > 1e-4) stop(sprintf("no reproduzco el bootstrap del script 15 (diferencia %.5f)", dif))
cat(sprintf("control discovery OK: reproduce discovery_conditional_contrast.csv (dif max %.1e)\n", dif))

hic <- rd("holdout_ci.csv")[, .(algoritmo, espec = especificacion, dR2 = dR2_vs_S2, lo = IC95_inf, hi = IC95_sup)]
lic <- rd("lc_ci.csv")[,      .(algoritmo, espec = especificacion, dR2, lo = IC95_inf, hi = IC95_sup)]
cic <- rd("guangzhou_ci.csv")[,   .(algoritmo, espec = especificacion, dR2, lo = IC95_inf, hi = IC95_sup)]
rlh <- rd("recalibration_holdout_lc.csv")
rh  <- rlh[cohorte == "holdout", .(algoritmo, espec = especificacion, dR2 = dR2_recal, lo = IC_inf, hi = IC_sup)]
rl  <- rlh[cohorte == "LC",      .(algoritmo, espec = especificacion, dR2 = dR2_recal, lo = IC_inf, hi = IC_sup)]
rc  <- rd("guangzhou_recalibration.csv")[, .(algoritmo, espec = especificacion, dR2 = dR2_recal, lo = IC95_inf, hi = IC95_sup)]

pega <- function(d, coh, tipo) merge(d, SC[, .(espec, clave)], by = "espec")[, `:=`(cohorte = coh, tipo = tipo)][]
dr <- rbindlist(list(
  disc[, `:=`(cohorte = COH[1], tipo = "raw")][],
  pega(hic, COH[2], "raw"), pega(lic, COH[3], "raw"), pega(cic, COH[4], "raw"),
  pega(rh,  COH[2], "recal"), pega(rl, COH[3], "recal"), pega(rc, COH[4], "recal")), fill = TRUE)
dr[, espec := NULL]
stopifnot(dr[, .N, by = .(algoritmo, cohorte, tipo)][, all(N == 4)])

ha <- rd("holdout_score_association.csv")
ass <- rbindlist(list(
  ha[conjunto == "discovery"][, cohorte := COH[1]],
  ha[conjunto == "holdout"][,   cohorte := COH[2]],
  rd("lc_score_association.csv")[, cohorte := COH[3]],
  rd("guangzhou_score_association.csv")[, cohorte := COH[4]]), fill = TRUE)
setnames(ass, "score", "clave")
NN <- c(1430, 351, 495, 492); names(NN) <- COH
stopifnot(all(ass[, n == NN[cohorte]]))

larga <- merge(dcast(dr, algoritmo + cohorte + clave ~ tipo, value.var = c("dR2", "lo", "hi")),
               ass[, .(cohorte, clave, n, beta, beta_lo = IC95_inf, beta_hi = IC95_sup, P)],
               by = c("cohorte", "clave"), all.x = TRUE)
larga[, `:=`(cohorte = factor(cohorte, levels = COH), clave = factor(clave, levels = SC$clave))]
setorder(larga, algoritmo, clave, cohorte)
fwrite(larga, file.path(P$res, "table2_results.csv"))

f4 <- function(x) ifelse(is.na(x), NA, sprintf("%+.4f", x))
ci <- function(e, l, h) ifelse(is.na(e), "&mdash;",
  paste0(f4(e), " (", f4(l), ", ", f4(h), ")", ifelse(l > 0 | h < 0, " *", "")))
fp <- function(p) ifelse(is.na(p), "&mdash;", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
nn <- setNames(c(1430, 351, 495, 492), COH)

tabla_md <- function(a) {
  d <- larga[algoritmo == a]
  d[, `:=`(Score = SC$etiqueta[match(clave, SC$clave)],
           beta_txt = ci(beta, beta_lo, beta_hi), P_txt = fp(P),
           raw_txt = ci(dR2_raw, lo_raw, hi_raw),
           recal_txt = ifelse(cohorte == COH[1], "n/a", ci(dR2_recal, lo_recal, hi_recal)))]
  d[clave == "QMRS", `:=`(beta_txt = "&mdash;", P_txt = "&mdash;")]
  out <- c("| Score | Cohort | n | &beta; per SD (95% CI) | P | &Delta;R&sup2; raw (95% CI) | &Delta;R&sup2; recalibrated (95% CI) |",
           "|:--|:--|--:|:--|--:|:--|:--|")
  for (i in seq_len(nrow(d))) out <- c(out, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
      if (i == 1 || d$clave[i] != d$clave[i-1]) d$Score[i] else "", d$cohorte[i],
      format(nn[as.character(d$cohorte[i])], big.mark = ","), d$beta_txt[i], d$P_txt[i], d$raw_txt[i], d$recal_txt[i]))
  out
}
notas <- c(
  "",
  "&beta;: change in total hip BMD (g/cm&sup2;) per SD of the score, from a linear model adjusted for the clinical covariates and the first two genetic principal components (Guangzhou: the six covariates of the reduced model; no genetic principal components exist there and sex is constant). Scores were standardized with discovery parameters in the holdout and within cohort in LC and Guangzhou. No single &beta; exists for the combined row.",
  "&Delta;R&sup2;: gain over the base model (clinical covariates + genetic principal components) with the score added; 95% CIs from 2,000 bootstrap resamples of participants. * CI excludes zero.",
  "Modeling discovery: out-of-fold predictions from one 5-fold split (seed 2026). PRS thresholds were selected in this same set, so the PRS rows are optimistic here; the repeated 5 &times; 5 cross-validation used for selection gave 0.0043 (LOS) and 0.0042 (Qian), with P < 0.007 against 150 phenotype permutations for both (Figure 2). Recalibration does not apply to out-of-fold predictions in the training set.",
  "Recalibrated &Delta;R&sup2;: computed after a linear recalibration fitted in the target cohort itself; it separates ranking from calibration and is an upper bound, not a deployable estimate (Figure 4). In a miscalibrated cohort, any added predictor that shifts predictions toward the cohort mean raises raw R&sup2; without adding information; this is why raw and recalibrated estimates diverge in Guangzhou, where the base model is offset by &minus;0.074 g/cm&sup2;.",
  "Guangzhou: base model reduced to the six predictors available in that cohort (age, weight, height, BMI, alcohol, exercise), trained in the modeling discovery set.",
  "MRS: the pre-specified MRS (balances at FDR < 0.05) selected zero balances in the discovery set; the MRS shown is the secondary, threshold-based score (P < 0.05).")

md <- c("# Table 2. Association and predictive contribution of each score across the four cohorts",
        "", "Elastic Net (primary algorithm).", "", tabla_md("ElasticNet"), notas,
        "", "", "# Supplementary Table S_. Same as Table 2, XGBoost", "", tabla_md("XGBoost"),
        "", "Footnotes as in Table 2. &beta; comes from a linear model and does not depend on the algorithm.")

writeLines(md, file.path(P$docs, "TABLE2.md"))

cat("\n=== Tabla 2, Elastic Net ===\n")
print(larga[algoritmo == "ElasticNet", .(clave, cohorte, beta, P, dR2_raw, lo_raw, hi_raw, dR2_recal, lo_recal, hi_recal)],
      row.names = FALSE)
cat("\nguardado en results/table2_results.csv y results/TABLE2.md\n")
cat("===== FIN OK =====\n")
