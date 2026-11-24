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
COH <- c("Modeling discovery", "Internal holdout", "LC subset", "Guangzhou cohort")

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
larga[, n := unname(NN[as.character(cohorte)])]
setorder(larga, algoritmo, clave, cohorte)
fwrite(larga, file.path(P$res, "table2_results.csv"))

f4 <- function(x) ifelse(is.na(x), NA, sprintf("%+.4f", x))
ci <- function(e, l, h) ifelse(is.na(e), "&mdash;",
  paste0(f4(e), " (", f4(l), ", ", f4(h), ")", ifelse(l > 0 | h < 0, " *", "")))
fp <- function(p) ifelse(is.na(p), "&mdash;", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))

tabla_md <- function(a) {
  d <- larga[algoritmo == a]
  d[, `:=`(Score = SC$etiqueta[match(clave, SC$clave)],
           beta_txt = ci(beta, beta_lo, beta_hi), P_txt = fp(P),
           raw_txt = ci(dR2_raw, lo_raw, hi_raw),
           recal_txt = ifelse(cohorte == COH[1], "n/a", ci(dR2_recal, lo_recal, hi_recal)))]
  d[clave == "QMRS", `:=`(beta_txt = "&mdash;", P_txt = "&mdash;")]
  out <- c("| Score | Cohort | &beta; per SD (95% CI)<sup>a</sup> | P<sup>a</sup> | &Delta;R&sup2; raw (95% CI)<sup>b</sup> | &Delta;R&sup2; recalibrated (95% CI)<sup>c</sup> |",
           "|:--|:--|:--|--:|:--|:--|")
  for (i in seq_len(nrow(d))) out <- c(out, sprintf("| %s | %s | %s | %s | %s | %s |",
      if (i == 1 || d$clave[i] != d$clave[i-1]) d$Score[i] else "", d$cohorte[i],
      d$beta_txt[i], d$P_txt[i], d$raw_txt[i], d$recal_txt[i]))
  out
}
notas <- c(
  "",
  "Sample sizes were 1,430 in modeling discovery, 351 in the internal holdout, 495 in the LC Sciences subset, and 492 in Guangzhou. PRS, polygenic risk score. MRS, metagenomic risk score. LOS, Louisiana Osteoporosis Study. The external PRS uses UK Biobank weights. The full base model includes age, sex, weight, height, recorded race (African American, White, and Asian indicators), sun exposure, regular exercise, smoking, alcohol consumption, cheese, milk, yogurt, and the first two genetic principal components. The Guangzhou base model was trained in modeling discovery with age, weight, height, sex, Asian race, alcohol consumption, and regular exercise.",
  "<sup>a</sup> &beta; is the adjusted difference in total hip bone mineral density (BMD, g/cm&sup2;) per SD of the score. Cohort-specific linear regressions omit constant or redundant predictors. P values are from two-sided t tests. A dash denotes the absence of a single score coefficient and P value for the combined model.",
  "<sup>b</sup> &Delta;R&sup2; is the R&sup2; with the indicated score(s) minus the base-model R&sup2; in the same participants. Raw R&sup2; is 1 minus the ratio of squared prediction error to total sum of squares. Both &Delta;R&sup2; intervals use 2,000 paired bootstrap resamples. An asterisk indicates an interval that excludes zero.",
  "<sup>c</sup> Recalibrated R&sup2; is the R&sup2; from a least-squares regression of observed on predicted BMD with an intercept in the evaluation cohort. It equals the squared Pearson correlation. n/a denotes discovery estimates from out-of-fold predictions, where recalibration was not evaluated.",
  "Evaluation with this covariate specification was post hoc. The secondary MRS used balances with nominal P < 0.05.")

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
