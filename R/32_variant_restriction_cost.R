# ======================================================================
# 32_variant_restriction_cost.R — Cost of restricting each PRS to the Guangzhou variants
# #
# Delta R2 of the frozen, transferable and Guangzhou-restricted scores in modeling discovery (5 x 5-fold CV, identical folds).
# #
# Inputs : derived/prs_los, derived/prs_external, derived/prs_reduced
# Outputs: results/coverage_control_levels.csv, coverage_control_paired.csv
# Paper  : Table S12
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
rm(list = ls()); gc()
library(data.table); library(readxl)

source("config/paths.R")
P <- list(
  disc = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  prs  = file.path(BASE,"derived/prs_los"),
  qian = file.path(BASE,"derived/prs_external"),
  red  = file.path(BASE,"derived/prs_restricted"),
  pca  = file.path(BASE,"data/genotypes/los_qc_pca.eigenvec"),
  clin = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  bmd  = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  out  = file.path(BASE,"results")
)
dir.create(P$out, showWarnings = FALSE, recursive = TRUE)

disc <- trimws(readLines(P$disc))
cat("discovery BGI:", length(disc), "\n")

SC <- data.table(
  fuente = c("LOS","LOS","LOS","QIAN","QIAN","QIAN"),
  nivel  = c("1_congelado","2_transferible","3_reducido",
             "1_congelado","2_transferible","3_reducido"),
  arch   = c(file.path(P$prs,  "los_prs.PT_0.01.sscore"),
             file.path(P$red,  "los_weights_transferable.sscore"),
             file.path(P$red,  "los_weights_guangzhou_variants.sscore"),
             file.path(P$qian, "external_prs.PT_0.1.sscore"),
             file.path(P$red,  "external_weights_transferable.sscore"),
             file.path(P$red,  "external_weights_guangzhou_variants.sscore"))
)
falta <- SC$arch[!file.exists(SC$arch)]
if (length(falta)) stop("faltan archivos de score:\n  ", paste(falta, collapse = "\n  "))

leer_score <- function(f) {
  s <- fread(f); setnames(s, sub("^#", "", names(s)))
  data.table(IID = trimws(as.character(s[["IID"]])),
             PRS = as.numeric(s[["SCORE1_AVG"]]),
             ALLELE_CT = as.numeric(s[["ALLELE_CT"]]))
}
scores <- rbindlist(lapply(seq_len(nrow(SC)), function(i) {
  x <- leer_score(SC$arch[i]); x[, `:=`(fuente = SC$fuente[i], nivel = SC$nivel[i])]; x
}))

chk <- scores[, .(n = .N, variantes = round(mean(ALLELE_CT)/2),
                  sd_score = sd(PRS, na.rm = TRUE), na = sum(is.na(PRS))),
              by = .(fuente, nivel)][order(fuente, nivel)]
cat("\n=== los seis scores, en los 6,732 de LOS.final ===\n"); print(chk)
if (any(is.na(chk$sd_score) | chk$sd_score < 1e-12)) stop("algun score es constante")

pca <- fread(P$pca, header = FALSE)
setnames(pca, c("FID","IID", paste0("PC", seq_len(ncol(pca)-2))))
pca <- pca[, .(IID = trimws(as.character(IID)), PC1, PC2)]

cl <- fread(P$clin); cl[, IID := trimws(as.character(Subject_ID))]
CLIN <- c("Age","gender_dummy","Weight_kg","Height_cm","Race_African","Race_Caucasian","Race_Asian",
          "Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
          "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")
cl <- unique(cl[, c("IID", CLIN), with = FALSE], by = "IID")

b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID = trimws(as.character(b$Sample_ID)),
                        y = suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by = "IID")

d <- Reduce(function(x,y) merge(x, y, by = "IID"),
            list(data.table(IID = disc), ph, cl, pca))
d <- d[complete.cases(d)]
cat("\ndiscovery con covariables completas:", nrow(d), "de", length(disc), "\n")

X <- model.matrix(~ ., data = as.data.frame(d[, c(CLIN,"PC1","PC2"), with = FALSE]))[, -1, drop = FALSE]
colnames(X) <- make.names(colnames(X), unique = TRUE)
base <- data.table(IID = d$IID, y = d$y, as.data.table(X))
vars <- setdiff(names(base), c("IID","y"))
cat("predictores en el bloque base:", length(vars), "\n")

f_base <- as.formula(paste("y ~", paste(vars, collapse = " + ")))
f_prs  <- as.formula(paste("y ~", paste(c(vars,"PRS"), collapse = " + ")))

set.seed(2026)
K <- 5; R <- 5; n <- nrow(base)
reference <- which(rowSums(X[, c("Race_AfricanYes","Race_CaucasianYes","Race_AsianYes"), drop = FALSE]) == 0L)
stopifnot(length(reference) >= 2L)
folds <- lapply(seq_len(R), function(r) {
  fold <- sample(rep(seq_len(K), length.out = n))
  if (length(unique(fold[reference])) == 1L) {
    swap <- c(reference[1L], which(fold != fold[reference[1L]])[1L]); fold[swap] <- rev(fold[swap]) }
  fold
})

r2_oof <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
cv_r2 <- function(dat, form, fold) {
  oof <- numeric(nrow(dat))
  for (k in sort(unique(fold))) {
    tr <- fold != k
    m  <- lm(form, data = dat[tr])
    oof[!tr] <- predict(m, newdata = dat[!tr])
  }
  r2_oof(dat$y, oof)
}

res <- rbindlist(lapply(seq_len(nrow(SC)), function(i) {
  s  <- scores[fuente == SC$fuente[i] & nivel == SC$nivel[i], .(IID, PRS)]
  dp <- merge(base, s, by = "IID", sort = FALSE)
  stopifnot(nrow(dp) == nrow(base))
  dp[, PRS := as.numeric(scale(PRS))]
  idx <- match(dp$IID, base$IID)
  rbindlist(lapply(seq_len(R), function(r) {
    fo <- folds[[r]][idx]
    data.table(fuente = SC$fuente[i], nivel = SC$nivel[i], rep = r,
               R2_base = cv_r2(dp, f_base, fo),
               R2_prs  = cv_r2(dp, f_prs,  fo))
  }))
}))
res[, dR2 := R2_prs - R2_base]

tab <- merge(res[, .(R2_base = mean(R2_base), R2_prs = mean(R2_prs),
                     dR2 = mean(dR2), sd = sd(dR2)), by = .(fuente, nivel)],
             chk[, .(fuente, nivel, variantes)], by = c("fuente","nivel"))[order(fuente, nivel)]
cat("\n=== DeltaR2 por nivel de recorte (CV 5x5, discovery) ===\n")
print(tab, digits = 4)

w <- dcast(res, fuente + rep ~ nivel, value.var = "dR2")
pares <- rbindlist(lapply(c("LOS","QIAN"), function(fu) {
  x <- w[fuente == fu]
  rbind(
    data.table(fuente = fu, comparacion = "congelado -> transferible",
               dif = mean(x$`2_transferible` - x$`1_congelado`),
               sd  = sd(x$`2_transferible` - x$`1_congelado`)),
    data.table(fuente = fu, comparacion = "transferible -> reducido",
               dif = mean(x$`3_reducido` - x$`2_transferible`),
               sd  = sd(x$`3_reducido` - x$`2_transferible`)),
    data.table(fuente = fu, comparacion = "congelado -> reducido",
               dif = mean(x$`3_reducido` - x$`1_congelado`),
               sd  = sd(x$`3_reducido` - x$`1_congelado`))
  )
}))
pares[, en_sd := dif / sd]
pares[, veredicto := fifelse(abs(en_sd) < 1, "equivalentes",
                      fifelse(abs(en_sd) > 2 & dif < 0, "el recorte SI cuesta",
                      fifelse(abs(en_sd) > 2 & dif > 0, "el recorte MEJORA (revisar)",
                              "no concluyente")))]
cat("\n=== diferencias pareadas entre niveles (misma fold, misma repeticion) ===\n")
print(pares, digits = 4)

cors <- rbindlist(lapply(c("LOS","QIAN"), function(fu) {
  m <- dcast(scores[fuente == fu & IID %in% base$IID], IID ~ nivel, value.var = "PRS")
  data.table(fuente = fu,
             r_cong_trans = cor(m$`1_congelado`,   m$`2_transferible`),
             r_trans_red  = cor(m$`2_transferible`, m$`3_reducido`),
             r_cong_red   = cor(m$`1_congelado`,   m$`3_reducido`))
}))
cat("\n=== correlacion entre niveles del score (discovery) ===\n"); print(cors, digits = 4)
cat("una r de exactamente 1 significaria que el recorte no hizo nada\n")

cat("\n=== LECTURA ===\n")
for (fu in c("LOS","QIAN")) {
  v <- pares[fuente == fu & comparacion == "transferible -> reducido"]
  cat(sprintf("%s: transferible -> reducido  dif = %+.5f (%.1f sd)  -> %s\n",
              fu, v$dif, v$en_sd, v$veredicto))
  if (v$veredicto == "equivalentes")
    cat("    la cobertura queda DESCARTADA como explicacion de lo que salga en Guangzhou\n")
  if (v$veredicto == "el recorte SI cuesta")
    cat("    un dR2 nulo en Guangzhou sera atribuible en parte a la cobertura; reportar ambas causas\n")
}

fwrite(tab,   file.path(P$out, "coverage_control_levels.csv"))
fwrite(pares, file.path(P$out, "coverage_control_paired.csv"))
fwrite(cors,  file.path(P$out, "coverage_control_correlations.csv"))
fwrite(res,   file.path(P$out, "coverage_control_by_repeat.csv"))
cat("\nguardado en results/\n")
cat("\n===== FIN OK =====\n")
