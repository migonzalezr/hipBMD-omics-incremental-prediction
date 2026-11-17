# ======================================================================
# 05a_prs_threshold_los.R — PRS threshold selection, LOS weights
#
# Selects the P-value threshold of the LOS-derived PRS by the cross-validated increment in explained variance over the base model, in the modeling discovery set only (5 x 5-fold CV, seed 2026).
#
# Inputs : derived/prs_los/los_prs.PT_*.sscore (plink/03_prs_los.sh)
# Outputs: results/prs_los_threshold_selection.csv; prs_los_threshold_frozen.txt
# Paper  : Figure 2a; Table S3
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)

source("config/paths.R")
P <- list(
  disc = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  prs  = file.path(BASE,"derived/prs_los"),
  pca  = file.path(BASE,"data/genotypes/los_qc_pca.eigenvec"),
  clin = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  bmd  = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  out  = file.path(BASE,"results")
)
dir.create(P$out, showWarnings = FALSE, recursive = TRUE)

disc <- trimws(readLines(P$disc))
cat("discovery BGI:", length(disc), "\n")

arch <- list.files(P$prs, pattern = "^los_prs\\..*\\.sscore$", full.names = TRUE)
stopifnot(length(arch) > 0)

leer_score <- function(f) {
  s <- fread(f); setnames(s, sub("^#", "", names(s)))
  iid <- grep("^IID$", names(s), value = TRUE)[1]
  sc  <- grep("^SCORE1_AVG$", names(s), value = TRUE)[1]
  ac  <- grep("^ALLELE_CT$", names(s), value = TRUE)[1]
  data.table(IID = as.character(s[[iid]]), PRS = as.numeric(s[[sc]]),
             ALLELE_CT = as.numeric(s[[ac]]),
             PT = sub("\\.sscore$", "", sub("^los_prs\\.", "", basename(f))))
}
scores <- rbindlist(lapply(arch, leer_score))

chk <- scores[, .(n = .N, variantes = round(mean(ALLELE_CT)/2),
                  sd = sd(PRS, na.rm = TRUE), na = sum(is.na(PRS))), by = PT][order(variantes)]
cat("\n=== PRS por umbral ===\n"); print(chk)

malos <- chk[is.na(sd) | sd < 1e-12 | na > 0, PT]
if (length(malos)) {
  cat("\ndescartados (sin varianza o con NA):", paste(malos, collapse = " "), "\n")
  scores <- scores[!PT %in% malos]
}
umbrales <- sort(unique(scores$PT))
cat("entran a la comparacion:", paste(umbrales, collapse = " "), "\n")

pca <- fread(P$pca, header = FALSE)
setnames(pca, c("FID","IID", paste0("PC", seq_len(ncol(pca)-2))))
pca <- pca[, .(IID = trimws(as.character(IID)), PC1, PC2)]

cl <- fread(P$clin); cl[, IID := trimws(as.character(Subject_ID))]
CLIN <- c("Age","gender_dummy","Weight_kg","Height_cm","BMI",
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
cat("clase de cada covariable clinica:\n")
print(sapply(d[, CLIN, with = FALSE], function(x) paste(class(x)[1], "niveles:", uniqueN(x))))

X <- model.matrix(~ ., data = as.data.frame(d[, c(CLIN,"PC1","PC2"), with = FALSE]))[, -1, drop = FALSE]
colnames(X) <- make.names(colnames(X), unique = TRUE)
base <- data.table(IID = d$IID, y = d$y, as.data.table(X))
vars <- setdiff(names(base), c("IID","y"))
cat("predictores en el bloque base:", length(vars), "\n")

f_base <- as.formula(paste("y ~", paste(vars, collapse = " + ")))
f_prs  <- as.formula(paste("y ~", paste(c(vars,"PRS"), collapse = " + ")))

set.seed(2026)
K <- 5; R <- 5; n <- nrow(base)
folds <- lapply(seq_len(R), function(r) sample(rep(seq_len(K), length.out = n)))

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

res <- rbindlist(lapply(umbrales, function(pt) {
  dp <- merge(base, scores[PT == pt, .(IID, PRS)], by = "IID", sort = FALSE)
  stopifnot(nrow(dp) == nrow(base))
  dp[, PRS := as.numeric(scale(PRS))]
  idx <- match(dp$IID, base$IID)
  rbindlist(lapply(seq_len(R), function(r) {
    fo <- folds[[r]][idx]
    data.table(PT = pt, rep = r,
               R2_base = cv_r2(dp, f_base, fo),
               R2_prs  = cv_r2(dp, f_prs,  fo))
  }))
}))
res[, dR2 := R2_prs - R2_base]

tab <- merge(res[, .(R2_base = mean(R2_base), R2_prs = mean(R2_prs),
                     dR2 = mean(dR2), sd = sd(dR2)), by = PT],
             chk[, .(PT, variantes)], by = "PT")[order(-dR2)]
cat("\n=== DeltaR2 del PRS por umbral (CV 5x5, discovery) ===\n")
print(tab, digits = 4)

ganador <- tab$PT[1]
cat("\n>>> UMBRAL ELEGIDO:", ganador,
    "| dR2 =", round(tab$dR2[1], 4), "+/-", round(tab$sd[1], 4),
    "| variantes:", tab$variantes[1], "<<<\n")
cat("    elegido SOLO con discovery; holdout, LC y China intactos\n")
if (nrow(tab) > 1)
  cat("    ventaja sobre el segundo:", round(tab$dR2[1] - tab$dR2[2], 5),
      "| sd entre repeticiones:", round(tab$sd[1], 5), "\n")
if (tab$dR2[1] > 0.10)
  cat("\n*** ALTO: dR2 > 0.10 supera la regla de paro. Auditar por fuga. ***\n")

fwrite(tab, file.path(P$out, "prs_los_threshold_selection.csv"))
writeLines(c(sprintf("selected_threshold=%s", ganador),
             sprintf("dR2_cv=%.5f", tab$dR2[1]),
             sprintf("sd=%.5f", tab$sd[1]),
             sprintf("variantes=%d", tab$variantes[1]),
             sprintf("n_discovery=%d", nrow(base)),
             sprintf("semilla=2026"), sprintf("fecha=%s", Sys.Date())),
           file.path(P$out, "prs_los_threshold_frozen.txt"))
cat("\nguardado en results/\n")
