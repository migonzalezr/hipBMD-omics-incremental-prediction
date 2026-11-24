# ======================================================================
# 06_prs_permutation.R — PRS permutation test
#
# Permutes each selected score 150 times within the discovery set and recomputes the cross-validated increment.
#
# Inputs : selected PRS; discovery data
# Outputs: results/prs_permutations.csv; prs_permutations_summary.csv
# Paper  : Figure 2c
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)

source("config/paths.R")
P <- list(
  disc = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  prs  = file.path(BASE,"derived/prs_los/los_prs.PT_0.01.sscore"),
  prsq = file.path(BASE,"derived/prs_external/external_prs.PT_0.1.sscore"),
  pca  = file.path(BASE,"data/genotypes/los_qc_pca.eigenvec"),
  clin = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  bmd  = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  out  = file.path(BASE,"results")
)
falta <- unlist(P)[!file.exists(unlist(P))]
if (length(falta)) stop("no encuentro:\n  ", paste(falta, collapse="\n  "))

NPERM <- 150; SEED <- 2026; K <- 5; R <- 5

leer_s <- function(f,nm){ s <- fread(file=f); setnames(s, sub("^#","",names(s)))
  o <- data.table(IID=trimws(as.character(s$IID)), v=as.numeric(s$SCORE1_AVG)); setnames(o,"v",nm); o }

disc <- trimws(readLines(P$disc))
pca <- fread(file=P$pca, header=FALSE)
setnames(pca, c("FID","IID", paste0("PC", seq_len(ncol(pca)-2))))
pca <- pca[, .(IID=trimws(as.character(IID)), PC1, PC2)]
cl <- fread(file=P$clin); cl[, IID := trimws(as.character(Subject_ID))]
CLIN <- c("Age","gender_dummy","Weight_kg","Height_cm","Race_African","Race_Caucasian","Race_Asian",
          "Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
          "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")
cl <- unique(cl[, c("IID",CLIN), with=FALSE], by="IID")
b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)),
                        y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")

d <- Reduce(function(x,z) merge(x,z,by="IID"),
            list(data.table(IID=disc), ph, cl, pca,
                 leer_s(P$prs,"PRS"), leer_s(P$prsq,"PRS_qian")))
d <- d[complete.cases(d)]
cat("discovery:", nrow(d), "\n")

X <- model.matrix(~ ., data=as.data.frame(d[, c(CLIN,"PC1","PC2"), with=FALSE]))[,-1,drop=FALSE]
colnames(X) <- make.names(colnames(X), unique=TRUE)
base <- data.table(IID=d$IID, y=d$y, as.data.table(X))
vars <- setdiff(names(base), c("IID","y"))
f_base <- as.formula(paste("y ~", paste(vars, collapse=" + ")))
f_prs  <- as.formula(paste("y ~", paste(c(vars,"PRS"), collapse=" + ")))

set.seed(SEED)
n <- nrow(base)
# Keep at least one member of the rare reference group in each training fold.
reference <- which(rowSums(X[, c("Race_AfricanYes", "Race_CaucasianYes",
                                 "Race_AsianYes"), drop = FALSE]) == 0L)
stopifnot(length(reference) >= 2L)
folds <- lapply(seq_len(R), function(r) {
  fold <- sample(rep(seq_len(K), length.out = n))
  if (length(unique(fold[reference])) == 1L) {
    swap <- c(reference[1L], which(fold != fold[reference[1L]])[1L])
    fold[swap] <- rev(fold[swap])
  }
  fold
})

r2_oof <- function(y,yhat) 1 - sum((y-yhat)^2)/sum((y-mean(y))^2)
cv_r2 <- function(dat, form, fold) {
  oof <- numeric(nrow(dat))
  for (k in sort(unique(fold))) {
    tr <- fold != k
    oof[!tr] <- predict(lm(form, data=dat[tr]), newdata=dat[!tr], rankdeficient = "NAwarn")
  }
  r2_oof(dat$y, oof)
}

r2b <- vapply(seq_len(R), function(r) cv_r2(base, f_base, folds[[r]]), numeric(1))
cat("R2 base por repeticion:", paste(round(r2b,4), collapse=" "), "\n\n")

dR2_de <- function(sc) {
  dp <- copy(base); dp[, PRS := as.numeric(scale(sc))]
  mean(vapply(seq_len(R), function(r) cv_r2(dp, f_prs, folds[[r]]) - r2b[r], numeric(1)))
}

SCORES <- list(LOS_PT0.01 = d$PRS, QIAN_PT0.1 = d$PRS_qian)
nulas <- list(); resumen <- data.table()
for (nm in names(SCORES)) {
  sc  <- SCORES[[nm]]
  obs <- dR2_de(sc)
  cat(sprintf("%s: dR2 observado = %.5f -- corriendo %d permutaciones...\n", nm, obs, NPERM))
  set.seed(SEED)
  t0 <- Sys.time()
  nul <- vapply(seq_len(NPERM), function(i) {
    if (i %% 25 == 0) cat("   ", i, "/", NPERM, "\n")
    dR2_de(sample(sc))
  }, numeric(1))
  cat(sprintf("   %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

  pemp <- (1 + sum(nul >= obs)) / (1 + NPERM)
  nulas[[nm]] <- data.table(score=nm, perm=seq_len(NPERM), dR2_nulo=nul)
  resumen <- rbind(resumen, data.table(
    score = nm, dR2_observado = round(obs,5),
    nula_media = round(mean(nul),5), nula_sd = round(sd(nul),5),
    nula_p95 = round(quantile(nul,0.95),5), nula_max = round(max(nul),5),
    n_nulas_mayores = sum(nul >= obs), n_perm = NPERM,
    p_empirico = round(pemp,4),
    z = round((obs-mean(nul))/sd(nul),2)))
}

cat("\n=== test de permutacion del dR2 del PRS (discovery, n =", nrow(base), ") ===\n")
print(resumen, row.names=FALSE)
cat("\nLa media nula es negativa por construccion: anadir un predictor sin senal\n",
    "a un R2 fuera de muestra cuesta varianza sin aportar ajuste.\n")

fwrite(rbindlist(nulas), file.path(P$out,"prs_permutations.csv"))
fwrite(resumen,          file.path(P$out,"prs_permutations_summary.csv"))
cat("\nguardado en results/\n")
cat("===== FIN OK =====\n")
