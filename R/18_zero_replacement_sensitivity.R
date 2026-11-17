# ======================================================================
# 18_zero_replacement_sensitivity.R — Zero-replacement sensitivity
#
# Repeats the species-level associations under four zero-replacement rules.
#
# Inputs : discovery species and clinical data
# Outputs: results/zero_replacement_summary.csv; zero_replacement_species.csv
# Paper  : Table S9
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table); library(readxl); library(compositions); library(zCompositions)
})

source("config/paths.R")
P <- list(
  mgs   = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  mast  = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  bmd   = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  clinB = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  disc  = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  ref   = file.path(BASE,"results/species_association.csv"),
  out   = file.path(BASE,"results")
)
falta <- unlist(P)[!file.exists(unlist(P))]
if (length(falta)) stop("no encuentro:\n  ", paste(falta, collapse="\n  "))

PREV <- 0.01
CLIN <- c("Age","gender_dummy","Weight_kg")
sp_std <- function(v){ s <- sub(".*(s__[^\\.]+).*","\\1",v); s[!grepl("s__",v)] <- NA; s }

b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)),
                        y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
disc <- unique(trimws(readLines(P$disc)))
mgs <- fread(file=P$mgs); setnames(mgs,1,"Metagenome_ID")
mgs <- fix_metagenome_ids(mgs)
pteB <- as.data.table(read_excel(P$mast, sheet="bgi"))
setnames(pteB,"BGI_1794","Subject_ID")
pteB <- pteB[, .(Metagenome_ID=trimws(as.character(Metagenome_ID)), IID=trimws(as.character(Subject_ID)))]
dtB <- merge(pteB[IID %in% disc], mgs, by="Metagenome_ID")
clB <- fread(file=P$clinB); clB[, IID := trimws(as.character(Subject_ID))]
clB <- unique(clB[, c("IID",CLIN), with=FALSE], by="IID")
mD <- Reduce(function(x,z) merge(x,z,by="IID"), list(dtB[,.(IID)], ph, clB))
mD <- mD[complete.cases(mD)]
dtB <- dtB[IID %in% mD$IID]; setkey(dtB,"IID"); setkey(mD,"IID"); dtB <- dtB[mD$IID]
stopifnot(identical(dtB$IID, mD$IID))
RAD <- as.matrix(dtB[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RAD) <- dtB$IID
colnames(RAD) <- sp_std(colnames(RAD)); RAD <- RAD[, !is.na(colnames(RAD)), drop=FALSE]

n_total <- ncol(RAD)
sp_var  <- colnames(RAD)[apply(RAD,2,var,na.rm=TRUE) > 0]
sp428   <- sp_var[colMeans(RAD[,sp_var,drop=FALSE] != 0) >= PREV]
stopifnot(length(sp428) == 428)
X0 <- RAD[, sp428, drop=FALSE]
cat(sprintf("discovery: %d individuos\n", nrow(X0)))
cat(sprintf("especies en el archivo: %d | con varianza > 0: %d | prevalencia >= 1%%: %d\n",
            n_total, length(sp_var), length(sp428)))
cat(sprintf("fraccion de ceros en la matriz 1,430 x 428: %.1f%%\n", 100*mean(X0 == 0)))
cat(sprintf("especies con > 80%% de ceros (las que multRepl borraria por defecto): %d\n",
            sum(colMeans(X0 == 0) > 0.8)))
mn_glob <- min(X0[X0 != 0])
cat(sprintf("minimo no nulo de toda la matriz: %.3g\n", mn_glob))

cierre <- function(M) (M / rowSums(M)) * 100
regla_A <- function(X) { X[X == 0] <- 0.5 * min(X[X != 0]); cierre(X) }
regla_B <- function(X) {
  dl <- apply(X, 2, function(v) min(v[v != 0]))
  R  <- multRepl(X, label = 0, dl = dl, frac = 0.65, z.warning = 1, z.delete = FALSE)
  R  <- as.matrix(R)
  if (ncol(R) != ncol(X)) stop("multRepl elimino especies: el universo de 428 no se conserva")
  cierre(R)
}
regla_C <- function(X) {
  for (j in seq_len(ncol(X))) { v <- X[, j]; v[v == 0] <- 0.5 * min(v[v != 0]); X[, j] <- v }
  cierre(X)
}
regla_D <- function(X) { Pp <- X / rowSums(X); cierre(Pp + min(Pp[Pp != 0])) }
REGLAS <- list(
  A_global_half_min        = regla_A,
  B_multiplicative         = regla_B,
  C_per_species_half_min   = regla_C,
  D_constant_added_to_all  = regla_D)

cov_sp <- data.table(y=mD$y, Age=mD$Age, sexo=mD$gender_dummy, Weight_kg=mD$Weight_kg)
mwas <- function(Xr) {
  Z <- as.matrix(compositions::clr(Xr))
  r <- rbindlist(lapply(colnames(Z), function(nm){
    d <- cbind(cov_sp, sp=Z[,nm])
    m0 <- summary(lm(y ~ sp, data=d))$coefficients
    m1 <- summary(lm(y ~ sp + Age + sexo + Weight_kg, data=d))$coefficients
    data.table(especie=nm, beta_crudo=m0["sp","Estimate"], p_crudo=m0["sp","Pr(>|t|)"],
               beta_ajust=m1["sp","Estimate"], p_ajust=m1["sp","Pr(>|t|)"])
  }))
  r[, `:=`(FDR_crudo=p.adjust(p_crudo,"BH"), FDR_ajust=p.adjust(p_ajust,"BH"))][]
}
res <- rbindlist(lapply(names(REGLAS), function(nm) {
  cat(sprintf("  regla %s ...\n", nm))
  mwas(REGLAS[[nm]](X0))[, regla := nm][]
}))

ref <- fread(file = P$ref)
chk <- merge(res[regla == "A_global_half_min"], ref, by = "especie", suffixes = c("", "_ref"))
stopifnot(nrow(chk) == 428)
dif <- chk[, max(abs(c(beta_crudo - beta_crudo_ref, beta_ajust - beta_ajust_ref,
                       p_crudo - p_crudo_ref, p_ajust - p_ajust_ref)))]
if (dif > 1e-10) stop(sprintf("la regla A no reproduce species_association.csv (dif %.2e)", dif))
cat(sprintf("\ncontrol OK: la regla A reproduce species_association.csv (dif max %.1e)\n", dif))

A <- res[regla == "A_global_half_min", .(especie, la = -log10(p_ajust), lc = -log10(p_crudo))]
resumen <- res[, {
  m <- merge(.SD[, .(especie, a2 = -log10(p_ajust), c2 = -log10(p_crudo))], A, by = "especie")
  .(n_FDR05_crudo = sum(FDR_crudo < 0.05), FDR_min_crudo = signif(min(FDR_crudo), 3),
    n_FDR05_ajust = sum(FDR_ajust < 0.05), FDR_min_ajust = signif(min(FDR_ajust), 3),
    top_ajust = especie[which.min(p_ajust)],
    rho_vs_A_crudo = round(cor(m$c2, m$lc, method = "spearman"), 3),
    rho_vs_A_ajust = round(cor(m$a2, m$la, method = "spearman"), 3),
    FDR_ajust_S_australis = signif(FDR_ajust[grepl("australis", especie)], 3))
}, by = regla]
cat("\n=== sensibilidad al reemplazo de ceros (discovery, n = 1,430, 428 especies) ===\n")
print(resumen, row.names = FALSE)
cat("\nrho_vs_A: correlacion de Spearman de -log10(P) entre cada regla y la del manuscrito.\n")

fwrite(resumen, file.path(P$out, "zero_replacement_summary.csv"))
fwrite(res, file.path(P$out, "zero_replacement_species.csv"))
cat("\nguardado en results/zero_replacement_summary.csv y zero_replacement_species.csv\n")
cat("===== FIN OK =====\n")
