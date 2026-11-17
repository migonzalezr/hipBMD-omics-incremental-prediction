# ======================================================================
# 15_balance_covariate_decomposition.R — Covariate decomposition of the balance signal
#
# Adjusts the balance-BMD associations for each covariate separately and in pairs.
#
# Inputs : discovery species and clinical data
# Outputs: results/balance_decomposition_covariates.csv; balance_decomposition_correlations.csv
# Paper  : Figure 3
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)
library(MASS); library(compositions); library(pls)

source("config/paths.R")
P <- list(
  mgs  = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  mast = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  bmd  = file.path(BASE,"data/phenotypes/los_hip_bmd.xlsx"),
  clin = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  disc = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  fun  = file.path(BASE,"external/PLS-PBs"),
  out  = file.path(BASE,"results")
)
for (f in c("fBalChip_PLS.r","fBPMaxOrthNewChip_PLS.R",
            "fBalChipman_PLS.R","fBPUpChi_PLS.r")) source(file.path(P$fun, f))
SEED <- 2026; PREV <- 0.01; FDRCUT <- 0.05
CLIN5 <- c("Age","gender_dummy","Weight_kg","Height_cm","BMI")

disc <- unique(trimws(readLines(P$disc)))
mgs <- fread(file = P$mgs); setnames(mgs, 1, "Metagenome_ID")
mgs <- fix_metagenome_ids(mgs)
pte <- as.data.table(read_excel(P$mast, sheet = "bgi"))
setnames(pte, "BGI_1794", "Subject_ID")
pte <- pte[, .(Metagenome_ID = trimws(as.character(Metagenome_ID)),
               IID = trimws(as.character(Subject_ID)))]
dt <- merge(pte[IID %in% disc], mgs, by = "Metagenome_ID")
b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID = trimws(as.character(b$Sample_ID)),
                        y = suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
cl <- fread(file = P$clin); cl[, IID := trimws(as.character(Subject_ID))]
cl <- unique(cl[, c("IID", CLIN5), with = FALSE], by = "IID")
meta <- Reduce(function(x,z) merge(x,z,by="IID"), list(dt[,.(IID)], ph, cl))
meta <- meta[complete.cases(meta)]
dt <- dt[IID %in% meta$IID]
setkey(dt,"IID"); setkey(meta,"IID"); dt <- dt[meta$IID]
stopifnot(identical(dt$IID, meta$IID))
RA <- as.matrix(dt[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RA) <- dt$IID
y <- meta$y
cat(sprintf("discovery: %d individuos, %d especies crudas\n", nrow(RA), ncol(RA)))

sp <- colnames(RA)[apply(RA,2,var,na.rm=TRUE) > 0]
sp <- sp[colMeans(RA[,sp,drop=FALSE] != 0) >= PREV]
X  <- RA[, sp, drop=FALSE]; pseudo <- 0.5*min(X[X!=0])
X[X==0] <- pseudo; X <- (X/rowSums(X))*100
cat(sprintf("especies tras filtros: %d\n", ncol(X)))
set.seed(SEED)
bal <- fBalChip_PLS(X, y, version = "cov")
Bt  <- t(bal$bal)[sp, , drop=FALSE]
PB  <- as.matrix(compositions::clr(X)) %*% Bt
cat(sprintf("balances: %d\n", ncol(PB)))

cov0 <- data.table(y = y, meta[, CLIN5, with=FALSE])

correr <- function(rhs) {
  res <- rbindlist(lapply(colnames(PB), function(nm) {
    d <- cbind(cov0, pb = PB[, nm])
    co <- summary(lm(as.formula(paste("y ~ pb", rhs)), data = d))$coefficients
    if (!"pb" %in% rownames(co)) return(data.table(PB_ID=nm, p=NA_real_, beta=NA_real_))
    data.table(PB_ID = nm, p = co["pb","Pr(>|t|)"], beta = co["pb","Estimate"])
  }))
  res[, FDR := p.adjust(p, method="BH")][]
}

esq <- list(
  "sin ajuste"                = "",
  "solo edad"                 = "+ Age",
  "solo sexo"                 = "+ gender_dummy",
  "solo peso"                 = "+ Weight_kg",
  "solo talla"                = "+ Height_cm",
  "solo BMI"                  = "+ BMI",
  "edad + sexo"               = "+ Age + gender_dummy",
  "edad + peso"               = "+ Age + Weight_kg",
  "sexo + peso"               = "+ gender_dummy + Weight_kg",
  "peso + talla"              = "+ Weight_kg + Height_cm",
  "edad + sexo + peso (esq.B)"= "+ Age + gender_dummy + Weight_kg"
)

ref <- correr("")
top0 <- ref[which.min(p), PB_ID]
cat(sprintf("\nbalance top sin ajustar: %s (p = %.3g)\n", top0, min(ref$p)))

tab <- rbindlist(lapply(names(esq), function(k) {
  r <- correr(esq[[k]])
  data.table(ajuste = k, formula = ifelse(esq[[k]]=="", "y ~ pb", paste("y ~ pb", esq[[k]])),
             n_FDR05 = sum(r$FDR < FDRCUT, na.rm=TRUE),
             n_p05   = sum(r$p   < 0.05,  na.rm=TRUE),
             FDR_min = signif(min(r$FDR, na.rm=TRUE), 3),
             top     = r[which.min(p), PB_ID],
             p_del_top0  = signif(r[PB_ID==top0, p], 3),
             beta_del_top0 = signif(r[PB_ID==top0, beta], 3))
}))
cat("\n=== CUANTOS BALANCES PASAN, QUITANDO UNA COVARIABLE A LA VEZ ===\n")
print(tab, row.names = FALSE)
fwrite(tab, file.path(P$out, "balance_decomposition_covariates.csv"))

pb_top <- PB[, top0]
cors <- rbindlist(lapply(c("y", CLIN5), function(v) {
  x <- if (v=="y") y else meta[[v]]
  if (is.character(x) || is.factor(x)) x <- as.integer(factor(x))
  data.table(variable = v,
             r_con_balance_top = round(cor(pb_top, x, use="complete.obs"), 3),
             r_con_BMD         = round(cor(y, x, use="complete.obs"), 3))
}))
cat(sprintf("\n=== correlaciones del balance %s y del BMD con cada covariable ===\n", top0))
print(cors, row.names = FALSE)
fwrite(cors, file.path(P$out, "balance_decomposition_correlations.csv"))

cat("\nLectura:\n",
    " si 'solo peso' ya tumba el conteo a ~0 -> la senal ES masa corporal\n",
    " si lo tumba 'solo edad' -> es confusion por envejecimiento, otra historia\n",
    " si hace falta la combinacion -> no hay un unico responsable y hay que decirlo\n")
cat("\nguardado en results/\n")
cat("===== FIN OK =====\n")
