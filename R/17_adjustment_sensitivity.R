# ======================================================================
# 17_adjustment_sensitivity.R — Adjustment sensitivity
#
# Repeats the community-level analyses with alternative representations of body size (weight and height, BMI, or all three).
#
# Inputs : as 16
# Outputs: results/adjustment_sensitivity.csv
# Paper  : Table S10
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl); library(vegan)

source("config/paths.R")
P <- list(
  mgs   = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  lcsp  = file.path(BASE,"data/microbiome/los_lc_species_abundance.csv"),
  chsp  = file.path(BASE,"data/microbiome/guangzhou_species_abundance.csv"),
  chcl  = file.path(BASE,"data/phenotypes/guangzhou_clinical.csv"),
  mast  = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  bmd   = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  clinB = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  clinL = file.path(BASE,"data/phenotypes/los_lc_clinical_imputed.csv"),
  disc  = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  lc    = file.path(BASE,"derived/analysis_sets/ids_lc_subset.txt"),
  out   = file.path(BASE,"results")
)
falta <- unlist(P)[!file.exists(unlist(P))]
if (length(falta)) stop("no encuentro:\n  ", paste(falta, collapse="\n  "))

SEED <- 2026; PREV <- 0.01; NPERM <- 999
CLIN <- c("Age","gender_dummy","Weight_kg","Height_cm","BMI")
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
sp428 <- colnames(RAD)[apply(RAD,2,var,na.rm=TRUE) > 0]
sp428 <- sp428[colMeans(RAD[,sp428,drop=FALSE] != 0) >= PREV]

lc_id <- unique(trimws(readLines(P$lc)))
pteL <- as.data.table(read_excel(P$mast, sheet="lc")); setnames(pteL,"LC_495","IID")
pteL <- pteL[, .(IID=trimws(as.character(IID)), Metagenome_ID=trimws(as.character(Metagenome_ID)))]
tabL <- fread(file=P$lcsp); ML <- as.matrix(tabL[, -1, with=FALSE])
rownames(ML) <- trimws(as.character(tabL[[1]]))
RAL <- t(ML); colnames(RAL) <- sp_std(colnames(RAL))
RAL <- RAL[, !is.na(colnames(RAL)), drop=FALSE]; RAL <- RAL[, !duplicated(colnames(RAL)), drop=FALSE]
RAL <- RAL[rownames(RAL) %in% pteL$Metagenome_ID, , drop=FALSE]
rownames(RAL) <- pteL$IID[match(rownames(RAL), pteL$Metagenome_ID)]
clL <- fread(file=P$clinL); clL[, IID := trimws(as.character(Subject_ID))]
clL <- unique(clL[, c("IID",CLIN), with=FALSE], by="IID")
mL <- Reduce(function(x,z) merge(x,z,by="IID"),
             list(data.table(IID=intersect(rownames(RAL), lc_id)), ph, clL))
mL <- mL[complete.cases(mL)]; RAL <- RAL[mL$IID, , drop=FALSE]

tabC <- fread(file=P$chsp); MC <- as.matrix(tabC[, -1, with=FALSE])
rownames(MC) <- trimws(as.character(tabC[[1]]))
RAC <- t(MC); colnames(RAC) <- sp_std(colnames(RAC))
RAC <- RAC[, !is.na(colnames(RAC)), drop=FALSE]; RAC <- RAC[, !duplicated(colnames(RAC)), drop=FALSE]
clC <- fread(file=P$chcl); clC[, IID := trimws(as.character(ID))]
clC[, `:=`(Age=as.numeric(Age), gender_dummy="Female", Weight_kg=as.numeric(Weight_kg),
           Height_cm=as.numeric(Height_cm), BMI=as.numeric(BMI), y=as.numeric(HTOT_BMD))]
mC <- clC[IID %in% rownames(RAC), c("IID","y",CLIN), with=FALSE]
mC <- mC[complete.cases(mC)]; RAC <- RAC[mC$IID, , drop=FALSE]

comunes <- Reduce(intersect, list(sp428, colnames(RAL), colnames(RAC)))
renorm <- function(M, sp){ X <- M[, sp, drop=FALSE]; X[X<0] <- 0; X/pmax(rowSums(X), .Machine$double.eps) }
COH <- list(
  `Discovery (BGI)`          = list(X=renorm(RAD,comunes), m=as.data.table(mD)),
  `Validation 1 (LC)`        = list(X=renorm(RAL,comunes), m=as.data.table(mL)),
  `Validation 2 (Guangzhou)` = list(X=renorm(RAC,comunes), m=as.data.table(mC)))
cat(sprintf("discovery %d | LC %d | Guangzhou %d | %d especies comunes\n",
            nrow(RAD), nrow(RAL), nrow(RAC), length(comunes)))

vif_bloque <- function(d){
  V <- as.data.frame(d[, .(Weight_kg, Height_cm, BMI)])
  rbindlist(lapply(names(V), function(cc){
    r2 <- summary(lm(as.formula(paste(cc,"~ .")), data=V))$r.squared
    data.table(variable=cc, R2_por_las_otras_dos=r2, VIF=1/(1-r2))
  }))
}
cat("\n=== colinealidad peso / talla / BMI, por cohorte ===\n")
vifs <- rbindlist(lapply(names(COH), function(co) vif_bloque(COH[[co]]$m)[, cohorte := co]))
print(vifs[, .(cohorte, variable, R2=round(R2_por_las_otras_dos,4), VIF=round(VIF,1))], row.names=FALSE)

SETS <- list(
  A = c("Age","sexo","Weight_kg","Height_cm","BMI"),
  B = c("Age","sexo","Weight_kg","Height_cm"),
  C = c("Age","sexo","BMI"))
ETQ <- c(A="age + sex + weight + height + BMI", B="age + sex + weight + height",
         C="age + sex + BMI")

out <- list()
for (co in names(COH)) {
  X <- COH[[co]]$X; m <- copy(COH[[co]]$m)[, sexo := gender_dummy]
  cat(sprintf("\n--- %s ---\n", co))
  D <- vegan::vegdist(X, method="bray")
  shv <- as.numeric(scale(vegan::diversity(X, index="shannon")))
  for (s in names(SETS)) {
    cv <- SETS[[s]]
    if (uniqueN(m$sexo) < 2) cv <- setdiff(cv, "sexo")
    f <- as.formula(paste("D ~", paste(c(cv,"y"), collapse=" + ")))
    set.seed(SEED)
    a <- vegan::adonis2(f, data=m, permutations=NPERM, by="margin")

    dd <- copy(m)[, shannon := shv]
    for (vv in setdiff(c(cv,"y"),"sexo")) set(dd, j=vv, value=as.numeric(scale(dd[[vv]])))
    fit <- lm(as.formula(paste("shannon ~", paste(c(cv,"y"), collapse=" + "))), data=dd)
    ci <- confint(fit)["y",]
    out[[length(out)+1]] <- data.table(
      cohorte=co, conjunto=s, covariables=ETQ[[s]], n=nrow(X),
      permanova_R2=a["y","R2"], permanova_p=a["y","Pr(>F)"],
      shannon_beta=unname(coef(fit)["y"]), shannon_lo=ci[1], shannon_hi=ci[2],
      shannon_p=summary(fit)$coefficients["y",4])
    cat(sprintf("  %s  %-30s  PERMANOVA R2 = %.5f (p = %.3f)   Shannon beta = %+.4f (p = %.3f)\n",
                s, ETQ[[s]], a["y","R2"], a["y","Pr(>F)"],
                coef(fit)["y"], summary(fit)$coefficients["y",4]))
  }
}
res <- rbindlist(out)
fwrite(res, file.path(P$out,"adjustment_sensitivity.csv"))

cat("\n=== resumen: R2 del BMD ajustado segun el conjunto de covariables ===\n")
print(dcast(res, cohorte ~ conjunto, value.var="permanova_R2"), row.names=FALSE)
cat("\n=== p correspondientes ===\n")
print(dcast(res, cohorte ~ conjunto, value.var="permanova_p"), row.names=FALSE)

rango <- res[, .(dif_max_R2 = max(permanova_R2) - min(permanova_R2),
                 mismos_signos = uniqueN(permanova_p < 0.05) == 1), by=cohorte]
cat("\n=== veredicto ===\n"); print(rango, row.names=FALSE)
if (all(rango$mismos_signos)) {
  cat("Los tres conjuntos llevan a la MISMA conclusion en las tres cohortes.\n",
      "La eleccion de covariables de ajuste no cambia el resultado y se puede\n",
      "declarar asi en Metodos, con esta tabla como respaldo.\n", sep="")
} else {
  cat("ATENCION: la conclusion DEPENDE del conjunto de covariables. Hay que\n",
      "reportar el conjunto B como primario y discutir la discrepancia.\n", sep="")
}

cat("\nguardado en results/adjustment_sensitivity.csv\n")
cat("===== FIN OK =====\n")
