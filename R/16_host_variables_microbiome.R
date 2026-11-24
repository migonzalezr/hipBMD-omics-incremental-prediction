# ======================================================================
# 16_host_variables_microbiome.R — Host variables and community composition
#
# Variance in community composition and Shannon diversity explained by each host variable, per cohort.
#
# Inputs : species abundances; clinical data
# Outputs: results/host_variables_permanova.csv; host_variables_shannon.csv
# Paper  : Figure S2
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
cat(sprintf("universo pre-declarado: %d especies\n", length(sp428)))

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

cat(sprintf("discovery %d | LC %d | Guangzhou %d individuos\n", nrow(RAD), nrow(RAL), nrow(RAC)))

comunes <- Reduce(intersect, list(sp428, colnames(RAL), colnames(RAC)))
cat(sprintf("especies comunes a las tres cohortes: %d\n", length(comunes)))
renorm <- function(M, sp){ X <- M[, sp, drop=FALSE]; X[X<0] <- 0; X/pmax(rowSums(X), .Machine$double.eps) }

COH <- list(
  `Discovery (BGI)`          = list(X=renorm(RAD,comunes), m=mD),
  `LC subset`        = list(X=renorm(RAL,comunes), m=mL),
  `Guangzhou cohort` = list(X=renorm(RAC,comunes), m=mC))

VARS <- c(Age="Age", gender_dummy="Sex", Weight_kg="Weight", Height_cm="Height",
          BMI="BMI", y="Hip BMD")

res_pm <- list()
for (co in names(COH)) {
  X <- COH[[co]]$X; m <- as.data.table(COH[[co]]$m)
  m[, sexo := gender_dummy]
  cat(sprintf("\n--- %s: Bray-Curtis sobre %d x %d ---\n", co, nrow(X), ncol(X)))
  D <- vegan::vegdist(X, method="bray")
  for (v in names(VARS)) {
    vv <- if (v == "gender_dummy") "sexo" else v
    if (uniqueN(m[[vv]]) < 2) {
      res_pm[[length(res_pm)+1]] <- data.table(
        cohorte=co, variable=VARS[[v]], modelo="marginal", n=nrow(X),
        R2=NA_real_, p=NA_real_, estimable=FALSE)
      next
    }
    set.seed(SEED)
    a <- vegan::adonis2(as.formula(paste("D ~", vv)), data=m, permutations=NPERM)
    res_pm[[length(res_pm)+1]] <- data.table(
      cohorte=co, variable=VARS[[v]], modelo="marginal", n=nrow(X),
      R2=a$R2[1], p=a$`Pr(>F)`[1], estimable=TRUE)
    cat(sprintf("  %-10s R2 = %.5f   p = %.3f\n", VARS[[v]], a$R2[1], a$`Pr(>F)`[1]))
  }

  rhs <- c("Age", if (uniqueN(m$sexo) > 1) "sexo", "Weight_kg","Height_cm","BMI","y")
  set.seed(SEED)
  aj <- vegan::adonis2(as.formula(paste("D ~", paste(rhs, collapse=" + "))),
                       data=m, permutations=NPERM, by="margin")
  res_pm[[length(res_pm)+1]] <- data.table(
    cohorte=co, variable="Hip BMD (adjusted)", modelo="ajustado", n=nrow(X),
    R2=aj["y","R2"], p=aj["y","Pr(>F)"], estimable=TRUE)
  cat(sprintf("  %-10s R2 = %.5f   p = %.3f\n", "BMD ajust", aj["y","R2"], aj["y","Pr(>F)"]))
}
pm <- rbindlist(res_pm)
fwrite(pm, file.path(P$out,"host_variables_permanova.csv"))

res_sh <- list()
for (co in names(COH)) {
  X <- COH[[co]]$X; m <- as.data.table(COH[[co]]$m)
  d <- data.table(shannon = as.numeric(scale(vegan::diversity(X, index="shannon"))),
                  Age=m$Age, sexo=m$gender_dummy, Weight_kg=m$Weight_kg,
                  Height_cm=m$Height_cm, BMI=m$BMI, y=m$y)
  for (v in names(VARS)) {
    vv <- if (v == "gender_dummy") "sexo" else v
    if (uniqueN(d[[vv]]) < 2) {
      res_sh[[length(res_sh)+1]] <- data.table(
        cohorte=co, variable=VARS[[v]], modelo="marginal", n=nrow(d),
        beta=NA_real_, lo=NA_real_, hi=NA_real_, p=NA_real_, estimable=FALSE)
      next
    }
    dd <- copy(d); if (vv != "sexo") dd[[vv]] <- as.numeric(scale(dd[[vv]]))
    fit <- lm(as.formula(paste("shannon ~", vv)), data=dd)
    k <- grep(paste0("^", vv), names(coef(fit)), value=TRUE)[1]
    ci <- confint(fit)[k,]
    res_sh[[length(res_sh)+1]] <- data.table(
      cohorte=co, variable=VARS[[v]], modelo="marginal", n=nrow(d),
      beta=unname(coef(fit)[k]), lo=ci[1], hi=ci[2],
      p=summary(fit)$coefficients[k,4], estimable=TRUE)
  }
  rhs <- c("Age", if (uniqueN(d$sexo) > 1) "sexo", "Weight_kg","Height_cm","BMI","y")
  dd <- copy(d); for (vv in setdiff(rhs,"sexo")) dd[[vv]] <- as.numeric(scale(dd[[vv]]))
  fit <- lm(as.formula(paste("shannon ~", paste(rhs, collapse=" + "))), data=dd)
  ci <- confint(fit)["y",]
  res_sh[[length(res_sh)+1]] <- data.table(
    cohorte=co, variable="Hip BMD (adjusted)", modelo="ajustado", n=nrow(d),
    beta=unname(coef(fit)["y"]), lo=ci[1], hi=ci[2],
    p=summary(fit)$coefficients["y",4], estimable=TRUE)
}
sh <- rbindlist(res_sh)
fwrite(sh, file.path(P$out,"host_variables_shannon.csv"))

cat("\n=== PERMANOVA: R2 marginal de comunidad completa ===\n")
print(dcast(pm, variable ~ cohorte, value.var="R2")[
  match(c(unname(VARS),"Hip BMD (adjusted)"), variable)], row.names=FALSE)
cat("\n=== PERMANOVA: p ===\n")
print(dcast(pm, variable ~ cohorte, value.var="p")[
  match(c(unname(VARS),"Hip BMD (adjusted)"), variable)], row.names=FALSE)
cat("\n=== Shannon: beta estandarizada ===\n")
print(dcast(sh, variable ~ cohorte, value.var="beta")[
  match(c(unname(VARS),"Hip BMD (adjusted)"), variable)], row.names=FALSE)
cat("\n=== Shannon: p ===\n")
print(dcast(sh, variable ~ cohorte, value.var="p")[
  match(c(unname(VARS),"Hip BMD (adjusted)"), variable)], row.names=FALSE)

chk <- pm[cohorte=="Discovery (BGI)" & variable=="Hip BMD" & modelo=="marginal", R2]
if (length(chk) == 1 && abs(chk - 0.006265) > 1e-4)
  warning(sprintf("el R2 marginal del BMD en discovery dio %.6f, no el 0.006265 del script 24", chk))

cat("\nguardado en results/host_variables_permanova.csv y host_variables_shannon.csv\n")
cat("===== FIN OK =====\n")
