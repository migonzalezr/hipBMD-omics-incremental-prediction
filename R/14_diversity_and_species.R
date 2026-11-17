# ======================================================================
# 14_diversity_and_species.R — Diversity and species-level associations
#
# Shannon diversity, Bray-Curtis PERMANOVA and species-level associations with hip BMD before and after adjustment for age, sex and weight.
#
# Inputs : species abundances of the three cohorts
# Outputs: results/diversity_*.csv; species_association.csv
# Paper  : Figure 3; Tables S9-S10
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl); library(vegan); library(compositions)

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

sp428 <- colnames(RAD)[apply(RAD,2,var,na.rm=TRUE) > 0]
sp428 <- sp428[colMeans(RAD[,sp428,drop=FALSE] != 0) >= PREV]
cat(sprintf("UNIVERSO PRE-DECLARADO: %d especies (prevalencia >= %.0f%% en el discovery)\n",
            length(sp428), 100*PREV))

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
clC[, `:=`(Age=as.numeric(Age), gender_dummy="Female",
           Weight_kg=as.numeric(Weight_kg), y=as.numeric(HTOT_BMD))]
mC <- clC[IID %in% rownames(RAC), c("IID","y",CLIN), with=FALSE]
mC <- mC[complete.cases(mC)]; RAC <- RAC[mC$IID, , drop=FALSE]

cat(sprintf("discovery %d | LC %d | Guangzhou %d individuos\n", nrow(RAD), nrow(RAL), nrow(RAC)))

enLC <- intersect(sp428, colnames(RAL)); enCH <- intersect(sp428, colnames(RAC))
venn <- data.table(
  conjunto = c("discovery (universo)","presentes en LC","presentes en Guangzhou",
               "en las tres","solo discovery"),
  n = c(length(sp428), length(enLC), length(enCH),
        length(intersect(enLC, enCH)),
        length(setdiff(sp428, union(enLC, enCH)))))
venn[, pct_del_universo := round(100*n/length(sp428), 1)]
cat("\n=== solape de especies ===\n"); print(venn, row.names=FALSE)
fwrite(venn, file.path(P$out,"diversity_venn.csv"))

fwrite(data.table(especie=sp428, en_LC=sp428 %in% enLC, en_China=sp428 %in% enCH),
       file.path(P$out,"diversity_species_universe.csv"))

comunes <- Reduce(intersect, list(sp428, colnames(RAL), colnames(RAC)))
cat(sprintf("\nespecies comunes a las tres cohortes: %d\n", length(comunes)))
renorm <- function(M, sp){ X <- M[, sp, drop=FALSE]; X[X<0] <- 0; X/pmax(rowSums(X), .Machine$double.eps) }
sh <- rbindlist(list(
  data.table(IID=rownames(RAD), cohorte="Discovery (BGI)",
             shannon=vegan::diversity(renorm(RAD,comunes), index="shannon"),
             sexo=mD$gender_dummy, y=mD$y, Age=mD$Age, Weight_kg=mD$Weight_kg),
  data.table(IID=rownames(RAL), cohorte="Validation 1 (LC)",
             shannon=vegan::diversity(renorm(RAL,comunes), index="shannon"),
             sexo=mL$gender_dummy, y=mL$y, Age=mL$Age, Weight_kg=mL$Weight_kg),
  data.table(IID=rownames(RAC), cohorte="Validation 2 (Guangzhou)",
             shannon=vegan::diversity(renorm(RAC,comunes), index="shannon"),
             sexo=mC$gender_dummy, y=mC$y, Age=mC$Age, Weight_kg=mC$Weight_kg)))
fwrite(sh, file.path(P$out,"diversity_shannon.csv"))
cat("\n=== Shannon por cohorte ===\n")
print(sh[, .(n=.N, media=round(mean(shannon),3), sd=round(sd(shannon),3)), by=cohorte], row.names=FALSE)

w <- wilcox.test(shannon ~ sexo, data=sh[cohorte=="Discovery (BGI)"])
cat(sprintf("Shannon por sexo en el discovery: Wilcoxon p = %.3g\n", w$p.value))

asoc <- rbindlist(lapply(unique(sh$cohorte), function(co){
  d <- sh[cohorte==co]
  uno <- function(f, et){
    m <- summary(lm(as.formula(f), data=d))$coefficients
    data.table(cohorte=co, modelo=et, n=nrow(d),
               beta=signif(m["shannon","Estimate"],3),
               se=signif(m["shannon","Std. Error"],3),
               p=signif(m["shannon","Pr(>|t|)"],3))
  }
  rhs <- if (uniqueN(d$sexo) > 1) "+ Age + sexo + Weight_kg" else "+ Age + Weight_kg"
  rbind(uno("y ~ shannon", "sin ajuste"),
        uno(paste("y ~ shannon", rhs), "ajustado por edad, sexo y peso"))
}))
cat("\n=== Shannon vs HTOT_BMD ===\n"); print(asoc, row.names=FALSE)
fwrite(asoc, file.path(P$out,"diversity_shannon_association.csv"))

Xall <- rbind(renorm(RAD,comunes), renorm(RAL,comunes), renorm(RAC,comunes))
meta_all <- data.table(IID=rownames(Xall),
                       cohorte=rep(c("Discovery (BGI)","Validation 1 (LC)","Validation 2 (Guangzhou)"),
                                   c(nrow(RAD),nrow(RAL),nrow(RAC))))
cat(sprintf("\nBray-Curtis sobre %d muestras x %d especies...\n", nrow(Xall), ncol(Xall)))
D <- vegan::vegdist(Xall, method="bray")
pc <- cmdscale(D, k=2, eig=TRUE)
ev <- round(100*pc$eig[1:2]/sum(pc$eig[pc$eig>0]), 1)
pcoa <- data.table(meta_all, PCo1=pc$points[,1], PCo2=pc$points[,2],
                   var_PCo1=ev[1], var_PCo2=ev[2])
fwrite(pcoa, file.path(P$out,"diversity_pcoa.csv"))
cat(sprintf("PCoA: PCo1 %.1f%%, PCo2 %.1f%% de la varianza\n", ev[1], ev[2]))

set.seed(SEED)
pm_coh <- vegan::adonis2(D ~ cohorte, data=meta_all, permutations=NPERM)
cat("\n=== PERMANOVA: diferencias entre cohortes ===\n"); print(pm_coh)

Dd <- vegan::vegdist(renorm(RAD,comunes), method="bray")
md <- data.table(y=mD$y, Age=mD$Age, sexo=mD$gender_dummy, Weight_kg=mD$Weight_kg)
set.seed(SEED); pm_cru <- vegan::adonis2(Dd ~ y, data=md, permutations=NPERM)
set.seed(SEED); pm_aju <- vegan::adonis2(Dd ~ Age + sexo + Weight_kg + y, data=md,
                                         permutations=NPERM, by="margin")
cat("\n=== PERMANOVA en el discovery: BMD sin ajustar ===\n"); print(pm_cru)
cat("\n=== PERMANOVA en el discovery: BMD ajustado (terminos marginales) ===\n"); print(pm_aju)

saca <- function(x, et){
  d <- as.data.table(as.data.frame(x), keep.rownames="termino")
  setnames(d, c("R2","Pr(>F)"), c("R2","p"), skip_absent=TRUE)
  d[, analisis := et][]
}
fwrite(rbindlist(list(saca(pm_coh,"cohortes"), saca(pm_cru,"BMD sin ajustar"),
                      saca(pm_aju,"BMD ajustado")), fill=TRUE),
       file.path(P$out,"diversity_permanova.csv"))

Xsp <- RAD[, sp428, drop=FALSE]
pseudo <- 0.5*min(Xsp[Xsp!=0]); Xsp[Xsp==0] <- pseudo; Xsp <- (Xsp/rowSums(Xsp))*100
Z <- as.matrix(compositions::clr(Xsp))
cov_sp <- data.table(y=mD$y, Age=mD$Age, sexo=mD$gender_dummy, Weight_kg=mD$Weight_kg)
esp <- rbindlist(lapply(colnames(Z), function(nm){
  d <- cbind(cov_sp, sp=Z[,nm])
  m0 <- summary(lm(y ~ sp, data=d))$coefficients
  m1 <- summary(lm(y ~ sp + Age + sexo + Weight_kg, data=d))$coefficients
  data.table(especie=nm,
             beta_crudo=m0["sp","Estimate"], p_crudo=m0["sp","Pr(>|t|)"],
             beta_ajust=m1["sp","Estimate"], p_ajust=m1["sp","Pr(>|t|)"])
}))
esp[, `:=`(FDR_crudo=p.adjust(p_crudo, "BH"), FDR_ajust=p.adjust(p_ajust, "BH"))]
setorder(esp, p_ajust)
cat(sprintf("\n=== asociacion por especie, universo de %d ===\n", nrow(esp)))
cat(sprintf("  sin ajustar : %d especies con FDR<0.05 (FDR minimo %.3g)\n",
            sum(esp$FDR_crudo<0.05), min(esp$FDR_crudo)))
cat(sprintf("  ajustado    : %d especies con FDR<0.05 (FDR minimo %.3g)\n",
            sum(esp$FDR_ajust<0.05), min(esp$FDR_ajust)))
cat("\ntop 10 por p ajustada:\n")
print(head(esp[, .(especie, beta_crudo=signif(beta_crudo,3), FDR_crudo=signif(FDR_crudo,3),
                   beta_ajust=signif(beta_ajust,3), FDR_ajust=signif(FDR_ajust,3))], 10),
      row.names=FALSE)

sa <- esp[grepl("australis", especie)]
if (nrow(sa)) { cat("\nStreptococcus australis, el que reportaba el manuscrito:\n"); print(sa, row.names=FALSE) }
fwrite(esp, file.path(P$out,"species_association.csv"))

cat("\nguardado en results/\n")
cat("===== FIN OK =====\n")

