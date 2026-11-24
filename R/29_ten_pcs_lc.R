# ======================================================================
# 29_ten_pcs_lc.R — Ten genetic principal components (LC subset)
# #
# Elastic Net models of script 12 with genetic PCs 1-10 instead of 1-2.
# #
# Inputs : as script 12
# Outputs: results/sensitivity_10pcs/lc/
# Paper  : Table S7 (panels B-C)
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
rm(list = ls()); gc()
NPC <- 10; PCV <- paste0('gPC', 1:NPC)
library(data.table); library(readxl)
library(MASS); library(compositions); library(pls)
library(glmnet); library(xgboost)

source("config/paths.R")
P <- list(
  mgs   = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  lcsp  = file.path(BASE,"data/microbiome/los_lc_species_abundance.csv"),
  mast  = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  bmd   = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  clinB = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  clinL = file.path(BASE,"data/phenotypes/los_lc_clinical_imputed.csv"),
  gpca  = file.path(BASE,"data/genotypes/los_qc_pca.eigenvec"),
  disc  = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  lc    = file.path(BASE,"derived/analysis_sets/ids_lc_subset.txt"),
  prs   = file.path(BASE,"derived/prs_los/los_prs.PT_0.01.sscore"),
  prsq  = file.path(BASE,"derived/prs_external/external_prs.PT_0.1.sscore"),
  rds   = file.path(BASE,"derived/mrs/mrs_frozen_crudo.rds"),
  ass   = file.path(BASE,"results/balance_association_crudo.csv"),
  fun   = file.path(BASE,"external/PLS-PBs"),
  out = file.path(BASE,"results/sensitivity_10pcs/lc")
)
for (f in c("fBalChip_PLS.r","fBPMaxOrthNewChip_PLS.R",
            "fBalChipman_PLS.R","fBPUpChi_PLS.r")) source(file.path(P$fun, f))
K <- 5; SEED <- 2026; PREV <- 0.01; PT_MRS <- 0.05; COBERTURA_MIN <- 0.70; B <- 2000
ALPHA <- 0.25
CLIN3 <- c("Age","gender_dummy","Weight_kg")
CLIN  <- c(CLIN3,"Height_cm","Race_African","Race_Caucasian","Race_Asian","Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
           "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")

leer_s <- function(f,nm){ s <- fread(f); setnames(s, sub("^#","",names(s)))
  o <- data.table(IID=trimws(as.character(s$IID)), v=as.numeric(s$SCORE1_AVG)); setnames(o,"v",nm); o }
cerrar <- function(M,ps){ M[M==0] <- ps; (M/rowSums(M))*100 }

disc <- unique(trimws(readLines(P$disc)))
mgs <- fread(file=P$mgs); setnames(mgs,1,"Metagenome_ID")
mgs <- fix_metagenome_ids(mgs)
pteB <- as.data.table(read_excel(P$mast, sheet="bgi"))
setnames(pteB,"BGI_1794","Subject_ID")
pteB <- pteB[, .(Metagenome_ID=trimws(as.character(Metagenome_ID)), IID=trimws(as.character(Subject_ID)))]
dtB <- merge(pteB[IID %in% disc], mgs, by="Metagenome_ID")

b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)),
                        y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
clB <- fread(file=P$clinB); clB[, IID := trimws(as.character(Subject_ID))]
clB <- unique(clB[, c("IID",CLIN), with=FALSE], by="IID")
gp <- fread(file=P$gpca, header=FALSE)
setnames(gp, c("FID","IID", paste0("gPC", seq_len(ncol(gp)-2))))
gp <- gp[, c("IID", PCV), with=FALSE]; gp[, IID := trimws(as.character(IID))]
prsA <- leer_s(P$prs,"PRS"); prsQ <- leer_s(P$prsq,"PRS_qian")

mD <- Reduce(function(x,z) merge(x,z,by="IID"), list(dtB[,.(IID)], ph, clB, gp, prsA, prsQ))
mD <- mD[complete.cases(mD)]
dtB <- dtB[IID %in% mD$IID]; setkey(dtB,"IID"); setkey(mD,"IID"); dtB <- dtB[mD$IID]
stopifnot(identical(dtB$IID, mD$IID))
RAD <- as.matrix(dtB[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RAD) <- dtB$IID
sp_std <- function(v) { s <- sub(".*(s__[^\\.]+).*","\\1",v); s[!grepl("s__",v)] <- NA; s }
colnames(RAD) <- sp_std(colnames(RAD))
RAD <- RAD[, !is.na(colnames(RAD)), drop=FALSE]
yD <- mD$y
cat(sprintf("discovery: %d individuos, %d especies con nombre\n", nrow(RAD), ncol(RAD)))

ajustar_prep <- function(R){
  s <- colnames(R)[apply(R,2,var,na.rm=TRUE)>0]
  s <- s[colMeans(R[,s,drop=FALSE]!=0) >= PREV]
  list(species=s, pseudo=0.5*min(R[,s,drop=FALSE][R[,s,drop=FALSE]!=0]))
}
aplicar_prep <- function(R,pr){
  M <- matrix(0,nrow(R),length(pr$species),dimnames=list(rownames(R),pr$species))
  p <- intersect(pr$species,colnames(R)); M[,p] <- R[,p,drop=FALSE]; cerrar(M,pr$pseudo)
}
pbs <- function(Xc,Bt) as.matrix(compositions::clr(Xc)) %*% Bt

set.seed(SEED)
q <- cut(yD, breaks=quantile(yD,probs=seq(0,1,length.out=K+1)), include.lowest=TRUE)
fold <- integer(length(yD))
for(lv in levels(q)){ i <- which(q==lv); fold[i] <- sample(rep(1:K,length.out=length(i))) }
MRS_D <- numeric(length(yD))
for (k in 1:K) {
  tr <- which(fold!=k); te <- which(fold==k)
  pr <- ajustar_prep(RAD[tr,,drop=FALSE]); Xtr <- aplicar_prep(RAD[tr,,drop=FALSE],pr)
  set.seed(SEED); bal <- fBalChip_PLS(Xtr, yD[tr], version="cov")
  Bt <- t(bal$bal)[pr$species,,drop=FALSE]; PBtr <- pbs(Xtr,Bt)
  pcm <- prcomp(as.matrix(compositions::clr(RAD[tr,,drop=FALSE]+1e-6)), center=TRUE, scale.=FALSE)
  cov <- data.table(y=yD[tr], mD[tr,CLIN3,with=FALSE], mPC1=pcm$x[,1], mPC2=pcm$x[,2])
  ss <- rbindlist(lapply(colnames(PBtr), function(nm){
    co <- summary(lm(y ~ pb + Age + gender_dummy + Weight_kg + mPC1 + mPC2,
                     data=cbind(cov, pb=PBtr[,nm])))$coefficients
    if(!"pb" %in% rownames(co)) return(data.table(PB_ID=nm,beta=NA_real_,p=NA_real_))
    data.table(PB_ID=nm, beta=co["pb","Estimate"], p=co["pb","Pr(>|t|)"]) }))
  s <- ss[!is.na(p) & p < PT_MRS]
  PBte <- pbs(aplicar_prep(RAD[te,,drop=FALSE],pr), Bt)
  MRS_D[fold==k] <- if (nrow(s)) as.numeric(PBte[, s$PB_ID, drop=FALSE] %*% s$beta) else 0
}
cat("MRS out-of-fold del discovery: listo\n")

lc_id <- unique(trimws(readLines(P$lc)))
pteL <- as.data.table(read_excel(P$mast, sheet="lc"))
setnames(pteL, "LC_495", "IID")
pteL <- pteL[, .(IID=trimws(as.character(IID)), Metagenome_ID=trimws(as.character(Metagenome_ID)))]

tab <- fread(file=P$lcsp)
sp_lc <- trimws(as.character(tab[[1]]))
M <- as.matrix(tab[, -1, with=FALSE])
rownames(M) <- sp_lc
RAL <- t(M)
RAL <- RAL[rownames(RAL) %in% pteL$Metagenome_ID, , drop=FALSE]
idx <- match(rownames(RAL), pteL$Metagenome_ID)
rownames(RAL) <- pteL$IID[idx]
cat(sprintf("LC: %d muestras con metagenoma, %d especies en la tabla\n", nrow(RAL), ncol(RAL)))

clL <- fread(file=P$clinL); clL[, IID := trimws(as.character(Subject_ID))]
clL <- unique(clL[, c("IID",CLIN), with=FALSE], by="IID")
mL <- Reduce(function(x,z) merge(x,z,by="IID"),
             list(data.table(IID=intersect(rownames(RAL), lc_id)), ph, clL, gp, prsA, prsQ))
mL <- mL[complete.cases(mL)]
RAL <- RAL[mL$IID, , drop=FALSE]
yL <- mL$y
cat(sprintf("LC con datos completos: %d de %d\n\n", nrow(mL), length(lc_id)))

fr  <- readRDS(P$rds)
rn <- sp_std(rownames(fr$Bt))
stopifnot(!any(is.na(rn)), !any(duplicated(rn)))
rownames(fr$Bt) <- rn
fr$prep$species <- sp_std(fr$prep$species)
names(fr$prep$pseudo) <- NULL
ass <- fread(file=P$ass)
sel <- ass[!is.na(p) & p < PT_MRS]
Bt_full <- fr$Bt[, sel$PB_ID, drop=FALSE]
esp_disc <- rownames(Bt_full)
presentes <- intersect(esp_disc, colnames(RAL))
cat(sprintf("especies del MRS presentes en LC: %d de %d (%.1f%%)\n",
            length(presentes), length(esp_disc), 100*length(presentes)/length(esp_disc)))
if (length(presentes) < 0.25*length(esp_disc))
  stop(sprintf(paste("Solo %d de %d especies cruzan. Esto es un problema de",
                     "formato de nombres, no un resultado biologico.",
                     "Comparar: rownames(fr$Bt)[1] y colnames(RAL)[1]."),
               length(presentes), length(esp_disc)))

cob <- rbindlist(lapply(colnames(Bt_full), function(nm) {
  bfull <- Bt_full[, nm]
  r0 <- sum(bfull > 0); s0 <- sum(bfull < 0)
  bp <- bfull[presentes]
  r1 <- sum(bp > 0); s1 <- sum(bp < 0)
  data.table(PB_ID=nm, n_especies=r0+s0, n_presentes=r1+s1,
             cobertura=round((r1+s1)/(r0+s0),3), num=r1, den=s1,
             evaluable=((r1+s1)/(r0+s0) >= COBERTURA_MIN) & r1>0 & s1>0)
}))
cat(sprintf("\nbalances seleccionados: %d | evaluables en LC (cobertura >= %.0f%%): %d\n",
            nrow(cob), 100*COBERTURA_MIN, sum(cob$evaluable)))
print(cob[order(-cobertura)], row.names=FALSE)
fwrite(cob, file.path(P$out,"lc_balance_coverage.csv"))

if (sum(cob$evaluable) == 0) {
  cat("\n*** NINGUN balance es evaluable en LC. El MRS no se puede transferir. ***\n")
  MRS_L <- rep(0, nrow(RAL))
} else {
  XL <- RAL[, presentes, drop=FALSE]
  XL <- cerrar(XL, fr$prep$pseudo)
  ZL <- as.matrix(compositions::clr(XL))
  ev <- cob[evaluable==TRUE, PB_ID]
  Bt_re <- sapply(ev, function(nm) {
    bp <- Bt_full[presentes, nm]
    r1 <- sum(bp>0); s1 <- sum(bp<0); v <- numeric(length(bp))
    v[bp>0] <-  sqrt(s1/(r1*(r1+s1))); v[bp<0] <- -sqrt(r1/(s1*(r1+s1))); v
  })
  rownames(Bt_re) <- presentes
  MRS_L <- as.numeric((ZL %*% Bt_re) %*% sel[PB_ID %in% ev, beta])
}

z <- function(v) (v-mean(v))/sd(v)

vars <- c(CLIN,PCV)
mAll <- rbind(mD[, vars, with=FALSE], mL[, vars, with=FALSE])
for (cn in vars) if (is.character(mAll[[cn]])) set(mAll, j=cn, value=factor(mAll[[cn]]))

const_LC   <- vars[sapply(vars, function(cn) uniqueN(mL[[cn]]) < 2)]
const_disc <- vars[sapply(vars, function(cn) uniqueN(mD[[cn]]) < 2)]
if (length(const_LC))
  cat("constantes dentro de LC (se conservan si varian en discovery):",
      paste(const_LC, collapse=", "), "\n")
quitar <- vars[sapply(vars, function(cn) uniqueN(mAll[[cn]]) < 2)]
if (length(quitar)) {
  cat("quitadas por ser constantes en AMBAS cohortes:", paste(quitar, collapse=", "), "\n")
  mAll[, (quitar) := NULL]
}

XA <- model.matrix(~ ., data=as.data.frame(mAll))[,-1,drop=FALSE]
colnames(XA) <- make.names(colnames(XA), unique=TRUE)
XA <- XA[, apply(XA,2,function(w) length(unique(w))>1), drop=FALSE]
XD <- XA[seq_len(nrow(mD)), , drop=FALSE]
XL <- XA[nrow(mD) + seq_len(nrow(mL)), , drop=FALSE]
stopifnot(identical(colnames(XD), colnames(XL)), nrow(XL) == nrow(mL))
cat(sprintf("\npredictores en la matriz comun: %d\n", ncol(XA)))
vclin <- setdiff(colnames(XA), PCV)

DD <- cbind(as.data.table(XD),
            data.table(PRS=z(mD$PRS), PRS_qian=z(mD$PRS_qian), MRS=z(MRS_D)))
DL <- cbind(as.data.table(XL),
            data.table(PRS=z(mL$PRS), PRS_qian=z(mL$PRS_qian),
                       MRS=if (sd(MRS_L)>0) z(MRS_L) else rep(0, length(MRS_L))))

asoc_cohorte <- function(Dx, yx, covs, etiqueta) {
  covs <- covs[sapply(covs, function(cn) uniqueN(Dx[[cn]]) > 1)]
  rbindlist(lapply(c("PRS","PRS_qian","MRS"), function(s) {
    if (uniqueN(Dx[[s]]) < 2)
      return(data.table(score=s, conjunto=etiqueta, n=length(yx), beta=NA_real_, se=NA_real_,
                        IC95_inf=NA_real_, IC95_sup=NA_real_, P=NA_real_))
    d <- data.table(y=yx, Dx[, c(covs, s), with=FALSE])
    fit <- lm(as.formula(paste("y ~", paste(c(covs, s), collapse="+"))), data=d)
    co <- summary(fit)$coefficients; ci <- confint(fit, s, level=0.95)
    data.table(score=s, conjunto=etiqueta, n=length(yx),
               beta=round(co[s,"Estimate"],5), se=round(co[s,"Std. Error"],5),
               IC95_inf=round(ci[1],5), IC95_sup=round(ci[2],5),
               P=signif(co[s,"Pr(>|t|)"],3))
  }))
}
assocL <- asoc_cohorte(DL, yL, c(vclin,PCV), "LC")
cat("\n=== LC: asociacion de cada score con HTOT_BMD (por DE del score) ===\n")
print(assocL, row.names=FALSE)
fwrite(assocL, file.path(P$out,"lc_score_association.csv"))

ESPEC <- list(
  S1_clinicas          = vclin,
  S2_base              = c(vclin,PCV),
  S3_base_PRS          = c(vclin,PCV,"PRS"),
  S3q_base_PRSqian     = c(vclin,PCV,"PRS_qian"),
  S4_base_MRS          = c(vclin,PCV,"MRS"),
  S5_base_PRS_MRS      = c(vclin,PCV,"PRS","MRS"),
  S5q_base_PRSqian_MRS = c(vclin,PCV,"PRS_qian","MRS")
)
fit_enet <- function(Xtr,ytr,Xte){ set.seed(SEED)
  cv <- cv.glmnet(Xtr,ytr,alpha=ALPHA,nfolds=10,standardize=TRUE)
  as.numeric(predict(cv,newx=Xte,s="lambda.min")) }
fit_xgb <- function(Xtr,ytr,Xte){
  dtr <- xgb.DMatrix(Xtr,label=ytr)
  pars <- list(objective="reg:squarederror", eta=0.05, max_depth=3,
               subsample=0.8, colsample_bytree=0.8, min_child_weight=5, nthread=2)
  set.seed(SEED)
  cv <- xgb.cv(pars,dtr,nrounds=400,nfold=5,early_stopping_rounds=30,verbose=0)
  nb <- tryCatch({ lg <- as.data.frame(cv$evaluation_log)
                   which.min(lg[[grep("^test.*mean$",names(lg),value=TRUE)[1]]]) },
                 error=function(e) NA_integer_)
  if (length(nb)!=1 || is.na(nb) || nb<1) nb <- 100
  set.seed(SEED); m <- xgb.train(pars,dtr,nrounds=nb,verbose=0)
  predict(m, xgb.DMatrix(Xte)) }

r2 <- function(o,p) 1 - sum((o-p)^2)/sum((o-mean(o))^2)
cat("\n>>> APLICANDO LOS MODELOS CONGELADOS A LC <<<\n\n")
PREDL <- data.table(IID=mL$IID, y=yL)
res <- data.table()
for (nm in names(ESPEC)) {
  v <- ESPEC[[nm]]
  Xtr <- as.matrix(DD[, v, with=FALSE]); Xte <- as.matrix(DL[, v, with=FALSE])
  pe <- fit_enet(Xtr,yD,Xte); px <- rep(NA_real_, nrow(Xte))
  PREDL[[paste0(nm,"__ElasticNet")]] <- pe
  PREDL[[paste0(nm,"__XGBoost")]]    <- px
  res <- rbind(res, data.table(especificacion=nm, algoritmo="ElasticNet", R2_LC=round(r2(yL,pe),4)),
                    data.table(especificacion=nm, algoritmo="XGBoost",    R2_LC=round(r2(yL,px),4)))
  cat(sprintf("  %-22s EN %.4f | XGB %.4f\n", nm, r2(yL,pe), r2(yL,px)))
}
base <- res[especificacion=="S2_base", .(algoritmo, R2_base=R2_LC)]
res <- merge(res, base, by="algoritmo"); res[, dR2_vs_S2 := round(R2_LC-R2_base,5)]
setorder(res, algoritmo, especificacion)
cat(sprintf("\n=== LC (n = %d) ===\n", nrow(mL)))
print(res[, .(algoritmo, especificacion, R2_LC, dR2_vs_S2)], row.names=FALSE)
fwrite(res, file.path(P$out,"lc_results.csv"))
fwrite(PREDL, file.path(P$out,"lc_predictions.csv"))

set.seed(SEED); nL <- nrow(PREDL)
idxb <- replicate(B, sample.int(nL,nL,replace=TRUE), simplify=FALSE)
cols <- setdiff(names(PREDL), c("IID","y"))
R2B <- sapply(cols, function(cn){ p <- PREDL[[cn]]
  vapply(idxb, function(i) r2(yL[i],p[i]), numeric(1)) })
ic <- rbindlist(lapply(c("ElasticNet","XGBoost"), function(a)
  rbindlist(lapply(names(ESPEC), function(nm){
    c1 <- paste0(nm,"__",a); c0 <- paste0("S2_base__",a)
    d <- R2B[,c1]-R2B[,c0]; q <- quantile(d,c(0.025,0.975), na.rm=TRUE)
    data.table(algoritmo=a, especificacion=nm,
               dR2=round(r2(yL,PREDL[[c1]])-r2(yL,PREDL[[c0]]),4),
               IC95_inf=round(q[1],4), IC95_sup=round(q[2],4),
               cruza_cero=(q[1]<=0 & q[2]>=0)) }))))
cat("\n=== LC con intervalos (bootstrap B = 2000) ===\n")
print(ic, row.names=FALSE)
fwrite(ic, file.path(P$out,"lc_ci.csv"))
cat("\nguardado en results/sensitivity_10pcs/lc/\n")
cat("===== FIN OK =====\n")
