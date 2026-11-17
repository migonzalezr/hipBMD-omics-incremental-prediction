# ======================================================================
# 13_transfer_guangzhou.R — Transfer to the Guangzhou cohort
#
# Applies the frozen scores and the reduced six-predictor base model to the Guangzhou cohort; reports calibration and recalibrated performance.
#
# Inputs : data/phenotypes/guangzhou_clinical.csv; data/microbiome/guangzhou_species_abundance.csv; derived/prs_guangzhou/*
# Outputs: results/guangzhou_predictions.csv; guangzhou_results.csv; guangzhou_calibration.csv
# Paper  : Table 2; Figure 4
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)
library(MASS); library(compositions); library(pls)
library(glmnet); library(xgboost)

source("config/paths.R")
P <- list(
  mgs   = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  chsp  = file.path(BASE,"data/microbiome/guangzhou_species_abundance.csv"),
  chcl  = file.path(BASE,"data/phenotypes/guangzhou_clinical.csv"),
  mast  = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  bmd   = file.path(BASE,"data/phenotypes/los_hip_bmd.xlsx"),
  clinB = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  disc  = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  prsC  = file.path(BASE,"derived/prs_guangzhou/guangzhou_prs_los_weights.txt"),
  prsCq = file.path(BASE,"derived/prs_guangzhou/guangzhou_prs_external_weights.txt"),
  prs   = file.path(BASE,"derived/prs_los/los_prs.PT_0.01.sscore"),
  prsq  = file.path(BASE,"derived/prs_external/external_prs.PT_0.1.sscore"),
  rds   = file.path(BASE,"derived/mrs/mrs_frozen_crudo.rds"),
  ass   = file.path(BASE,"results/balance_association_crudo.csv"),
  fun   = file.path(BASE,"external/PLS-PBs"),
  out   = file.path(BASE,"results")
)
falta <- unlist(P)[!file.exists(unlist(P))]
if (length(falta)) stop("no encuentro:\n  ", paste(falta, collapse="\n  "))
for (f in c("fBalChip_PLS.r","fBPMaxOrthNewChip_PLS.R",
            "fBalChipman_PLS.R","fBPUpChi_PLS.r")) source(file.path(P$fun, f))

K <- 5; SEED <- 2026; PREV <- 0.01; PT_MRS <- 0.05; COBERTURA_MIN <- 0.70; B <- 2000
ALPHA <- 0.25
CLIN3 <- c("Age","gender_dummy","Weight_kg")
CLIN  <- c(CLIN3,"Height_cm","BMI","Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
           "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")

RED <- c("Age","Weight_kg","Height_cm","BMI","Alcohol","Exercise")

leer_s  <- function(f,nm){ s <- fread(file=f); setnames(s, sub("^#","",names(s)))
  o <- data.table(IID=trimws(as.character(s$IID)), v=as.numeric(s$SCORE1_AVG)); setnames(o,"v",nm); o }
leer_ch <- function(f,nm){ s <- fread(file=f)
  o <- data.table(IID=trimws(as.character(s$IID)), v=as.numeric(s$SCORE_AVG)); setnames(o,"v",nm); o }
cerrar  <- function(M,ps){ M[M==0] <- ps; (M/rowSums(M))*100 }
sp_std  <- function(v){ s <- sub(".*(s__[^\\.]+).*","\\1",v); s[!grepl("s__",v)] <- NA; s }
r2  <- function(o,p) 1 - sum((o-p)^2)/sum((o-mean(o))^2)
rmse<- function(o,p) sqrt(mean((o-p)^2))
z   <- function(v) (v-mean(v))/sd(v)

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
prsA <- leer_s(P$prs,"PRS"); prsQ <- leer_s(P$prsq,"PRS_qian")

mD <- Reduce(function(x,zz) merge(x,zz,by="IID"), list(dtB[,.(IID)], ph, clB, prsA, prsQ))
mD <- mD[complete.cases(mD)]
dtB <- dtB[IID %in% mD$IID]; setkey(dtB,"IID"); setkey(mD,"IID"); dtB <- dtB[mD$IID]
stopifnot(identical(dtB$IID, mD$IID))
RAD <- as.matrix(dtB[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RAD) <- dtB$IID
colnames(RAD) <- sp_std(colnames(RAD)); RAD <- RAD[, !is.na(colnames(RAD)), drop=FALSE]
yD <- mD$y
cat(sprintf("discovery: %d individuos, %d especies\n", nrow(RAD), ncol(RAD)))

si <- function(v) as.integer(trimws(as.character(v)) %in% c("Yes","yes","Y","1","TRUE"))
mD[, `:=`(Alcohol  = si(Alcohol_Drinking_cat),
          Exercise = si(Exercise_Regular_cat),
          Milk     = si(Milk_Drink_cat),
          mujer    = as.integer(trimws(as.character(gender_dummy)) %in% c("Female","female","F")))]

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

clC <- fread(file=P$chcl); clC[, IID := trimws(as.character(ID))]
clC[, `:=`(Alcohol  = as.integer(Alcohol_cat),
           Exercise = as.integer(Exercise_cat),
           Milk     = as.integer(Milk_cat),
           y        = as.numeric(HTOT_BMD))]
prsC <- leer_ch(P$prsC,"PRS"); prsCq <- leer_ch(P$prsCq,"PRS_qian")

tabC <- fread(file=P$chsp)
spC  <- sp_std(trimws(as.character(tabC[[1]])))
MC   <- as.matrix(tabC[, -1, with=FALSE]); rownames(MC) <- as.character(tabC[[1]])
RAC  <- t(MC); colnames(RAC) <- spC
RAC  <- RAC[, !is.na(colnames(RAC)), drop=FALSE]
RAC  <- RAC[, !duplicated(colnames(RAC)), drop=FALSE]
cat(sprintf("Guangzhou: %d muestras, %d especies con nombre de especie\n", nrow(RAC), ncol(RAC)))

cas <- data.table(paso=character(), n=integer())
cas <- rbind(cas, data.table(paso="en el archivo clinico", n=nrow(clC)))
con_sp  <- intersect(clC$IID, rownames(RAC))
cas <- rbind(cas, data.table(paso="con tabla de especies", n=length(con_sp)))
con_p   <- intersect(con_sp, prsC$IID)
cas <- rbind(cas, data.table(paso="con PRS de LOS",  n=length(con_p)))
con_pq  <- intersect(con_p,  prsCq$IID)
cas <- rbind(cas, data.table(paso="con PRS de Qian", n=length(con_pq)))
tmp <- clC[IID %in% con_pq, c("IID","y",RED), with=FALSE]
for (cn in c("y",RED)) {
  nm <- sum(is.na(tmp[[cn]]))
  if (nm > 0) cat(sprintf("  AVISO: %s tiene %d ausentes entre los %d cruzados\n", cn, nm, nrow(tmp)))
}
cas <- rbind(cas, data.table(paso="sin ausentes en las covariables", n=sum(complete.cases(tmp))))
cat("\n=== de donde sale cada perdida ===\n"); print(cas, row.names=FALSE)

mC <- Reduce(function(x,zz) merge(x,zz,by="IID"),
             list(data.table(IID=con_pq),
                  clC[, c("IID","y",RED), with=FALSE], prsC, prsCq))
mC <- mC[complete.cases(mC)]
RAC <- RAC[mC$IID, , drop=FALSE]
yC <- mC$y
cat(sprintf("Guangzhou con datos completos: %d\n", nrow(mC)))
if (nrow(mC) < 450)
  stop(sprintf(paste("solo quedan %d de %d. Mirar la tabla 'de donde sale cada perdida'",
                     "de arriba: dice en que paso se van."), nrow(mC), nrow(clC)))

comp <- data.table(
  variable = c("n","edad media","edad sd","% mujeres","peso medio","BMI medio",
               "BMD media","BMD sd"),
  discovery = c(nrow(mD), round(mean(mD$Age),1), round(sd(mD$Age),1),
                round(100*mean(mD$mujer),1), round(mean(mD$Weight_kg),1),
                round(mean(mD$BMI),1), round(mean(yD),4), round(sd(yD),4)),
  guangzhou = c(nrow(mC), round(mean(mC$Age),1), round(sd(mC$Age),1), 100,
                round(mean(mC$Weight_kg),1), round(mean(mC$BMI),1),
                round(mean(yC),4), round(sd(yC),4)))
cat("\n=== composicion de las cohortes ===\n"); print(comp, row.names=FALSE)

fr <- readRDS(P$rds)
rn <- sp_std(rownames(fr$Bt))
stopifnot(!any(is.na(rn)), !any(duplicated(rn)))
rownames(fr$Bt) <- rn
fr$prep$species <- sp_std(fr$prep$species); names(fr$prep$pseudo) <- NULL

ass <- fread(file=P$ass); sel <- ass[!is.na(p) & p < PT_MRS]
Bt_full <- fr$Bt[, sel$PB_ID, drop=FALSE]
esp_disc <- rownames(Bt_full)
presentes <- intersect(esp_disc, colnames(RAC))
cat(sprintf("\nespecies del MRS presentes en Guangzhou: %d de %d (%.1f%%)\n",
            length(presentes), length(esp_disc), 100*length(presentes)/length(esp_disc)))
if (length(presentes) < 0.25*length(esp_disc))
  stop(sprintf(paste("Solo %d de %d especies cruzan. Es un problema de formato de",
                     "nombres, no un resultado. Comparar rownames(fr$Bt)[1] con colnames(RAC)[1]."),
               length(presentes), length(esp_disc)))

cob <- rbindlist(lapply(colnames(Bt_full), function(nm){
  bf <- Bt_full[, nm]; r0 <- sum(bf>0); s0 <- sum(bf<0)
  bp <- bf[presentes]; r1 <- sum(bp>0); s1 <- sum(bp<0)
  data.table(PB_ID=nm, n_especies=r0+s0, n_presentes=r1+s1,
             cobertura=round((r1+s1)/(r0+s0),3), num=r1, den=s1,
             evaluable=((r1+s1)/(r0+s0) >= COBERTURA_MIN) & r1>0 & s1>0) }))
cat(sprintf("balances seleccionados: %d | evaluables en Guangzhou: %d\n",
            nrow(cob), sum(cob$evaluable)))
print(cob[order(-cobertura)], row.names=FALSE)
fwrite(cob, file.path(P$out,"guangzhou_balance_coverage.csv"))

if (sum(cob$evaluable)==0) {
  cat("\n*** NINGUN balance evaluable en Guangzhou. El MRS no se transfiere. ***\n")
  MRS_C <- rep(0, nrow(RAC))
} else {
  XC <- cerrar(RAC[, presentes, drop=FALSE], fr$prep$pseudo)
  ZC <- as.matrix(compositions::clr(XC))
  ev <- cob[evaluable==TRUE, PB_ID]
  Bt_re <- sapply(ev, function(nm){
    bp <- Bt_full[presentes, nm]; r1 <- sum(bp>0); s1 <- sum(bp<0); v <- numeric(length(bp))
    v[bp>0] <-  sqrt(s1/(r1*(r1+s1))); v[bp<0] <- -sqrt(r1/(s1*(r1+s1))); v })
  rownames(Bt_re) <- presentes
  MRS_C <- as.numeric((ZC %*% Bt_re) %*% sel[PB_ID %in% ev, beta])
}

mAll <- rbind(mD[, RED, with=FALSE], mC[, RED, with=FALSE])
XA <- model.matrix(~ ., data=as.data.frame(mAll))[,-1,drop=FALSE]
colnames(XA) <- make.names(colnames(XA), unique=TRUE)
XA <- XA[, apply(XA,2,function(w) length(unique(w))>1), drop=FALSE]
XD <- XA[seq_len(nrow(mD)), , drop=FALSE]
XC2<- XA[nrow(mD) + seq_len(nrow(mC)), , drop=FALSE]
stopifnot(identical(colnames(XD), colnames(XC2)), nrow(XC2)==nrow(mC))
cat(sprintf("\npredictores del modelo reducido: %d (%s)\n",
            ncol(XA), paste(colnames(XA), collapse=", ")))

DD <- cbind(as.data.table(XD), data.table(PRS=z(mD$PRS), PRS_qian=z(mD$PRS_qian), MRS=z(MRS_D)))
DC <- cbind(as.data.table(XC2),
            data.table(PRS=z(mC$PRS), PRS_qian=z(mC$PRS_qian),
                       MRS=if (sd(MRS_C)>0) z(MRS_C) else rep(0,length(MRS_C))))
vclin <- colnames(XA)

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
stopifnot(length(yC) == nrow(DC))
assocC <- asoc_cohorte(DC, yC, vclin, "Guangzhou")
cat("\n=== R1.3 en Guangzhou: asociacion de cada score con HTOT_BMD (por DE del score) ===\n")
print(assocC, row.names=FALSE)
fwrite(assocC, file.path(P$out,"guangzhou_score_association.csv"))

ESPEC <- list(
  S2_base              = vclin,
  S3_base_PRS          = c(vclin,"PRS"),
  S3q_base_PRSqian     = c(vclin,"PRS_qian"),
  S4_base_MRS          = c(vclin,"MRS"),
  S5_base_PRS_MRS      = c(vclin,"PRS","MRS"),
  S5q_base_PRSqian_MRS = c(vclin,"PRS_qian","MRS")
)
fit_enet <- function(Xtr,ytr,Xte){ set.seed(SEED)
  cv <- cv.glmnet(Xtr,ytr,alpha=ALPHA,nfolds=10,standardize=TRUE)
  as.numeric(predict(cv,newx=Xte,s="lambda.min")) }
fit_xgb <- function(Xtr,ytr,Xte){
  dtr <- xgb.DMatrix(Xtr,label=ytr)
  pars <- list(objective="reg:squarederror", eta=0.05, max_depth=3,
               subsample=0.8, colsample_bytree=0.8, min_child_weight=5)
  set.seed(SEED)
  cv <- xgb.cv(pars,dtr,nrounds=400,nfold=5,early_stopping_rounds=30,verbose=0)
  nb <- tryCatch({ lg <- as.data.frame(cv$evaluation_log)
                   which.min(lg[[grep("^test.*mean$",names(lg),value=TRUE)[1]]]) },
                 error=function(e) NA_integer_)
  if (length(nb)!=1 || is.na(nb) || nb<1) nb <- 100
  set.seed(SEED); m <- xgb.train(pars,dtr,nrounds=nb,verbose=0)
  predict(m, xgb.DMatrix(Xte)) }

cat("\n>>> APLICANDO LOS MODELOS REDUCIDOS A GUANGZHOU <<<\n\n")
PREDC <- data.table(IID=mC$IID, y=yC)
res <- data.table()
for (nm in names(ESPEC)) {
  v <- ESPEC[[nm]]
  Xtr <- as.matrix(DD[, v, with=FALSE]); Xte <- as.matrix(DC[, v, with=FALSE])
  pe <- fit_enet(Xtr,yD,Xte); px <- fit_xgb(Xtr,yD,Xte)
  PREDC[[paste0(nm,"__ElasticNet")]] <- pe
  PREDC[[paste0(nm,"__XGBoost")]]    <- px
  res <- rbind(res,
    data.table(especificacion=nm, algoritmo="ElasticNet",
               R2=round(r2(yC,pe),4), RMSE=round(rmse(yC,pe),4)),
    data.table(especificacion=nm, algoritmo="XGBoost",
               R2=round(r2(yC,px),4), RMSE=round(rmse(yC,px),4)))
  cat(sprintf("  %-22s EN R2 %.4f RMSE %.4f | XGB R2 %.4f RMSE %.4f\n",
              nm, r2(yC,pe), rmse(yC,pe), r2(yC,px), rmse(yC,px)))
}
bs <- res[especificacion=="S2_base", .(algoritmo, R2_base=R2)]
res <- merge(res, bs, by="algoritmo"); res[, dR2_vs_S2 := round(R2-R2_base,5)]
setorder(res, algoritmo, especificacion)
cat(sprintf("\n=== GUANGZHOU (n = %d) ===\n", nrow(mC)))
print(res[, .(algoritmo, especificacion, R2, RMSE, dR2_vs_S2)], row.names=FALSE)

cv_sub <- function(idx, etiqueta) {
  if (length(idx) < 100) return(NULL)
  Xs <- as.matrix(DD[idx, vclin, with=FALSE]); ys <- yD[idx]
  set.seed(SEED); fo <- sample(rep(1:K, length.out=length(ys)))
  oe <- numeric(length(ys)); ox <- numeric(length(ys))
  for (k in 1:K) {
    tr <- fo!=k; te <- fo==k
    oe[te] <- fit_enet(Xs[tr,,drop=FALSE], ys[tr], Xs[te,,drop=FALSE])
    ox[te] <- fit_xgb (Xs[tr,,drop=FALSE], ys[tr], Xs[te,,drop=FALSE])
  }
  rbind(data.table(subgrupo=etiqueta, n=length(ys), algoritmo="ElasticNet",
                   R2=round(r2(ys,oe),4), RMSE=round(rmse(ys,oe),4)),
        data.table(subgrupo=etiqueta, n=length(ys), algoritmo="XGBoost",
                   R2=round(r2(ys,ox),4), RMSE=round(rmse(ys,ox),4)))
}
lo <- quantile(mC$Age,0.05); hi <- quantile(mC$Age,0.95)
banco <- rbindlist(list(
  cv_sub(seq_len(nrow(mD)), "BGI completo"),
  cv_sub(which(mD$mujer==1), "BGI mujeres"),
  cv_sub(which(mD$mujer==1 & mD$Age>=lo & mD$Age<=hi),
         sprintf("BGI mujeres %d-%d anos", round(lo), round(hi)))))
cat(sprintf("\n=== banco de referencia: mismo modelo de %d predictores, CV dentro de BGI ===\n",
            length(vclin)))
print(banco, row.names=FALSE)
cat("Guangzhou debe leerse contra 'BGI mujeres', no contra el BGI completo.\n")
fwrite(banco, file.path(P$out,"guangzhou_reference_bank.csv"))

res[, c("R2_ref_mujeres","RMSE_ref_mujeres") :=
      .(banco[subgrupo=="BGI mujeres", R2][match(algoritmo, banco[subgrupo=="BGI mujeres", algoritmo])],
        banco[subgrupo=="BGI mujeres", RMSE][match(algoritmo, banco[subgrupo=="BGI mujeres", algoritmo])])]
fwrite(res, file.path(P$out,"guangzhou_results.csv"))
fwrite(PREDC, file.path(P$out,"guangzhou_predictions.csv"))

set.seed(SEED); nC <- nrow(PREDC)
idxb <- replicate(B, sample.int(nC,nC,replace=TRUE), simplify=FALSE)
cols <- setdiff(names(PREDC), c("IID","y"))
R2B <- sapply(cols, function(cn){ p <- PREDC[[cn]]
  vapply(idxb, function(i) r2(yC[i],p[i]), numeric(1)) })
ic <- rbindlist(lapply(c("ElasticNet","XGBoost"), function(a)
  rbindlist(lapply(names(ESPEC), function(nm){
    c1 <- paste0(nm,"__",a); c0 <- paste0("S2_base__",a)
    d <- R2B[,c1]-R2B[,c0]; qq <- quantile(d,c(0.025,0.975))
    data.table(algoritmo=a, especificacion=nm,
               dR2=round(r2(yC,PREDC[[c1]])-r2(yC,PREDC[[c0]]),4),
               IC95_inf=round(qq[1],4), IC95_sup=round(qq[2],4),
               cruza_cero=(qq[1]<=0 & qq[2]>=0)) }))))
cat("\n=== Guangzhou con intervalos (bootstrap B = 2000) ===\n")
print(ic, row.names=FALSE)
fwrite(ic, file.path(P$out,"guangzhou_ci.csv"))

r2recal <- function(o,p){ ft <- lm(o ~ p); 1 - sum(residuals(ft)^2)/sum((o-mean(o))^2) }

cal <- rbindlist(lapply(cols, function(cn){
  pp <- PREDC[[cn]]; ft <- lm(yC ~ pp)
  data.table(modelo=cn,
             R2        = round(r2(yC,pp),4),
             RMSE      = round(rmse(yC,pp),4),
             sesgo     = round(mean(pp)-mean(yC),4),
             pendiente = round(coef(ft)[2],3),
             sd_pred   = round(sd(pp),4),
             r         = round(cor(yC,pp),4),
             R2_recal  = round(r2recal(yC,pp),4)) }))
cat("\n=== calibracion y discriminacion en Guangzhou ===\n"); print(cal, row.names=FALSE)
cat(sprintf("sd de la BMD observada en Guangzhou: %.4f\n", sd(yC)))
fwrite(cal, file.path(P$out,"guangzhou_calibration.csv"))

recal <- rbindlist(lapply(c("ElasticNet","XGBoost"), function(a)
  rbindlist(lapply(names(ESPEC), function(nm){
    c1 <- paste0(nm,"__",a); c0 <- paste0("S2_base__",a)
    dd <- vapply(idxb, function(i)
            r2recal(yC[i],PREDC[[c1]][i]) - r2recal(yC[i],PREDC[[c0]][i]), numeric(1))
    qq <- quantile(dd,c(0.025,0.975))
    data.table(algoritmo=a, especificacion=nm,
               R2_recal  = round(r2recal(yC,PREDC[[c1]]),4),
               dR2_recal = round(r2recal(yC,PREDC[[c1]])-r2recal(yC,PREDC[[c0]]),4),
               IC95_inf=round(qq[1],4), IC95_sup=round(qq[2],4),
               cruza_cero=(qq[1]<=0 & qq[2]>=0)) }))))
cat("\n=== dR2 sobre la escala RECALIBRADA (el que debe interpretarse) ===\n")
print(recal, row.names=FALSE)
fwrite(recal, file.path(P$out,"guangzhou_recalibration.csv"))

cat("\nguardado en results/\n")
cat("===== FIN OK =====\n")
