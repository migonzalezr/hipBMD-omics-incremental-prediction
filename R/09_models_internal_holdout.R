# ======================================================================
# 09_models_internal_holdout.R — Prediction models and internal holdout
#
# Fits the seven nested specifications (Elastic Net, alpha = 0.25; XGBoost) in the discovery set and applies the frozen models to the internal holdout.
#
# Inputs : clinical data, PCs, frozen PRS and MRS
# Outputs: results/holdout_predictions.csv; holdout_results.csv; holdout_score_association.csv
# Paper  : Table 2; Figure 4
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)
library(MASS); library(compositions); library(pls)
library(glmnet); library(xgboost)

source("config/paths.R")
P <- list(
  mgs  = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  mast = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  bmd  = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  clin = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  gpca = file.path(BASE,"data/genotypes/los_qc_pca.eigenvec"),
  disc = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  hold = file.path(BASE,"derived/analysis_sets/ids_internal_holdout.txt"),
  prs  = file.path(BASE,"derived/prs_los/los_prs.PT_0.01.sscore"),
  prsq = file.path(BASE,"derived/prs_external/external_prs.PT_0.1.sscore"),
  rds  = file.path(BASE,"derived/mrs/mrs_frozen_crudo.rds"),
  ass  = file.path(BASE,"results/balance_association_crudo.csv"),
  fun  = file.path(BASE,"external/PLS-PBs"),
  out  = file.path(BASE,"results")
)
for (f in c("fBalChip_PLS.r","fBPMaxOrthNewChip_PLS.R",
            "fBalChipman_PLS.R","fBPUpChi_PLS.r")) source(file.path(P$fun, f))

K <- 5; SEED <- 2026; PREV <- 0.01
PT_MRS <- 0.05
CLIN3 <- c("Age","gender_dummy","Weight_kg")
CLIN  <- c(CLIN3,"Height_cm","Race_African","Race_Caucasian","Race_Asian","Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
           "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")

disc_id <- unique(trimws(readLines(P$disc)))
hold_id <- unique(trimws(readLines(P$hold)))
stopifnot(length(intersect(disc_id, hold_id)) == 0)
cat(sprintf("discovery %d | holdout %d | disjuntos: OK\n", length(disc_id), length(hold_id)))

mgs <- fread(file = P$mgs); setnames(mgs, 1, "Metagenome_ID")
mgs <- fix_metagenome_ids(mgs)
pte <- as.data.table(read_excel(P$mast, sheet = "bgi"))
setnames(pte, "BGI_1794", "Subject_ID")
pte <- pte[, .(Metagenome_ID=trimws(as.character(Metagenome_ID)), IID=trimws(as.character(Subject_ID)))]
dt <- merge(pte[IID %in% c(disc_id, hold_id)], mgs, by="Metagenome_ID")

b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)),
                        y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
cl <- fread(file=P$clin); cl[, IID := trimws(as.character(Subject_ID))]
cl <- unique(cl[, c("IID",CLIN), with=FALSE], by="IID")
gp <- fread(file=P$gpca, header=FALSE)
setnames(gp, c("FID","IID", paste0("gPC", seq_len(ncol(gp)-2))))
gp <- gp[, .(IID=trimws(as.character(IID)), gPC1, gPC2)]

leer_sscore <- function(f, nombre) {
  s <- fread(f); setnames(s, sub("^#","",names(s)))
  out <- data.table(IID=trimws(as.character(s$IID)), v=as.numeric(s$SCORE1_AVG))
  setnames(out, "v", nombre); out
}
prs  <- leer_sscore(P$prs,  "PRS")
prsq <- leer_sscore(P$prsq, "PRS_qian")

meta <- Reduce(function(x,z) merge(x,z,by="IID"),
               list(dt[,.(IID)], ph, cl, gp, prs, prsq))
meta <- meta[complete.cases(meta)]
dt <- dt[IID %in% meta$IID]
setkey(dt,"IID"); setkey(meta,"IID"); dt <- dt[meta$IID]
stopifnot(identical(dt$IID, meta$IID))
RA <- as.matrix(dt[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RA) <- dt$IID
y  <- meta$y
es_disc <- meta$IID %in% disc_id
cat(sprintf("con datos completos -> discovery %d | holdout %d\n\n",
            sum(es_disc), sum(!es_disc)))

cerrar <- function(M,ps){ M[M==0] <- ps; (M/rowSums(M))*100 }
ajustar_prep <- function(R){
  s <- colnames(R)[apply(R,2,var,na.rm=TRUE)>0]
  s <- s[colMeans(R[,s,drop=FALSE]!=0) >= PREV]
  list(species=s, pseudo=0.5*min(R[,s,drop=FALSE][R[,s,drop=FALSE]!=0]))
}
aplicar_prep <- function(R,pr){
  M <- matrix(0,nrow(R),length(pr$species),dimnames=list(rownames(R),pr$species))
  p <- intersect(pr$species,colnames(R)); M[,p] <- R[,p,drop=FALSE]; cerrar(M,pr$pseudo)
}
pbscores <- function(Xc,Bt) as.matrix(compositions::clr(Xc)) %*% Bt

fr  <- readRDS(P$rds)
ass <- fread(file = P$ass)
sel <- ass[!is.na(p) & p < PT_MRS]
cat(sprintf("MRS congelado: %d balances con p < %.2g (de %d)\n", nrow(sel), PT_MRS, nrow(ass)))
stopifnot(nrow(sel) > 0, all(sel$PB_ID %in% colnames(fr$Bt)))

mrs_frozen <- function(idx) {
  PB <- pbscores(aplicar_prep(RA[idx,,drop=FALSE], fr$prep), fr$Bt)
  as.numeric(PB[, sel$PB_ID, drop=FALSE] %*% sel$beta)
}

idx_d <- which(es_disc)
set.seed(SEED)
q <- cut(y[idx_d], breaks=quantile(y[idx_d],probs=seq(0,1,length.out=K+1)), include.lowest=TRUE)
fold <- integer(length(idx_d))
for(lv in levels(q)){ i <- which(q==lv); fold[i] <- sample(rep(1:K,length.out=length(i))) }

mrs_oof <- rep(NA_real_, length(idx_d))
for (k in 1:K) {
  tr <- idx_d[fold!=k]; te <- idx_d[fold==k]
  pr <- ajustar_prep(RA[tr,,drop=FALSE]); Xtr <- aplicar_prep(RA[tr,,drop=FALSE], pr)
  set.seed(SEED)
  bal <- fBalChip_PLS(Xtr, y[tr], version="cov")
  Bt  <- t(bal$bal)[pr$species,,drop=FALSE]
  PBtr <- pbscores(Xtr, Bt)
  cov <- data.table(y=y[tr], meta[tr, CLIN3, with=FALSE])
  pcm <- prcomp(as.matrix(compositions::clr(RA[tr,,drop=FALSE]+1e-6)), center=TRUE, scale.=FALSE)
  cov[, `:=`(mPC1=pcm$x[,1], mPC2=pcm$x[,2])]
  ss <- rbindlist(lapply(colnames(PBtr), function(nm){
    co <- summary(lm(y ~ pb + Age + gender_dummy + Weight_kg + mPC1 + mPC2,
                     data=cbind(cov, pb=PBtr[,nm])))$coefficients
    if(!"pb" %in% rownames(co)) return(data.table(PB_ID=nm,beta=NA_real_,p=NA_real_))
    data.table(PB_ID=nm, beta=co["pb","Estimate"], p=co["pb","Pr(>|t|)"])
  }))
  s <- ss[!is.na(p) & p < PT_MRS]
  PBte <- pbscores(aplicar_prep(RA[te,,drop=FALSE], pr), Bt)
  mrs_oof[fold==k] <- if (nrow(s)) as.numeric(PBte[, s$PB_ID, drop=FALSE] %*% s$beta) else 0
  cat(sprintf("  fold %d: %d balances\n", k, nrow(s)))
}

MRS <- numeric(length(y))
MRS[idx_d]  <- mrs_oof
MRS[!es_disc] <- mrs_frozen(which(!es_disc))

estand <- function(v) { m <- mean(v[es_disc]); s <- sd(v[es_disc]); (v - m)/s }
meta[, `:=`(PRS = estand(PRS), PRS_qian = estand(PRS_qian))]
MRS <- estand(MRS)

Xc <- model.matrix(~ ., data=as.data.frame(meta[,c(CLIN,"gPC1","gPC2"),with=FALSE]))[,-1,drop=FALSE]
colnames(Xc) <- make.names(colnames(Xc), unique=TRUE)
Xc <- Xc[, apply(Xc,2,function(z) length(unique(z))>1), drop=FALSE]
vclin <- setdiff(colnames(Xc), c("gPC1","gPC2"))

D <- cbind(as.data.table(Xc),
           data.table(PRS=meta$PRS, PRS_qian=meta$PRS_qian, MRS=MRS))

ESPEC <- list(
  S1_clinicas          = vclin,
  S2_base              = c(vclin,"gPC1","gPC2"),
  S3_base_PRS          = c(vclin,"gPC1","gPC2","PRS"),
  S3q_base_PRSqian     = c(vclin,"gPC1","gPC2","PRS_qian"),
  S4_base_MRS          = c(vclin,"gPC1","gPC2","MRS"),
  S5_base_PRS_MRS      = c(vclin,"gPC1","gPC2","PRS","MRS"),

  S5q_base_PRSqian_MRS = c(vclin,"gPC1","gPC2","PRS_qian","MRS")
)

r2 <- function(obs,pred) 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)

elegir_alpha <- function(Xtr, ytr) {
  alphas <- c(0,0.25,0.5,0.75,1); best <- NULL
  for (a in alphas) {
    set.seed(SEED)
    cv <- cv.glmnet(Xtr, ytr, alpha=a, nfolds=10, standardize=TRUE)
    if (is.null(best) || min(cv$cvm) < best$cvm) best <- list(cvm=min(cv$cvm), alpha=a)
  }
  best$alpha
}
fit_enet <- function(Xtr,ytr,Xte,alpha){
  set.seed(SEED)
  cv <- cv.glmnet(Xtr, ytr, alpha=alpha, nfolds=10, standardize=TRUE)
  list(pred=as.numeric(predict(cv, newx=Xte, s="lambda.min")), alpha=alpha)
}
fit_xgb <- function(Xtr,ytr,Xte){
  dtr <- xgb.DMatrix(Xtr, label=ytr)
  pars <- list(objective="reg:squarederror", eta=0.05, max_depth=3,
               subsample=0.8, colsample_bytree=0.8, min_child_weight=5, nthread=2)
  set.seed(SEED)
  cvres <- xgb.cv(pars, dtr, nrounds=400, nfold=5, early_stopping_rounds=30, verbose=0)
  nb <- tryCatch({ lg <- as.data.frame(cvres$evaluation_log)
                   which.min(lg[[grep("^test.*mean$", names(lg), value=TRUE)[1]]]) },
                 error=function(e) NA_integer_)
  if (length(nb)!=1 || is.na(nb) || nb<1) nb <- 100
  set.seed(SEED)
  m <- xgb.train(pars, dtr, nrounds=nb, verbose=0)
  list(pred=predict(m, xgb.DMatrix(Xte)), nrounds=nb)
}

ALPHA <- elegir_alpha(as.matrix(D[es_disc, ESPEC$S2_base, with=FALSE]), y[es_disc])
cat(sprintf("\nalpha elegido sobre S2 en el discovery y FIJADO para todas: %.2f\n", ALPHA))

cat("\n>>> ABRIENDO EL HOLDOUT (una sola vez) <<<\n\n")
res <- data.table()
PRED <- data.table(IID = meta$IID[!es_disc], y = y[!es_disc])
for (nm in names(ESPEC)) {
  v <- ESPEC[[nm]]
  Xtr <- as.matrix(D[es_disc, v, with=FALSE]); Xte <- as.matrix(D[!es_disc, v, with=FALSE])
  ytr <- y[es_disc]; yte <- y[!es_disc]
  e <- fit_enet(Xtr,ytr,Xte,ALPHA); x <- fit_xgb(Xtr,ytr,Xte)
  res <- rbind(res,
    data.table(especificacion=nm, algoritmo="ElasticNet",
               R2_holdout=round(r2(yte,e$pred),4), detalle=sprintf("alpha=%.2f", e$alpha)),
    data.table(especificacion=nm, algoritmo="XGBoost",
               R2_holdout=round(r2(yte,x$pred),4), detalle=sprintf("nrounds=%d", x$nrounds)))
  PRED[[paste0(nm,"__ElasticNet")]] <- e$pred
  PRED[[paste0(nm,"__XGBoost")]]    <- x$pred
  cat(sprintf("  %-20s EN %.4f | XGB %.4f\n", nm, r2(yte,e$pred), r2(yte,x$pred)))
}
fwrite(PRED, file.path(P$out,"holdout_predictions.csv"))

base <- res[especificacion=="S2_base", .(algoritmo, R2_base=R2_holdout)]
res <- merge(res, base, by="algoritmo")
res[, dR2_vs_S2 := round(R2_holdout - R2_base, 5)]
setorder(res, algoritmo, especificacion)

cat("\n=== RESULTADOS EN EL HOLDOUT (n = ", sum(!es_disc), ") ===\n", sep="")
print(res[, .(algoritmo, especificacion, R2_holdout, dR2_vs_S2, detalle)], row.names=FALSE)
fwrite(res, file.path(P$out,"holdout_results.csv"))

vb <- c(vclin,"gPC1","gPC2")
assoc <- rbindlist(lapply(c("PRS","PRS_qian","MRS"), function(s) {
  rbindlist(lapply(c("discovery","holdout"), function(cj) {
    i <- if (cj=="discovery") es_disc else !es_disc
    d <- data.table(y=y[i], D[i, c(vb,s), with=FALSE])
    fit <- lm(as.formula(paste("y ~", paste(c(vb,s), collapse="+"))), data=d)
    co <- summary(fit)$coefficients; ci <- confint(fit, s, level=0.95)
    data.table(score=s, conjunto=cj, n=sum(i),
               beta=round(co[s,"Estimate"],5), se=round(co[s,"Std. Error"],5),
               IC95_inf=round(ci[1],5), IC95_sup=round(ci[2],5),
               P=signif(co[s,"Pr(>|t|)"],3))
  }))
}))
cat("\n=== R1.3: asociacion de cada score con HTOT_BMD (por DE del score) ===\n")
print(assoc, row.names=FALSE)
fwrite(assoc, file.path(P$out,"holdout_score_association.csv"))

cat("\nNOTA OBLIGATORIA PARA EL MANUSCRITO:\n",
    " el MRS pre-especificado (FDR<0.05) selecciono CERO balances en el\n",
    " discovery, de modo que no existe y no entra en esta tabla. El MRS\n",
    " que aparece aqui es el secundario por umbral (Enmienda 4).\n")
cat("\nguardado en results/holdout_results.csv\n")
cat("===== FIN OK =====\n")
