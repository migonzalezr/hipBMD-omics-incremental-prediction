# ======================================================================
# 19_model_interpretation.R — Model interpretation
#
# Permutation importance on out-of-sample R2, SHAP values (XGBoost) and standardized coefficients (Elastic Net) for the frozen models in the holdout.
#
# Inputs : frozen models; holdout data
# Outputs: results/interpretation_*.csv
# Paper  : Figure S6; Table S14
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
K <- 5; SEED <- 2026; PREV <- 0.01; PT_MRS <- 0.05; ALPHA <- 0.25; NPERM <- 50
CLIN3 <- c("Age","gender_dummy","Weight_kg")
CLIN  <- c(CLIN3,"Height_cm","BMI","Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
           "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")

disc_id <- unique(trimws(readLines(P$disc))); hold_id <- unique(trimws(readLines(P$hold)))
stopifnot(length(intersect(disc_id, hold_id)) == 0)
mgs <- fread(file=P$mgs); setnames(mgs,1,"Metagenome_ID")
mgs <- fix_metagenome_ids(mgs)
pte <- as.data.table(read_excel(P$mast, sheet="bgi"))
setnames(pte,"BGI_1794","Subject_ID")
pte <- pte[, .(Metagenome_ID=trimws(as.character(Metagenome_ID)), IID=trimws(as.character(Subject_ID)))]
dt <- merge(pte[IID %in% c(disc_id,hold_id)], mgs, by="Metagenome_ID")
b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)),
                        y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
cl <- fread(file=P$clin); cl[, IID := trimws(as.character(Subject_ID))]
cl <- unique(cl[, c("IID",CLIN), with=FALSE], by="IID")
gp <- fread(file=P$gpca, header=FALSE)
setnames(gp, c("FID","IID", paste0("gPC", seq_len(ncol(gp)-2))))
gp <- gp[, .(IID=trimws(as.character(IID)), gPC1, gPC2)]
leer_sscore <- function(f,nm){ s <- fread(f); setnames(s, sub("^#","",names(s)))
  o <- data.table(IID=trimws(as.character(s$IID)), v=as.numeric(s$SCORE1_AVG)); setnames(o,"v",nm); o }
meta <- Reduce(function(x,z) merge(x,z,by="IID"),
               list(dt[,.(IID)], ph, cl, gp, leer_sscore(P$prs,"PRS"), leer_sscore(P$prsq,"PRS_qian")))
meta <- meta[complete.cases(meta)]
dt <- dt[IID %in% meta$IID]; setkey(dt,"IID"); setkey(meta,"IID"); dt <- dt[meta$IID]
stopifnot(identical(dt$IID, meta$IID))
RA <- as.matrix(dt[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RA) <- dt$IID
y <- meta$y; es_disc <- meta$IID %in% disc_id
cat(sprintf("discovery %d | holdout %d\n", sum(es_disc), sum(!es_disc)))

cerrar <- function(M,ps){ M[M==0] <- ps; (M/rowSums(M))*100 }
aplicar_prep <- function(R,pr){
  M <- matrix(0,nrow(R),length(pr$species),dimnames=list(rownames(R),pr$species))
  p <- intersect(pr$species,colnames(R)); M[,p] <- R[,p,drop=FALSE]; cerrar(M,pr$pseudo) }
pbscores <- function(Xc,Bt) as.matrix(compositions::clr(Xc)) %*% Bt
ajustar_prep <- function(R){
  s <- colnames(R)[apply(R,2,var,na.rm=TRUE)>0]
  s <- s[colMeans(R[,s,drop=FALSE]!=0) >= PREV]
  list(species=s, pseudo=0.5*min(R[,s,drop=FALSE][R[,s,drop=FALSE]!=0])) }
fr <- readRDS(P$rds); ass <- fread(file=P$ass); sel <- ass[!is.na(p) & p < PT_MRS]
stopifnot(nrow(sel) > 0, all(sel$PB_ID %in% colnames(fr$Bt)))

mrs_frozen <- function(idx){
  PB <- pbscores(aplicar_prep(RA[idx,,drop=FALSE], fr$prep), fr$Bt)
  as.numeric(PB[, sel$PB_ID, drop=FALSE] %*% sel$beta) }

idx_d <- which(es_disc)
set.seed(SEED)
q <- cut(y[idx_d], breaks=quantile(y[idx_d], probs=seq(0,1,length.out=K+1)), include.lowest=TRUE)
fold <- integer(length(idx_d))
for (lv in levels(q)) { i <- which(q==lv); fold[i] <- sample(rep(1:K, length.out=length(i))) }
mrs_oof <- rep(NA_real_, length(idx_d))
for (k in 1:K) {
  tr <- idx_d[fold!=k]; te <- idx_d[fold==k]
  pr <- ajustar_prep(RA[tr,,drop=FALSE]); Xtr <- aplicar_prep(RA[tr,,drop=FALSE], pr)
  set.seed(SEED); bal <- fBalChip_PLS(Xtr, y[tr], version="cov")
  Bt <- t(bal$bal)[pr$species,,drop=FALSE]; PBtr <- pbscores(Xtr, Bt)
  pcm <- prcomp(as.matrix(compositions::clr(RA[tr,,drop=FALSE]+1e-6)), center=TRUE, scale.=FALSE)
  cov <- data.table(y=y[tr], meta[tr, CLIN3, with=FALSE], mPC1=pcm$x[,1], mPC2=pcm$x[,2])
  ss <- rbindlist(lapply(colnames(PBtr), function(nm){
    co <- summary(lm(y ~ pb + Age + gender_dummy + Weight_kg + mPC1 + mPC2,
                     data=cbind(cov, pb=PBtr[,nm])))$coefficients
    if(!"pb" %in% rownames(co)) return(data.table(PB_ID=nm, beta=NA_real_, p=NA_real_))
    data.table(PB_ID=nm, beta=co["pb","Estimate"], p=co["pb","Pr(>|t|)"]) }))
  s <- ss[!is.na(p) & p < PT_MRS]
  PBte <- pbscores(aplicar_prep(RA[te,,drop=FALSE], pr), Bt)
  mrs_oof[fold==k] <- if (nrow(s)) as.numeric(PBte[, s$PB_ID, drop=FALSE] %*% s$beta) else 0
  cat(sprintf("  fold %d: %d balances\n", k, nrow(s)))
}
MRS <- numeric(length(y))
MRS[idx_d]    <- mrs_oof
MRS[!es_disc] <- mrs_frozen(which(!es_disc))

estand <- function(v){ m <- mean(v[es_disc]); s <- sd(v[es_disc]); (v-m)/s }
meta[, `:=`(PRS=estand(PRS), PRS_qian=estand(PRS_qian))]; MRS <- estand(MRS)

Xc <- model.matrix(~ ., data=as.data.frame(meta[,c(CLIN,"gPC1","gPC2"),with=FALSE]))[,-1,drop=FALSE]
colnames(Xc) <- make.names(colnames(Xc), unique=TRUE)
Xc <- Xc[, apply(Xc,2,function(z) length(unique(z))>1), drop=FALSE]
vclin <- setdiff(colnames(Xc), c("gPC1","gPC2"))
D <- cbind(as.data.table(Xc), data.table(PRS=meta$PRS, PRS_qian=meta$PRS_qian, MRS=MRS))

ESPEC <- list(S5q_base_PRSqian_MRS = c(vclin,"gPC1","gPC2","PRS_qian","MRS"),
              S5_base_PRS_MRS      = c(vclin,"gPC1","gPC2","PRS","MRS"))
r2 <- function(o,p) 1 - sum((o-p)^2)/sum((o-mean(o))^2)

COEF <- data.table(); PERM <- data.table(); SHAP <- data.table()
for (nm in names(ESPEC)) {
  v   <- ESPEC[[nm]]
  Xtr <- as.matrix(D[es_disc,  v, with=FALSE]); ytr <- y[es_disc]
  Xte <- as.matrix(D[!es_disc, v, with=FALSE]); yte <- y[!es_disc]

  set.seed(SEED)
  cvg <- cv.glmnet(Xtr, ytr, alpha=ALPHA, nfolds=10, standardize=TRUE)
  cf  <- as.matrix(coef(cvg, s="lambda.min"))
  sdx <- apply(Xtr, 2, sd)
  co  <- data.table(especificacion=nm, variable=rownames(cf)[-1],
                    coef=as.numeric(cf[-1,1]))
  co[, coef_estandarizado := coef * sdx[variable]]
  COEF <- rbind(COEF, co)
  pe <- as.numeric(predict(cvg, newx=Xte, s="lambda.min"))

  dtr  <- xgb.DMatrix(Xtr, label=ytr)
  pars <- list(objective="reg:squarederror", eta=0.05, max_depth=3,
               subsample=0.8, colsample_bytree=0.8, min_child_weight=5)
  set.seed(SEED)
  cvx <- xgb.cv(pars, dtr, nrounds=400, nfold=5, early_stopping_rounds=30, verbose=0)
  nb <- tryCatch({ lg <- as.data.frame(cvx$evaluation_log)
                   which.min(lg[[grep("^test.*mean$", names(lg), value=TRUE)[1]]]) },
                 error=function(e) NA_integer_)
  if (length(nb)!=1 || is.na(nb) || nb<1) nb <- 100
  set.seed(SEED); mx <- xgb.train(pars, dtr, nrounds=nb, verbose=0)
  px <- predict(mx, xgb.DMatrix(Xte))
  cat(sprintf("%s: R2 holdout  EN %.4f | XGB %.4f\n", nm, r2(yte,pe), r2(yte,px)))

  if (nm == "S5q_base_PRSqian_MRS" && abs(r2(yte,pe) - 0.3928) > 0.01)
    warning(sprintf("S5q EN da %.4f y 09_models_internal_holdout.R daba 0.3928: revisar", r2(yte,pe)))

  base_en <- r2(yte, pe); base_xg <- r2(yte, px)
  set.seed(SEED)
  pm <- rbindlist(lapply(v, function(cn){
    de <- dx <- numeric(NPERM)
    for (b_ in seq_len(NPERM)) {
      Xp <- Xte; Xp[,cn] <- sample(Xp[,cn])
      de[b_] <- base_en - r2(yte, as.numeric(predict(cvg, newx=Xp, s="lambda.min")))
      dx[b_] <- base_xg - r2(yte, predict(mx, xgb.DMatrix(Xp)))
    }
    data.table(especificacion=nm, variable=cn,
               perdida_R2_EN=mean(de),  sd_EN=sd(de),
               perdida_R2_XGB=mean(dx), sd_XGB=sd(dx))
  }))
  PERM <- rbind(PERM, pm)

  sv <- predict(mx, xgb.DMatrix(Xte), predcontrib=TRUE)
  colnames(sv) <- c(colnames(Xte), "BIAS")
  SHAP <- rbind(SHAP, data.table(especificacion=nm, variable=colnames(Xte),
                                 shap_medio_abs=colMeans(abs(sv[, colnames(Xte), drop=FALSE]))))
}

setorder(PERM, especificacion, -perdida_R2_EN)
cat("\n=== importancia por permutacion (perdida de R2 en el holdout) ===\n")
print(PERM[especificacion=="S5q_base_PRSqian_MRS"][1:12], row.names=FALSE, digits=3)
cat("\n=== coeficientes estandarizados del Elastic Net, |mayores| ===\n")
print(COEF[especificacion=="S5q_base_PRSqian_MRS"][order(-abs(coef_estandarizado))][1:12],
      row.names=FALSE, digits=3)
cat("\n=== SHAP medio absoluto (XGBoost) ===\n")
print(SHAP[especificacion=="S5q_base_PRSqian_MRS"][order(-shap_medio_abs)][1:12],
      row.names=FALSE, digits=3)
cat("\nPRS y MRS en las tres vistas:\n")
print(merge(merge(COEF[variable %in% c("PRS_qian","MRS","PRS"), .(especificacion,variable,coef_estandarizado)],
                  PERM[, .(especificacion,variable,perdida_R2_EN,perdida_R2_XGB)], by=c("especificacion","variable")),
            SHAP, by=c("especificacion","variable")), row.names=FALSE, digits=3)

fwrite(COEF, file.path(P$out,"interpretation_coefficients.csv"))
fwrite(PERM, file.path(P$out,"interpretation_permutation.csv"))
fwrite(SHAP, file.path(P$out,"interpretation_shap.csv"))
cat("\nguardado en results/\n")
cat("===== FIN OK =====\n")
