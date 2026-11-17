# ======================================================================
# 11_discovery_cv_increments.R — Discovery-set increments
#
# Out-of-fold predictions and conditional contrasts in the modeling discovery set.
#
# Inputs : discovery data, frozen scores
# Outputs: results/discovery_predictions_oof.csv; discovery_conditional_contrast.csv
# Paper  : Table 2
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
  prs  = file.path(BASE,"derived/prs_los/los_prs.PT_0.01.sscore"),
  prsq = file.path(BASE,"derived/prs_external/external_prs.PT_0.1.sscore"),
  fun  = file.path(BASE,"external/PLS-PBs"),
  out  = file.path(BASE,"results")
)
for (f in c("fBalChip_PLS.r","fBPMaxOrthNewChip_PLS.R",
            "fBalChipman_PLS.R","fBPUpChi_PLS.r")) source(file.path(P$fun, f))
K <- 5; SEED <- 2026; PREV <- 0.01; PT_MRS <- 0.05; B <- 2000
CLIN3 <- c("Age","gender_dummy","Weight_kg")
CLIN  <- c(CLIN3,"Height_cm","BMI","Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
           "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")

disc <- unique(trimws(readLines(P$disc)))
mgs <- fread(file = P$mgs); setnames(mgs, 1, "Metagenome_ID")
mgs <- fix_metagenome_ids(mgs)
pte <- as.data.table(read_excel(P$mast, sheet = "bgi"))
setnames(pte, "BGI_1794", "Subject_ID")
pte <- pte[, .(Metagenome_ID=trimws(as.character(Metagenome_ID)), IID=trimws(as.character(Subject_ID)))]
dt <- merge(pte[IID %in% disc], mgs, by="Metagenome_ID")
b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)),
                        y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
cl <- fread(file=P$clin); cl[, IID := trimws(as.character(Subject_ID))]
cl <- unique(cl[, c("IID",CLIN), with=FALSE], by="IID")
gp <- fread(file=P$gpca, header=FALSE)
setnames(gp, c("FID","IID", paste0("gPC", seq_len(ncol(gp)-2))))
gp <- gp[, .(IID=trimws(as.character(IID)), gPC1, gPC2)]
leer_s <- function(f,nm){ s <- fread(f); setnames(s, sub("^#","",names(s)))
  o <- data.table(IID=trimws(as.character(s$IID)), v=as.numeric(s$SCORE1_AVG)); setnames(o,"v",nm); o }
meta <- Reduce(function(x,z) merge(x,z,by="IID"),
               list(dt[,.(IID)], ph, cl, gp, leer_s(P$prs,"PRS"), leer_s(P$prsq,"PRS_qian")))
meta <- meta[complete.cases(meta)]
dt <- dt[IID %in% meta$IID]
setkey(dt,"IID"); setkey(meta,"IID"); dt <- dt[meta$IID]
stopifnot(identical(dt$IID, meta$IID))
RA <- as.matrix(dt[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RA) <- dt$IID
y <- meta$y; n <- length(y)
cat(sprintf("discovery: %d individuos (el holdout NO se toca)\n", n))

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
pbs <- function(Xc,Bt) as.matrix(compositions::clr(Xc)) %*% Bt

set.seed(SEED)
q <- cut(y, breaks=quantile(y,probs=seq(0,1,length.out=K+1)), include.lowest=TRUE)
fold <- integer(n); for(lv in levels(q)){ i <- which(q==lv); fold[i] <- sample(rep(1:K,length.out=length(i))) }

MRS <- numeric(n)
for (k in 1:K) {
  tr <- which(fold!=k); te <- which(fold==k)
  pr <- ajustar_prep(RA[tr,,drop=FALSE]); Xtr <- aplicar_prep(RA[tr,,drop=FALSE], pr)
  set.seed(SEED); bal <- fBalChip_PLS(Xtr, y[tr], version="cov")
  Bt <- t(bal$bal)[pr$species,,drop=FALSE]; PBtr <- pbs(Xtr,Bt)
  pcm <- prcomp(as.matrix(compositions::clr(RA[tr,,drop=FALSE]+1e-6)), center=TRUE, scale.=FALSE)
  cov <- data.table(y=y[tr], meta[tr,CLIN3,with=FALSE], mPC1=pcm$x[,1], mPC2=pcm$x[,2])
  ss <- rbindlist(lapply(colnames(PBtr), function(nm){
    co <- summary(lm(y ~ pb + Age + gender_dummy + Weight_kg + mPC1 + mPC2,
                     data=cbind(cov, pb=PBtr[,nm])))$coefficients
    if(!"pb" %in% rownames(co)) return(data.table(PB_ID=nm,beta=NA_real_,p=NA_real_))
    data.table(PB_ID=nm, beta=co["pb","Estimate"], p=co["pb","Pr(>|t|)"]) }))
  s <- ss[!is.na(p) & p < PT_MRS]
  PBte <- pbs(aplicar_prep(RA[te,,drop=FALSE], pr), Bt)
  MRS[te] <- if (nrow(s)) as.numeric(PBte[, s$PB_ID, drop=FALSE] %*% s$beta) else 0
  cat(sprintf("  MRS fold %d: %d balances\n", k, nrow(s)))
}

z <- function(v) (v-mean(v))/sd(v)
Xc <- model.matrix(~ ., data=as.data.frame(meta[,c(CLIN,"gPC1","gPC2"),with=FALSE]))[,-1,drop=FALSE]
colnames(Xc) <- make.names(colnames(Xc), unique=TRUE)
Xc <- Xc[, apply(Xc,2,function(w) length(unique(w))>1), drop=FALSE]
vclin <- setdiff(colnames(Xc), c("gPC1","gPC2"))
D <- cbind(as.data.table(Xc), data.table(PRS=z(meta$PRS), PRS_qian=z(meta$PRS_qian), MRS=z(MRS)))

ESPEC <- list(
  S2_base              = c(vclin,"gPC1","gPC2"),
  S3_base_PRS          = c(vclin,"gPC1","gPC2","PRS"),
  S3q_base_PRSqian     = c(vclin,"gPC1","gPC2","PRS_qian"),
  S4_base_MRS          = c(vclin,"gPC1","gPC2","MRS"),
  S5_base_PRS_MRS      = c(vclin,"gPC1","gPC2","PRS","MRS"),
  S5q_base_PRSqian_MRS = c(vclin,"gPC1","gPC2","PRS_qian","MRS")
)

fit_enet <- function(Xtr,ytr,Xte,a){ set.seed(SEED)
  cv <- cv.glmnet(Xtr,ytr,alpha=a,nfolds=10,standardize=TRUE)
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

set.seed(SEED)
cvb <- cv.glmnet(as.matrix(D[, ESPEC$S2_base, with=FALSE]), y, alpha=0.25, nfolds=10)
ALPHA <- 0.25
cat(sprintf("\nalpha fijo = %.2f (igual que en el holdout)\n", ALPHA))

PRED <- data.table(IID=meta$IID, y=y, fold=fold)
for (nm in names(ESPEC)) {
  v <- ESPEC[[nm]]; pe <- px <- numeric(n)
  for (k in 1:K) {
    tr <- fold!=k; te <- fold==k
    Xtr <- as.matrix(D[tr, v, with=FALSE]); Xte <- as.matrix(D[te, v, with=FALSE])
    pe[te] <- fit_enet(Xtr, y[tr], Xte, ALPHA)
    px[te] <- fit_xgb (Xtr, y[tr], Xte)
  }
  PRED[[paste0(nm,"__ElasticNet")]] <- pe
  PRED[[paste0(nm,"__XGBoost")]]    <- px
  cat(sprintf("  %-22s EN %.4f | XGB %.4f\n", nm,
              1-sum((y-pe)^2)/sum((y-mean(y))^2), 1-sum((y-px)^2)/sum((y-mean(y))^2)))
}
fwrite(PRED, file.path(P$out,"discovery_predictions_oof.csv"))

r2 <- function(o,p) 1 - sum((o-p)^2)/sum((o-mean(o))^2)
set.seed(SEED); idx <- replicate(B, sample.int(n,n,replace=TRUE), simplify=FALSE)
R2B <- sapply(setdiff(names(PRED),c("IID","y","fold")),
              function(cn){ p <- PRED[[cn]]; vapply(idx, function(i) r2(y[i],p[i]), numeric(1)) })

CONTR <- list(
  "dR2(MRS | PRS LOS)"  = c("S5_base_PRS_MRS",      "S3_base_PRS"),
  "dR2(MRS | PRS Qian)" = c("S5q_base_PRSqian_MRS", "S3q_base_PRSqian"),
  "dR2(PRS Qian | MRS)" = c("S5q_base_PRSqian_MRS", "S4_base_MRS"),
  "dR2(MRS)"            = c("S4_base_MRS",          "S2_base"),
  "dR2(PRS Qian)"       = c("S3q_base_PRSqian",     "S2_base")
)
cond <- rbindlist(lapply(c("ElasticNet","XGBoost"), function(a)
  rbindlist(lapply(names(CONTR), function(nm){
    c1 <- paste0(CONTR[[nm]][1],"__",a); c0 <- paste0(CONTR[[nm]][2],"__",a)
    d <- R2B[,c1]-R2B[,c0]; ic <- quantile(d, c(0.025,0.975))
    data.table(algoritmo=a, contraste=nm,
               dR2=round(r2(y,PRED[[c1]])-r2(y,PRED[[c0]]),4),
               IC95_inf=round(ic[1],4), IC95_sup=round(ic[2],4),
               cruza_cero=(ic[1]<=0 & ic[2]>=0)) }))))

cat(sprintf("\n=== CONTRASTES EN EL DISCOVERY (n = %d, fuera de fold) ===\n", n))
print(cond, row.names=FALSE)
fwrite(cond, file.path(P$out,"discovery_conditional_contrast.csv"))

cat("\n--- comparacion con el holdout (n = 351) ---\n")
cat("  holdout, XGBoost, dR2(MRS | PRS Qian) = +0.0109  IC [0.0011, 0.0208]\n")
cat("  holdout, ElasticNet, mismo contraste  = +0.0011  IC [-0.0036, 0.0053]\n")
cat("\nSi aqui el IC cruza cero, aquel resultado era ruido.\n")
cat("===== FIN OK =====\n")
