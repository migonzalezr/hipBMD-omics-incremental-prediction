# ======================================================================
# 08_mrs_threshold.R — MRS threshold curve
#
# Evaluates the cross-fitted MRS over P-value thresholds, analogous to the PRS threshold selection.
#
# Inputs : as 07
# Outputs: results/mrs_threshold_curve.csv
# Paper  : Experimental Section (MRS)
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)
library(MASS); library(compositions); library(pls)

source("config/paths.R")
P <- list(
  mgs  = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  mast = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  bmd  = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  clin = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  gpca = file.path(BASE,"data/genotypes/los_qc_pca.eigenvec"),
  disc = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  fun  = file.path(BASE,"external/PLS-PBs"),
  out  = file.path(BASE,"results")
)
dir.create(P$out, showWarnings = FALSE, recursive = TRUE)
for (f in c("fBalChip_PLS.r","fBPMaxOrthNewChip_PLS.R",
            "fBalChipman_PLS.R","fBPUpChi_PLS.r")) source(file.path(P$fun, f))

K <- 5; SEED <- 2026; PREV <- 0.01
UMBRALES <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)
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
meta <- Reduce(function(x,z) merge(x,z,by="IID"), list(dt[,.(IID)], ph, cl, gp))
meta <- meta[complete.cases(meta)]
dt <- dt[IID %in% meta$IID]
setkey(dt,"IID"); setkey(meta,"IID"); dt <- dt[meta$IID]
stopifnot(identical(dt$IID, meta$IID))
RA <- as.matrix(dt[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RA) <- dt$IID
y <- meta$y
cat(sprintf("discovery: %d individuos, %d especies crudas\n", nrow(RA), ncol(RA)))

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

bwas_training <- function(idx, residual){
  Rtr <- RA[idx,,drop=FALSE]
  pr  <- ajustar_prep(Rtr); Xtr <- aplicar_prep(Rtr,pr)
  pc  <- prcomp(as.matrix(compositions::clr(Rtr+1e-6)), center=TRUE, scale.=FALSE)
  cov <- data.table(y=y[idx], meta[idx,CLIN3,with=FALSE], mPC1=pc$x[,1], mPC2=pc$x[,2])
  obj <- if (residual)
    residuals(lm(y ~ Age + gender_dummy + Weight_kg + mPC1 + mPC2, data=cov)) else y[idx]
  bal <- fBalChip_PLS(Xtr, obj, version="cov")
  Bt  <- t(bal$bal)[pr$species,,drop=FALSE]
  PB  <- pbscores(Xtr,Bt)
  ss  <- rbindlist(lapply(colnames(PB), function(nm){
    co <- summary(lm(y ~ pb + Age + gender_dummy + Weight_kg + mPC1 + mPC2,
                     data=cbind(cov, pb=PB[,nm])))$coefficients
    if(!"pb" %in% rownames(co)) return(data.table(PB_ID=nm,beta=NA_real_,p=NA_real_))
    data.table(PB_ID=nm, beta=co["pb","Estimate"], p=co["pb","Pr(>|t|)"])
  }))
  list(prep=pr, Bt=Bt, ss=ss[!is.na(p)])
}

cat("\n=== la base de balances es ortogonal? ===\n")
set.seed(SEED)
f0 <- bwas_training(seq_len(nrow(RA)), FALSE)
G  <- crossprod(f0$Bt)
fo <- max(abs(G - diag(nrow(G))))
PBf <- pbscores(aplicar_prep(RA, f0$prep), f0$Bt)
Cc <- cor(PBf); diag(Cc) <- 0
cat(sprintf("  desviacion maxima de B'B respecto de la identidad: %.3g\n", fo))
cat(sprintf("  correlacion absoluta maxima entre puntajes de balances: %.3f\n", max(abs(Cc))))
cat(sprintf("  correlacion absoluta mediana: %.3f\n", median(abs(Cc[upper.tri(Cc)]))))
cat(if (fo < 1e-6) "  -> base ORTONORMAL: sumar betas por umbral es limpio, no hace falta clumping\n"
    else "  -> NO ortonormal: los balances comparten especies, hay doble conteo. Ver nota final\n")

set.seed(SEED)
q <- cut(y, breaks=quantile(y,probs=seq(0,1,length.out=K+1)), include.lowest=TRUE)
fold <- integer(length(y))
for(lv in levels(q)){ i <- which(q==lv); fold[i] <- sample(rep(1:K,length.out=length(i))) }

Xd <- model.matrix(~ ., data=as.data.frame(meta[,c(CLIN,"gPC1","gPC2"),with=FALSE]))[,-1,drop=FALSE]
colnames(Xd) <- make.names(colnames(Xd), unique=TRUE)
Xd <- Xd[, apply(Xd,2,function(z) length(unique(z))>1), drop=FALSE]
r2_oof <- function(d,f){
  pr <- rep(NA_real_,nrow(d))
  for(k in 1:K){ tr <- fold!=k; pr[!tr] <- predict(lm(f,data=d[tr]), newdata=d[!tr]) }
  1 - sum((d$y-pr)^2)/sum((d$y-mean(d$y))^2)
}
dbase <- data.table(y=y, as.data.table(Xd)); v <- setdiff(names(dbase),"y")
R0 <- r2_oof(dbase, as.formula(paste("y ~", paste(v,collapse="+"))))
cat(sprintf("\nR2 base (clinicas + gPC1 + gPC2), out-of-fold: %.4f\n", R0))

curva <- data.table()
for (modo in c("crudo","residual")) {
  cat(sprintf("\n=== variante: %s ===\n", modo))
  M <- matrix(0, nrow(RA), length(UMBRALES),
              dimnames=list(NULL, paste0("PT_",UMBRALES)))
  nbal <- matrix(0, K, length(UMBRALES))
  for (k in 1:K) {
    tr <- which(fold!=k); te <- which(fold==k)
    set.seed(SEED)
    fk <- bwas_training(tr, modo=="residual")
    PBte <- pbscores(aplicar_prep(RA[te,,drop=FALSE], fk$prep), fk$Bt)
    for (j in seq_along(UMBRALES)) {
      s <- fk$ss[p < UMBRALES[j]]
      nbal[k,j] <- nrow(s)
      if (nrow(s)) M[te,j] <- as.numeric(PBte[, s$PB_ID, drop=FALSE] %*% s$beta)
    }
    cat(sprintf("  fold %d: balances por umbral %s\n", k, paste(nbal[k,], collapse="/")))
  }
  for (j in seq_along(UMBRALES)) {
    d <- copy(dbase)
    if (sd(M[,j]) > 0) {
      d[, MRS := as.numeric(scale(M[,j]))]
      R1 <- r2_oof(d, as.formula(paste("y ~", paste(c(v,"MRS"),collapse="+"))))
    } else R1 <- R0
    curva <- rbind(curva, data.table(variante=modo, PT=UMBRALES[j],
                                     balances_medio=round(mean(nbal[,j]),1),
                                     R2_base=round(R0,4), R2_mas_MRS=round(R1,4),
                                     dR2=round(R1-R0,5)))
  }
  print(curva[variante==modo, .(PT, balances_medio, dR2)], row.names=FALSE)
}

cat("\n=== CURVA COMPLETA ===\n"); print(curva, row.names=FALSE)
fwrite(curva, file.path(P$out,"mrs_threshold_curve.csv"))

mejor <- curva[which.max(dR2)]
cat(sprintf("\nmejor: variante %s, PT %.4g, %.1f balances, dR2 %.5f\n",
            mejor$variante, mejor$PT, mejor$balances_medio, mejor$dR2))
cat(sprintf("referencia: el PRS propio de LOS dio dR2 = 0.00431 con PT_0.01\n"))
cat("\nNOTA: si la base NO resulto ortonormal, los balances comparten especies y\n",
    "la suma ponderada cuenta la misma senal varias veces. En ese caso el\n",
    "analogo del clumping seria quedarse con balances no anidados, o cambiar\n",
    "la suma por una regresion penalizada sobre los puntajes.\n")
cat("\nguardado en results/mrs_threshold_curve.csv\n")
