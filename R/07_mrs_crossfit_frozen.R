# ======================================================================
# 07_mrs_crossfit_frozen.R — Metagenomic risk score: cross-fitting and frozen model
#
# Builds principal balances (PLS) and the balance association scan within training folds (cross-fitting), then fits and freezes the full-discovery model for projection to other sets.
#
# Inputs : data/microbiome/los_bgi_species_abundance.csv; clinical data; external/PLS-PBs functions
# Outputs: derived/mrs/mrs_frozen_*.rds; results/balance_association_*.csv; mrs_crossfit_variants.csv
# Paper  : Experimental Section (MRS); Table S4
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
  mrs  = file.path(BASE,"derived/mrs"),
  out  = file.path(BASE,"results")
)
for (d in c(P$mrs,P$out)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
for (f in c("fBalChip_PLS.r","fBPMaxOrthNewChip_PLS.R",
            "fBalChipman_PLS.R","fBPUpChi_PLS.r")) source(file.path(P$fun, f))

K <- 5; SEED <- 2026; PREV <- 0.01; FDRCUT <- 0.05
CLIN3 <- c("Age","gender_dummy","Weight_kg")
CLIN  <- c(CLIN3,"Height_cm","BMI","Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
           "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")

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
cl <- unique(cl[, c("IID", CLIN), with = FALSE], by = "IID")
gp <- fread(file = P$gpca, header = FALSE)
setnames(gp, c("FID","IID", paste0("gPC", seq_len(ncol(gp)-2))))
gp <- gp[, .(IID = trimws(as.character(IID)), gPC1, gPC2)]
meta <- Reduce(function(x,z) merge(x,z,by="IID"), list(dt[,.(IID)], ph, cl, gp))
meta <- meta[complete.cases(meta)]
dt <- dt[IID %in% meta$IID]
setkey(dt,"IID"); setkey(meta,"IID"); dt <- dt[meta$IID]
stopifnot(identical(dt$IID, meta$IID))
RA <- as.matrix(dt[, !c("IID","Metagenome_ID"), with=FALSE]); rownames(RA) <- dt$IID
y <- meta$y
cat(sprintf("discovery: %d individuos, %d especies crudas\n\n", nrow(RA), ncol(RA)))

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
aplicar_mpca <- function(R,pc){
  Z <- as.matrix(compositions::clr(R+1e-6))[,rownames(pc$rotation),drop=FALSE]
  sweep(Z,2,pc$center,"-") %*% pc$rotation[,1:2,drop=FALSE]
}
pbscores <- function(Xc,Bt) as.matrix(compositions::clr(Xc)) %*% Bt

ajustar_mrs <- function(idx, residual){
  Rtr <- RA[idx,,drop=FALSE]
  pr  <- ajustar_prep(Rtr); Xtr <- aplicar_prep(Rtr,pr)
  pc  <- prcomp(as.matrix(compositions::clr(Rtr+1e-6)), center=TRUE, scale.=FALSE)
  cov <- data.table(y=y[idx], meta[idx,CLIN3,with=FALSE], mPC1=pc$x[,1], mPC2=pc$x[,2])
  obj <- if (residual)
    residuals(lm(y ~ Age + gender_dummy + Weight_kg + mPC1 + mPC2, data=cov)) else y[idx]
  bal <- fBalChip_PLS(Xtr, obj, version="cov")
  Bt  <- t(bal$bal)[pr$species,,drop=FALSE]
  PB  <- pbscores(Xtr,Bt)
  as  <- rbindlist(lapply(colnames(PB), function(nm){
    co <- summary(lm(y ~ pb + Age + gender_dummy + Weight_kg + mPC1 + mPC2,
                     data=cbind(cov, pb=PB[,nm])))$coefficients
    if(!"pb" %in% rownames(co)) return(data.table(PB_ID=nm,beta=NA_real_,p=NA_real_))
    data.table(PB_ID=nm, beta=co["pb","Estimate"], p=co["pb","Pr(>|t|)"])
  }))[, FDR := p.adjust(p,method="BH")][]
  sel <- as[!is.na(FDR) & FDR < FDRCUT][order(FDR)]
  list(prep=pr, pc=pc, Bt=Bt, assoc=as, sel=sel)
}
proyectar <- function(fit, idx){
  if(!nrow(fit$sel)) return(rep(0, length(idx)))
  PB <- pbscores(aplicar_prep(RA[idx,,drop=FALSE], fit$prep), fit$Bt)
  as.numeric(PB[, fit$sel$PB_ID, drop=FALSE] %*% fit$sel$beta)
}
especies_de <- function(f){
  if(!nrow(f$sel)) return(character(0))
  B <- f$Bt[,f$sel$PB_ID,drop=FALSE]; rownames(B)[apply(B!=0,1,any)]
}

set.seed(SEED)
q <- cut(y, breaks=quantile(y,probs=seq(0,1,length.out=K+1)), include.lowest=TRUE)
fold <- integer(length(y))
for(lv in levels(q)){ i <- which(q==lv); fold[i] <- sample(rep(1:K,length.out=length(i))) }

Xd <- model.matrix(~ ., data=as.data.frame(meta[,c(CLIN,"gPC1","gPC2"),with=FALSE]))[,-1,drop=FALSE]
colnames(Xd) <- make.names(colnames(Xd), unique=TRUE)
r2_oof <- function(d, f){
  pr <- rep(NA_real_, nrow(d))
  for(k in 1:K){ tr <- fold!=k; pr[!tr] <- predict(lm(f,data=d[tr]), newdata=d[!tr]) }
  1 - sum((d$y-pr)^2)/sum((d$y-mean(d$y))^2)
}

res <- data.table(); oof <- list()
for (modo in c("crudo","residual")) {
  cat(sprintf("=== variante: %s ===\n", modo))
  rsd <- modo == "residual"
  set.seed(SEED)
  full <- ajustar_mrs(seq_len(nrow(RA)), rsd)
  cat(sprintf("  ajuste completo: %d especies, %d balances, %d con FDR<%.2f\n",
              length(full$prep$species), ncol(full$Bt), nrow(full$sel), FDRCUT))
  if (nrow(full$sel))
    print(head(full$sel[, .(PB_ID, beta=signif(beta,3), FDR=signif(FDR,3))], 15), row.names=FALSE)

  m <- rep(NA_real_, nrow(RA)); nsel <- integer(K); esp <- list()
  for(k in 1:K){
    tr <- which(fold!=k); te <- which(fold==k)
    fk <- ajustar_mrs(tr, rsd)
    m[te] <- proyectar(fk, te); nsel[k] <- nrow(fk$sel); esp[[k]] <- especies_de(fk)
    cat(sprintf("  fold %d: %d seleccionados, %d especies\n", k, nsel[k], length(esp[[k]])))
  }
  oof[[modo]] <- m
  jac <- function(a,b) if(length(union(a,b))==0) NA_real_ else length(intersect(a,b))/length(union(a,b))
  J <- outer(1:K,1:K,Vectorize(function(i,j) jac(esp[[i]],esp[[j]])))

  d <- data.table(y=y, as.data.table(Xd))
  v <- setdiff(names(d),"y")
  R0 <- r2_oof(d, as.formula(paste("y ~",paste(v,collapse="+"))))
  if (sd(m, na.rm=TRUE) > 0) {
    d[, MRS := as.numeric(scale(m))]
    R1 <- r2_oof(d, as.formula(paste("y ~",paste(c(v,"MRS"),collapse="+"))))
  } else R1 <- R0
  cat(sprintf("  R2 base %.4f | R2 base+MRS %.4f | dR2 %.5f\n\n", R0, R1, R1-R0))

  res <- rbind(res, data.table(
    variante=modo, n_sel_full=nrow(full$sel), n_sel_folds=paste(nsel,collapse="/"),
    folds_con_senal=sum(nsel>0), jaccard_medio=round(mean(J[upper.tri(J)],na.rm=TRUE),3),
    R2_base=round(R0,4), R2_mas_MRS=round(R1,4), dR2=round(R1-R0,5)))
  saveRDS(list(prep=full$prep, Bt=full$Bt, sel=full$sel, pc_center=full$pc$center,
               pc_rotation=full$pc$rotation[,1:2], modo=modo, semilla=SEED, fecha=Sys.Date()),
          file.path(P$mrs, sprintf("mrs_frozen_%s.rds", modo)))
  fwrite(full$assoc, file.path(P$out, sprintf("balance_association_%s.csv", modo)))
}

cat("=== LAS DOS VARIANTES ===\n"); print(res, row.names=FALSE)
fwrite(res, file.path(P$out,"mrs_crossfit_variants.csv"))
fwrite(data.table(IID=meta$IID, fold=fold,
                  MRS_crudo=oof[["crudo"]], MRS_residual=oof[["residual"]]),
       file.path(P$out,"mrs_oof_variants.csv"))
cat("\nPrimario = 'crudo' (pre-especificado). 'residual' = Enmienda 3, reportar como tal.\n")
cat("guardado en derived/mrs/ y results/\n")
