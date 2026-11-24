# ======================================================================
# 30_reduced_model_holdout.R — Reduced seven-predictor base model in LOS
# #
# Trains the Guangzhou base model (age, sex, weight, height, Asian race indicator, alcohol, regular exercise) in modeling discovery and evaluates it in the internal holdout and subgroups.
# #
# Inputs : clinical and phenotype files
# Outputs: results/reduced_model_holdout.csv
# Paper  : Table S7 (panel D)
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
suppressMessages({library(data.table); library(readxl); library(glmnet)})
source("config/paths.R")
A <- file.path(BASE,"derived/analysis_sets")
disc <- unique(trimws(readLines(file.path(A,"ids_modeling_discovery.txt")))); hold <- unique(trimws(readLines(file.path(A,"ids_internal_holdout.txt"))))
b <- as.data.frame(read_excel(file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx")))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)), y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
cl <- fread(file=file.path(BASE,"data/phenotypes/los_clinical_imputed.csv")); cl[, IID:=trimws(as.character(Subject_ID))]
cl <- unique(cl, by="IID")
si <- function(v) as.integer(trimws(as.character(v)) %in% c("Yes","yes","Y","1","TRUE"))
cl[, `:=`(Alcohol=si(Alcohol_Drinking_cat), Exercise=si(Exercise_Regular_cat), mujer=as.integer(trimws(as.character(gender_dummy)) %in% c("Female","female","F")), Asia=si(Race_Asian))]
RED <- c("Age","Weight_kg","Height_cm","mujer","Asia","Alcohol","Exercise")
m <- merge(ph, cl[, c("IID",RED), with=FALSE], by="IID"); m <- m[complete.cases(m)]
ch <- fread(file=file.path(BASE,"data/phenotypes/guangzhou_clinical.csv"))
agec <- as.numeric(ch$Age); lo <- quantile(agec,0.05,na.rm=TRUE); hi <- quantile(agec,0.95,na.rm=TRUE)
cat(sprintf("ventana de edad Guangzhou (p5-p95): %.1f-%.1f\n", lo, hi))
D <- m[IID %in% disc]; H <- m[IID %in% hold]
cat(sprintf("discovery %d | holdout %d (del total %d)\n", nrow(D), nrow(H), length(hold)))
r2 <- function(o,p) 1-sum((o-p)^2)/sum((o-mean(o))^2)
X <- function(d) as.matrix(d[, RED, with=FALSE])
set.seed(2026); cv <- cv.glmnet(X(D), D$y, alpha=0.25, nfolds=10, standardize=TRUE)
ols <- lm(y ~ ., data=D[, c("y",RED), with=FALSE])
sub <- list("todos"=rep(TRUE,nrow(H)), "mujeres"=H$mujer==1, "mujeres edad Guangzhou"=H$mujer==1 & H$Age>=lo & H$Age<=hi)
set.seed(2026)
res <- rbindlist(lapply(names(sub), function(s){ i <- sub[[s]]; hh <- H[i]
  pe <- as.numeric(predict(cv, newx=X(hh), s="lambda.min")); po <- predict(ols, newdata=hh)
  bt <- replicate(2000, {k <- sample.int(nrow(hh), replace=TRUE); r2(hh$y[k], pe[k])})
  rc <- summary(lm(hh$y ~ pe))$r.squared
  data.table(conjunto="holdout", subgrupo=s, n=nrow(hh), R2_EN=round(r2(hh$y,pe),4), IC_inf=round(quantile(bt,.025),4), IC_sup=round(quantile(bt,.975),4), R2_OLS=round(r2(hh$y,po),4), R2_recal_EN=round(rc,4))}))
print(res, row.names=FALSE)
fwrite(res, file.path(BASE,"results/reduced_model_holdout.csv"))
