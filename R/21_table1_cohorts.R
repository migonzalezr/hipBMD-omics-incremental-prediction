# ======================================================================
# 21_table1_cohorts.R — Table 1
#
# Characteristics of the five analytic cohorts, with analytic sample sizes taken from the prediction files.
#
# Inputs : clinical data; prediction files
# Outputs: results/table1_cohorts.csv; TABLE1.md
# Paper  : Table 1
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)

source("config/paths.R")
AR   <- file.path(BASE,"derived/analysis_sets")
P <- list(
  disc  = file.path(AR,"ids_modeling_discovery.txt"),
  gkeep = file.path(AR,"gwas_keep.txt"),
  gphen = file.path(AR,"gwas_pheno.txt"),
  gcov  = file.path(AR,"gwas_covar.txt"),
  clinB = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  clinL = file.path(BASE,"data/phenotypes/los_lc_clinical_imputed.csv"),
  chcl  = file.path(BASE,"data/phenotypes/guangzhou_clinical.csv"),
  predC = file.path(BASE,"results/guangzhou_predictions.csv"),
  predL = file.path(BASE,"results/lc_predictions.csv"),
  predH = file.path(BASE,"results/holdout_predictions.csv"),
  bmd   = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
  out   = file.path(BASE,"results"),
  docs  = file.path(BASE,"results")
)
falta <- unlist(P)[!file.exists(unlist(P))]
if (length(falta)) stop("no encuentro:\n  ", paste(falta, collapse="\n  "))

CLIN  <- c("Age","gender_dummy","Weight_kg","Height_cm","BMI",
           "Exercise_Regular_cat","Smoke_cat","Sun_Exposure_cat",
           "Alcohol_Drinking_cat","Cheese_Eat_cat","Milk_Drink_cat","Yogurt_Eat_cat")
RAZAS <- c("Race_Caucasian","Race_African","Race_Asian","Race_Other")
DASH  <- "\u2014"

b  <- as.data.frame(read_excel(P$bmd))
ph <- unique(data.table(IID=trimws(as.character(b$Sample_ID)),
                        y=suppressWarnings(as.numeric(b$HTOT_BMD)))[!is.na(y)], by="IID")
clB <- fread(file=P$clinB); clB[, IID := trimws(as.character(Subject_ID))]
cols_raza <- intersect(RAZAS, names(clB))
clB <- unique(clB[, c("IID",CLIN,cols_raza), with=FALSE], by="IID")
clL <- fread(file=P$clinL); clL[, IID := trimws(as.character(Subject_ID))]
clL <- unique(clL[, c("IID", intersect(c(CLIN,cols_raza), names(clL))), with=FALSE], by="IID")
chC <- fread(file=P$chcl); chC[, IID := trimws(as.character(ID))]

ids_disc <- trimws(readLines(P$disc))
ids_hold <- fread(file=P$predH)$IID
ids_lc   <- fread(file=P$predL)$IID
ids_ch   <- fread(file=P$predC)$IID

gk <- fread(file=P$gkeep); setnames(gk, sub("^#","",names(gk)))
gp <- fread(file=P$gphen); setnames(gp, sub("^#","",names(gp)))
gc_<- fread(file=P$gcov ); setnames(gc_, sub("^#","",names(gc_)))
gw <- Reduce(function(x,z) merge(x,z,by="IID"),
             list(data.table(IID=trimws(as.character(gk$IID))),
                  gp[, .(IID=trimws(as.character(IID)), y=as.numeric(HTOT_BMD))],
                  gc_[, .(IID=trimws(as.character(IID)), Age=as.numeric(Age),
                          Sex=as.integer(Sex), Weight_kg=as.numeric(Weight))]))
gw <- gw[complete.cases(gw)]
cat(sprintf("GWAS: %d elegibles en keep -> %d con fenotipo y covariables completas\n",
            nrow(gk), nrow(gw)))
cat(sprintf("  (debe coincidir con el OBS_CT que reporto plink; si no, revisar)\n"))

solapa <- length(intersect(gw$IID, c(ids_disc, ids_hold, ids_lc)))
cat(sprintf("  solapamiento GWAS vs cohortes de modelado: %d (debe ser 0)\n", solapa))
if (solapa > 0) stop("el GWAS comparte individuos con las cohortes de modelado")

sn  <- function(v) sprintf("%.2f \u00b1 %.2f", mean(v,na.rm=TRUE), sd(v,na.rm=TRUE))

sn3 <- function(v) sprintf("%.3f \u00b1 %.3f", mean(v,na.rm=TRUE), sd(v,na.rm=TRUE))
pc <- function(v, si=c("Yes","yes","1","TRUE")) {
  v <- trimws(as.character(v)); ok <- !is.na(v) & v!="" & v!="NA"
  if (!sum(ok)) return(DASH)
  sprintf("%d (%.0f%%)", sum(v %in% si), 100*sum(v %in% si)/sum(ok))
}
FILAS <- c("n","Total hip BMD (g/cm2)","Age (years)","Female",
           "Weight (kg)","Height (cm)","BMI (kg/m2)","Regular exercise","Smoking",
           "Sun exposure","Alcohol","Cheese","Milk","Yogurt",
           "Caucasian","African American","Asian","Other")

col_los <- function(ids, cl, extra=NULL) {
  d <- merge(merge(data.table(IID=trimws(as.character(ids))), ph, by="IID"), cl, by="IID")
  v <- c(as.character(nrow(d)), sn3(d$y), sn(d$Age),
         pc(d$gender_dummy, c("Female","female","F")),
         sn(d$Weight_kg), sn(d$Height_cm), sn(d$BMI),
         pc(d$Exercise_Regular_cat), pc(d$Smoke_cat), pc(d$Sun_Exposure_cat),
         pc(d$Alcohol_Drinking_cat), pc(d$Cheese_Eat_cat),
         pc(d$Milk_Drink_cat), pc(d$Yogurt_Eat_cat))
  for (rz in RAZAS)
    v <- c(v, if (rz %in% names(d)) pc(d[[rz]], c("1","TRUE","Yes")) else DASH)
  v
}
col_gwas <- function(d) {
  cl <- merge(d[, .(IID)], clB, by="IID", all.x=TRUE)
  setkey(cl,"IID"); setkey(d,"IID"); cl <- cl[d$IID]
  v <- c(as.character(nrow(d)), sn3(d$y), sn(d$Age),
         sprintf("%d (%.0f%%)", sum(d$Sex==2), 100*mean(d$Sex==2)),
         sn(d$Weight_kg), sn(cl$Height_cm), sn(cl$BMI),
         pc(cl$Exercise_Regular_cat), pc(cl$Smoke_cat), pc(cl$Sun_Exposure_cat),
         pc(cl$Alcohol_Drinking_cat), pc(cl$Cheese_Eat_cat),
         pc(cl$Milk_Drink_cat), pc(cl$Yogurt_Eat_cat))
  for (rz in RAZAS) v <- c(v, if (rz %in% names(cl)) pc(cl[[rz]], c("1","TRUE","Yes")) else DASH)
  v
}
col_china <- function(ids) {
  d <- merge(data.table(IID=trimws(as.character(ids))), chC, by="IID")
  d[, y := as.numeric(HTOT_BMD)]
  v <- c(as.character(nrow(d)), sn3(d$y), sn(d$Age), sprintf("%d (100%%)", nrow(d)),
         sn(d$Weight_kg), sn(d$Height_cm), sn(d$BMI),
         pc(d$Exercise_cat,"1"), pc(d$Smoke_cat,"1"), DASH,
         pc(d$Alcohol_cat,"1"), DASH, pc(d$Milk_cat,"1"), pc(d$Yogurt_cat,"1"))
  c(v, rep(DASH, length(RAZAS)))
}

tt <- data.table(
  Characteristic                   = FILAS,
  `GWAS discovery`                 = col_gwas(gw),
  `Modeling discovery (BGI)`       = col_los(ids_disc, clB),
  `Internal holdout (BGI)`         = col_los(ids_hold, clB),
  `Validation 1: LC Sciences`      = col_los(ids_lc,   clL),
  `Validation 2: Guangzhou`        = col_china(ids_ch))

cat("\n=== TABLE 1 \u2014 analytic cohorts ===\n"); print(tt, row.names=FALSE)
fwrite(tt, file.path(P$out,"table1_cohorts.csv"))

md <- c("# Table 1. Characteristics of the analytic cohorts", "",
        "Values are mean \u00b1 SD; lifestyle variables are n (%) answering yes.", "",
        "**Denominator.** Percentages are computed over **non-missing** responses, not over the",
        "full cohort. This differs from the convention used in the originally submitted Table 1",
        "and explains, for example, why milk consumption in Guangzhou reads 37% here and 29%",
        sprintf("there. %s = variable not collected in that cohort.", DASH), "",
        "**Two values in Guangzhou are constants, not errors.** No participant reported smoking",
        "(0%), and every participant who answered the yogurt question reported some consumption",
        "(100%), because the response option \"rarely\" counts as yes. Both variables were",
        "therefore uninformative in that cohort and were excluded from the reduced model.", "",
        paste0("| ", paste(names(tt), collapse=" | "), " |"),
        paste0("|", paste(rep("---", ncol(tt)), collapse="|"), "|"),
        apply(tt, 1, function(r) paste0("| ", paste(r, collapse=" | "), " |")), "",
        "**Disjointness.** The GWAS discovery cohort shares no participant with any modeling",
        "cohort, by design: it comprises LOS participants without metagenomic data, so PRS",
        "weights are never estimated in the individuals used to evaluate the score. The script",
        "verifies this overlap is zero and halts otherwise. Sample sizes are analytic",
        "(individuals actually entering each model), not recruitment counts.", "",
        "**The GWAS cohort is not demographically representative of the modeling cohort.**",
        "It is 43% female versus 66%, 47% African American versus 28%, and 5.1 years younger,",
        "with mean BMD 0.05 g/cm2 higher. This follows from the metagenomic subsample not",
        "being drawn at random, and it means PRS weights are estimated under a different",
        "ancestry mix than the one they are applied to. Because GWAS effect sizes depend on",
        "ancestry-specific linkage disequilibrium, an ancestry mismatch exists within LOS,",
        "not only between LOS and the external cohorts.")
writeLines(md, file.path(P$docs,"TABLE1.md"))
cat("\nguardado results/table1_cohorts.csv y results/TABLE1.md\n")
cat("===== FIN OK =====\n")

cs <- readRDS(file.path(BASE,"derived/analysis_sets/cohort_final.rds"))
str(cs); print(names(cs))

print(sapply(cs, function(x) if (is.character(x)||is.factor(x)) length(unique(x)) else NA))
