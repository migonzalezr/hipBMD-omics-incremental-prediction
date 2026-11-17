# ======================================================================
# 01_cohort_census.R — Cohort census
#
# Assigns every genotyped LOS participant to the GWAS pool or the evaluation set according to availability of gut metagenomes and hip BMD.
#
# Inputs : data/genotypes/los_qc.fam; data/phenotypes/los_hip_bmd.xlsx; data/microbiome/metagenome_sample_ids.xlsx; data/microbiome/metagenome_sample_filter.xlsx
# Outputs: derived/analysis_sets/cohort_census.rds
# Paper  : Results 2.1; Figure 1a
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)

source("config/paths.R")
P <- list(
  fam        = file.path(BASE,"data/genotypes/los_qc.fam"),
  mgs_ids    = file.path(BASE,"data/microbiome/metagenome_sample_ids.xlsx"),
  mgs_filter = file.path(BASE,"data/microbiome/metagenome_sample_filter.xlsx"),
  bmd        = file.path(BASE,"data/phenotypes/los_hip_bmd.xlsx"),
  pca        = file.path(BASE,"data/genotypes/los_qc_pca.eigenvec"),
  clin       = file.path(BASE,"data/phenotypes/los_clinical_imputed.csv"),
  kin        = file.path(BASE,"derived/kinship/los_king.kin0"),
  out        = file.path(BASE,"derived/analysis_sets")
)
dir.create(P$out, showWarnings = FALSE, recursive = TRUE)

yes <- function(x) toupper(trimws(as.character(x))) %in% c("YES","Y","TRUE","1")

leer_plataforma <- function(sheet) {
  d <- as.data.frame(read_excel(P$mgs_ids, sheet = sheet))[, 1:2]
  unique(trimws(as.character(d[[2]][!is.na(d[[2]]) & d[[2]] != ""])))
}

leer_bmd <- function() {
  b <- as.data.frame(read_excel(P$bmd))
  d <- data.table(IID = trimws(as.character(b$Sample_ID)),
                  BMD = suppressWarnings(as.numeric(b$HTOT_BMD)))
  unique(d[!is.na(BMD)], by = "IID")
}

censo <- function() {
  gid     <- as.character(fread(P$fam, header = FALSE)$V2)
  bgi_ids <- leer_plataforma("bgi")
  lc_ids  <- leer_plataforma("lc")
  bmd     <- leer_bmd()
  mgs_any <- intersect(union(bgi_ids, lc_ids), gid)
  evalset <- intersect(mgs_any, bmd$IID)
  list(gid = gid, bgi_ids = bgi_ids, lc_ids = lc_ids, bmd = bmd,
       mgs_any = mgs_any, evalset = evalset,
       gwas0   = intersect(setdiff(gid, mgs_any), bmd$IID),
       BGI     = intersect(bgi_ids, evalset),
       LC      = setdiff(evalset, intersect(bgi_ids, evalset)))
}

C <- censo()
cat("universo (los_qc.fam) :", length(C$gid), "\n")
cat("con metagenomica         :", length(C$mgs_any), "| con BMD:", length(C$evalset), "\n")
cat("solo WGS                 :", length(setdiff(C$gid, C$mgs_any)),
    "| con BMD:", length(C$gwas0), "\n")
cat("BGI                      :", length(C$BGI), "| LC:", length(C$LC), "\n")
stopifnot(length(C$mgs_any) + length(setdiff(C$gid, C$mgs_any)) == length(C$gid))

flt <- as.data.frame(read_excel(P$mgs_filter, sheet = "sample_sizes"))[, 1:5]
names(flt) <- c("LOS_ID","excl_gwas","bmd_flag","mgs_flag","all3_flag")
flt$LOS_ID <- trimws(as.character(flt$LOS_ID))
flt <- flt[!is.na(flt$LOS_ID) & flt$LOS_ID != "", ]

cat("\n--- cascada de exclusiones (tabla maestra) ---\n")
print(table(trimws(as.character(flt$excl_gwas)), useNA = "ifany"))
flag_mgs <- intersect(flt$LOS_ID[yes(flt$mgs_flag)], C$gid)
cat("\nquien tiene MGS, por dos caminos independientes:\n")
cat("  mapeo de plataforma:", length(C$mgs_any),
    "| bandera mgs_flag:", length(flag_mgs),
    "| diferencias:", length(setdiff(C$mgs_any, flag_mgs)),
                     length(setdiff(flag_mgs, C$mgs_any)), "\n")

cat("\n", length(C$gwas0), "+", length(C$evalset), "=",
    length(C$gwas0) + length(C$evalset),
    " (debe ser el n de participantes con HTOT_BMD entre los",
    length(C$gid), ")\n")

saveRDS(C, file.path(P$out, "cohort_census.rds"))
