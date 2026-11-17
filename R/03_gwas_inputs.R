# ======================================================================
# 03_gwas_inputs.R — GWAS inputs
#
# Writes the keep, phenotype and covariate files (age, sex, weight, PC1-PC10) for the independent LOS GWAS.
#
# Inputs : derived/analysis_sets/ids_gwas_sample.txt; data/genotypes/los_qc_pca.eigenvec; data/phenotypes/los_clinical_imputed.csv
# Outputs: derived/analysis_sets/gwas_keep.txt, gwas_pheno.txt, gwas_covar.txt
# Paper  : Experimental Section (GWAS)
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)

source("config/paths.R")
P <- list(
  fam        = file.path(BASE,"data/genotypes/los_qc.fam"),
  mgs_ids    = file.path(BASE,"data/microbiome/metagenome_sample_ids.xlsx"),
  mgs_filter = file.path(BASE,"data/microbiome/metagenome_sample_filter.xlsx"),
  bmd        = file.path(BASE,"data/phenotypes/los_dxa_bmd.xlsx"),
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

ids <- trimws(readLines(file.path(P$out, "ids_gwas_sample.txt")))
cat("IID del GWAS discovery:", length(ids), "\n")

pca <- fread(P$pca, header = FALSE)
setnames(pca, c("FID","IID", paste0("PC", seq_len(ncol(pca) - 2))))
pca[, IID := trimws(as.character(IID))][, FID := NULL]

cl <- fread(P$clin)
cl[, Subject_ID := trimws(as.character(Subject_ID))]
cl <- unique(cl[, .(IID = Subject_ID, Age, Sex = gender_dummy, Weight = Weight_kg)],
             by = "IID")

ph <- leer_bmd(); setnames(ph, "BMD", "HTOT_BMD")

d <- Reduce(function(x, y) merge(x, y, by = "IID", all.x = TRUE),
            list(data.table(IID = ids), ph, cl, pca))
pcs <- paste0("PC", 1:10)
cat("\n--- faltantes ---\n")
for (v in c("HTOT_BMD","Age","Sex","Weight","PC1","PC10"))
  cat(sprintf("  %-10s %d\n", v, sum(is.na(d[[v]]))))

d <- d[complete.cases(d[, c("HTOT_BMD","Age","Sex","Weight", pcs), with = FALSE])]
cat("\n>>> N FINAL DEL GWAS:", nrow(d), "de", length(ids), "<<<\n")

d[, Sex := fifelse(tolower(trimws(as.character(Sex))) %in% c("male","m","1"),   1L,
           fifelse(tolower(trimws(as.character(Sex))) %in% c("female","f","2"), 2L,
                   NA_integer_))]
stopifnot(!any(is.na(d$Sex)))
cat("sexo (1=Male, 2=Female):\n"); print(table(d$Sex))

d[, FID := 0L]
keep <- d[, .(FID, IID)];           setnames(keep, "FID", "#FID")
phe  <- d[, .(FID, IID, HTOT_BMD)]; setnames(phe,  "FID", "#FID")
cov  <- d[, c("FID","IID","Age","Sex","Weight", pcs), with = FALSE]
setnames(cov, "FID", "#FID")
stopifnot(!any(sapply(cov, is.character)))

fwrite(keep, file.path(P$out,"gwas_keep.txt"),  sep = "\t")
fwrite(phe,  file.path(P$out,"gwas_pheno.txt"), sep = "\t")
fwrite(cov,  file.path(P$out,"gwas_covar.txt"), sep = "\t")
cat("\nescritos:", nrow(keep), "filas |", ncol(cov), "columnas en covar\n")
