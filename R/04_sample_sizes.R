# ======================================================================
# 04_sample_sizes.R — Sample-size cascade
#
# Rebuilds the participant flow from primary sources with consistency checks at every step.
#
# Inputs : primary data and frozen ID lists
# Outputs: derived/analysis_sets/sample_sizes.csv
# Paper  : Figure 1a
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================

rm(list = ls()); gc()
library(data.table); library(readxl)

source("config/paths.R")
P <- list(
  fam  = file.path(BASE, "data/genotypes/los_qc.fam"),
  bmd  = file.path(BASE, "data/phenotypes/los_hip_bmd.xlsx"),
  king = file.path(BASE, "derived/kinship/los_king.kin0"),
  ar   = file.path(BASE, "derived/analysis_sets")
)
N_GUANGZHOU <- 499

lee <- function(f) unique(trimws(readLines(file.path(P$ar, f))))

fam <- fread(file = P$fam, header = FALSE, colClasses = "character")
universo <- unique(trimws(fam$V2))

bmd <- as.data.table(read_excel(P$bmd))
bmd[, Sample_ID := trimws(as.character(Sample_ID))]
con_bmd <- intersect(universo, bmd[!is.na(HTOT_BMD), Sample_ID])

ev   <- lee("ids_evaluation_set.txt")
lcc  <- lee("ids_lc_subset.txt")
gdis <- lee("ids_gwas_sample.txt")
bdis <- lee("ids_modeling_discovery.txt")
bhol <- lee("ids_internal_holdout.txt")
keep <- fread(file = file.path(P$ar, "gwas_keep.txt"), colClasses = "character")
ganal <- unique(trimws(keep[[2]]))

bgi  <- setdiff(ev, lcc)
cand <- setdiff(con_bmd, ev)

stopifnot(
  length(universo) == 6732,
  all(ev %in% con_bmd),
  length(cand) + length(ev) == length(con_bmd),
  all(ganal %in% gdis),
  all(bdis %in% bgi), all(bhol %in% bgi),
  length(intersect(bdis, bhol)) == 0,
  length(intersect(bgi, lcc))  == 0,
  length(intersect(ganal, ev)) == 0
)

kin <- fread(file = P$king)
setnames(kin, sub("^#", "", names(kin)))
en_par <- unique(c(trimws(kin$IID1), trimws(kin$IID2)))
sin_cov <- setdiff(gdis, ganal)
cat("excluidos del GWAS por covariables incompletas:", length(sin_cov), "\n")
cat("  de ellos, con algun par de parentesco:", sum(sin_cov %in% en_par),
    " (debe ser 0)\n")
stopifnot(sum(sin_cov %in% en_par) == 0)

podados_bgi <- setdiff(bgi, union(bdis, bhol))
cat("podados de BGI por parentesco con LC:", length(podados_bgi), "\n")
cat("  de ellos, con algun par de parentesco:", sum(podados_bgi %in% en_par),
    " (debe ser", length(podados_bgi), ")\n")
stopifnot(all(podados_bgi %in% en_par))

tab <- data.table(
  nivel = c(0,1,2,3,4,2,3,4,5,5,3,0),
  bloque = c(
    "Universo post-QC (los_qc.fam)",
    "Con HTOT_BMD",
    "Sin metagenomica: candidatos al GWAS",
    "  tras poda de parentesco (KING)",
    "  analizados en el GWAS",
    "Con metagenomica: evaluacion",
    "BGI",
    "  tras poda de parentesco con LC",
    "    BGI discovery",
    "    BGI holdout",
    "LC",
    "Guangzhou (cohorte externa)"),
  block_en = c(
    "Post-QC genotyped universe",
    "With total hip BMD",
    "Without metagenomics: GWAS candidates",
    "  after kinship pruning (KING)",
    "  analysed in the GWAS",
    "With metagenomics: evaluation set",
    "BGI platform",
    "  after kinship pruning against LC",
    "    BGI discovery",
    "    BGI holdout",
    "LC platform",
    "Guangzhou (external cohort)"),
  n = c(
    length(universo), length(con_bmd), length(cand), length(gdis), length(ganal),
    length(ev), length(bgi), length(union(bdis, bhol)),
    length(bdis), length(bhol), length(lcc), N_GUANGZHOU),
  derivacion = c(
    "filas de los_qc.fam",
    "universo con HTOT_BMD no missing",
    "con BMD menos el conjunto de evaluacion",
    sprintf("-%d por parentesco de tercer grado o mas cercano", length(cand) - length(gdis)),
    sprintf("-%d sin peso/talla/IMC (outliers antropometricos)", length(sin_cov)),
    "con genotipo, metagenomica y BMD",
    "evaluacion menos LC",
    sprintf("-%d en pares de parentesco con LC", length(podados_bgi)),
    "particion ~80% de BGI",
    "particion ~20% de BGI; se abre una sola vez",
    "segunda plataforma de secuenciacion, misma cohorte LOS",
    "constante declarada; no deriva de estos archivos")
)

stopifnot(
  tab[bloque == "    BGI discovery", n] + tab[bloque == "    BGI holdout", n] ==
    tab[bloque == "  tras poda de parentesco con LC", n],
  tab[bloque == "BGI", n] + tab[bloque == "LC", n] ==
    tab[bloque == "Con metagenomica: evaluacion", n]
)

out <- file.path(P$ar, "sample_sizes.csv")
if (file.exists(out))
  file.rename(out, file.path(P$ar, sprintf("sample_sizes_previous_%s.csv", Sys.Date())))
fwrite(tab, out)

cat("\n=== CASCADA DE TAMANOS MUESTRALES ===\n")
print(tab[, .(bloque, n, derivacion)], row.names = FALSE)
cat("\nguardado en derived/analysis_sets/sample_sizes.csv\n")
cat("todos los chequeos pasaron\n")
