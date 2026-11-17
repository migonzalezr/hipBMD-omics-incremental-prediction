# ======================================================================
# 02_kinship_firewall.R — Kinship firewall
#
# Uses KING kinship (>= 0.0442) to remove relatives across analysis boundaries (GWAS vs evaluation, BGI vs LC) and splits BGI participants 80:20 into modeling discovery and internal holdout by family components (seed 2026).
#
# Inputs : derived/kinship/los_king.kin0 (plink/01_king.sh); cohort census inputs
# Outputs: derived/analysis_sets/ids_*.txt; cohort_final.rds
# Paper  : Results 2.1; Figure 1a; Experimental Section
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

C <- censo()

k <- fread(P$kin); setnames(k, sub("^#", "", names(k)))
i1 <- grep("^(IID1|ID1)$", names(k), value = TRUE)[1]
i2 <- grep("^(IID2|ID2)$", names(k), value = TRUE)[1]
kc <- grep("KINSHIP", names(k), ignore.case = TRUE, value = TRUE)[1]
k  <- data.table(A = as.character(k[[i1]]), B = as.character(k[[i2]]),
                 PHI = as.numeric(k[[kc]]))
cat("pares con phi >= 0.0442:", nrow(k), "\n")
print(table(cut(k$PHI, c(0.0442,0.0884,0.177,0.354,Inf),
                labels = c("3er","2do","1er","dup/MZ"), right = FALSE)))

cross <- k[(A %in% C$gwas0 & B %in% C$evalset) | (B %in% C$gwas0 & A %in% C$evalset)]
drop_gwas <- unique(c(cross$A[cross$A %in% C$gwas0], cross$B[cross$B %in% C$gwas0]))
gwas <- setdiff(C$gwas0, drop_gwas)
cat("\n(a) GWAS<->eval:", nrow(cross), "pares |", length(drop_gwas),
    "removidos |", length(C$gwas0), "->", length(gwas), "\n")

bl <- k[(A %in% C$BGI & B %in% C$LC) | (B %in% C$BGI & A %in% C$LC)]
drop_bgi <- unique(c(bl$A[bl$A %in% C$BGI], bl$B[bl$B %in% C$BGI]))
BGI <- setdiff(C$BGI, drop_bgi)
cat("(b) BGI<->LC   :", nrow(bl), "pares |", length(drop_bgi),
    "removidos de BGI | LC intacta:", length(C$LC), "\n")

intra <- k[A %in% BGI & B %in% BGI]
comp_of <- function(ids, a, b) {
  comp <- setNames(seq_along(ids), ids)
  repeat {
    cambio <- FALSE
    for (i in seq_along(a)) {
      x <- a[i]; y <- b[i]
      if (comp[[x]] != comp[[y]]) {
        m <- min(comp[[x]], comp[[y]]); o <- max(comp[[x]], comp[[y]])
        comp[comp == o] <- m; cambio <- TRUE
      }
    }
    if (!cambio) break
  }
  comp
}
comp   <- comp_of(BGI, intra$A, intra$B)
grupos <- split(names(comp), comp)
cat("componentes familiares:", length(grupos),
    "| mayor:", max(lengths(grupos)), "\n")

set.seed(2026)
grupos   <- grupos[sample(length(grupos))]
objetivo <- round(0.80 * length(BGI))
disc <- character(0)
for (g in grupos) if (length(disc) < objetivo) disc <- c(disc, g)
hold <- setdiff(BGI, disc)

chk <- function(x, y) nrow(k[(A %in% x & B %in% y) | (B %in% x & A %in% y)])
cat("\nsplit: discovery", length(disc), "| holdout", length(hold), "\n")
cat("cruces (todos deben ser 0) -> GWAS<->eval:", chk(gwas, C$evalset),
    "| disc<->hold:", chk(disc, hold), "| disc<->LC:", chk(disc, C$LC), "\n")
stopifnot(chk(gwas, C$evalset) == 0, chk(disc, hold) == 0, chk(disc, C$LC) == 0,
          length(intersect(disc, hold)) == 0)

bd <- setNames(C$bmd$BMD, C$bmd$IID)
cat("BMD media -> discovery", round(mean(bd[disc]), 4),
    "| holdout", round(mean(bd[hold]), 4), "\n")

writeLines(gwas,      file.path(P$out,"ids_gwas_sample.txt"))
writeLines(C$evalset, file.path(P$out,"ids_evaluation_set.txt"))
writeLines(disc,      file.path(P$out,"ids_modeling_discovery.txt"))
writeLines(hold,      file.path(P$out,"ids_internal_holdout.txt"))
writeLines(C$LC,      file.path(P$out,"ids_lc_subset.txt"))
saveRDS(list(gwas0 = C$gwas0, gwas = gwas, drop_gwas = drop_gwas,
             drop_bgi = drop_bgi, evalset = C$evalset, BGI = BGI, LC = C$LC,
             discovery = disc, holdout = hold, comp = comp, kin = k,
             semilla = 2026, fecha = Sys.time()),
        file.path(P$out,"cohort_final.rds"))
cat("\nLISTAS CONGELADAS. No se vuelven a tocar.\n")
