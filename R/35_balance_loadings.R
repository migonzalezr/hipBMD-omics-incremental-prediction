# ======================================================================
# 35_balance_loadings.R — Species composition of principal balances 1 and 136
#
# Loadings, side (numerator/denominator) and detection prevalence in the modeling discovery set of the
# species that define principal balances 1 and 136 of the frozen balance basis.
#
# Inputs : derived/mrs/mrs_frozen_crudo.rds, species abundance table, metagenome ID map, discovery IDs
# Outputs: results/balance1_balance136_species_loadings.csv
# Paper  : Figure S8
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table); library(readxl)
source("config/paths.R")
P <- list(
  mgs  = file.path(BASE,"data/microbiome/los_bgi_species_abundance.csv"),
  mast = file.path(BASE,"data/microbiome/metagenome_id_map.xlsx"),
  disc = file.path(BASE,"derived/analysis_sets/ids_modeling_discovery.txt"),
  rds  = file.path(BASE,"derived/mrs/mrs_frozen_crudo.rds"),
  out  = file.path(BASE,"results")
)
disc_id <- unique(trimws(readLines(P$disc)))
mgs <- fread(file = P$mgs); setnames(mgs, 1, "Metagenome_ID")
mgs <- fix_metagenome_ids(mgs)
pte <- as.data.table(read_excel(P$mast, sheet = "bgi"))
setnames(pte, "BGI_1794", "Subject_ID")
pte <- pte[, .(Metagenome_ID=trimws(as.character(Metagenome_ID)), IID=trimws(as.character(Subject_ID)))]
dt <- merge(pte[IID %in% disc_id], mgs, by="Metagenome_ID")
dt <- dt[IID %in% disc_id]
stopifnot(nrow(dt) == length(disc_id))
fr <- readRDS(P$rds)
part <- function(tax, tag) sub(paste0("^.*", tag, "__([^.]+).*$"), "\\1", tax)
rows <- rbindlist(lapply(c(1L, 136L), function(b) {
  v <- fr$Bt[, paste0("bal", b)]; v <- v[v != 0]
  ab <- as.matrix(dt[, names(v), with = FALSE])
  data.table(balance = b, species = part(names(v), "s"), kingdom = part(names(v), "k"),
             phylum = part(names(v), "p"), side = ifelse(v < 0, "Denominator", "Numerator"),
             loading = unname(v), n_detected = colSums(ab > 0), n_total = nrow(ab),
             prevalence_percent = 100 * colSums(ab > 0) / nrow(ab), taxonomy = names(v))
}))
fwrite(rows, file.path(P$out, "balance1_balance136_species_loadings.csv"))
cat("guardado en results/balance1_balance136_species_loadings.csv\n")
