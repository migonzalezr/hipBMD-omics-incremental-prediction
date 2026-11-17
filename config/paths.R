# Project root. Run the scripts from the repository root, or set HIPBMD_ROOT
# to the directory that contains data/, derived/, results/, figures/ and external/.
BASE <- Sys.getenv("HIPBMD_ROOT", unset = normalizePath("."))

# Sample-label corrections for the metagenome table, kept with the (private) data:
# data/microbiome/metagenome_id_corrections.csv with columns `old` and `new`.
fix_metagenome_ids <- function(d, col = "Metagenome_ID") {
  d[, (col) := trimws(get(col))]
  f <- file.path(BASE, "data/microbiome/metagenome_id_corrections.csv")
  if (file.exists(f)) {
    fx <- data.table::fread(f)
    m <- match(d[[col]], fx$old)
    d[!is.na(m), (col) := fx$new[m[!is.na(m)]]]
  }
  d[]
}
