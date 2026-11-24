# ======================================================================
# 33_lookup_li2025.R — Lookup of the species reported by Li et al. (2025)
# #
# Association of the reported species with hip BMD in the 428-species universe, before and after adjustment.
# #
# Inputs : results/species_association.csv
# Outputs: results/lookup_li2025.csv
# Paper  : Table S14
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
rm(list = ls())
library(data.table)
source("config/paths.R")
RES  <- file.path(BASE, "results")
esp  <- fread(file = file.path(RES, "species_association.csv"))

LI <- data.table(
  especie_li = c("Lachnospira eligens", "Blautia wexlerae", "Lachnospiraceae bacterium AM48_27BH",
                 "Clostridium sp. AF20_17LB", "SGB4664", "Firmicutes bacterium AF16_15",
                 "Ruminococcus torques", "Dorea longicatena", "Bacteroides vulgatus"),
  fuente     = c(rep("Li 2025, FDR < 0.05", 8), "Lin 2023"),
  signo_li   = c("+", "+", "+", "+", "+", "+", "-", "-", "-"),
  patron     = c("eligens", "wexlerae", "AM48", "AF20", "SGB4664", "AF16",
                 "Ruminococcus_torques", "Dorea_longicatena", "Bacteroides_vulgatus"))
out <- rbindlist(lapply(seq_len(nrow(LI)), function(i) {
  m <- esp[grepl(LI$patron[i], especie, ignore.case = TRUE)]
  if (!nrow(m)) return(cbind(LI[i, .(especie_li, fuente, signo_li)],
                             data.table(especie_nuestra = "not in the 428-species universe")))
  cbind(LI[i, .(especie_li, fuente, signo_li)],
        m[, .(especie_nuestra = especie, beta_crudo = signif(beta_crudo, 3), FDR_crudo = signif(FDR_crudo, 3),
              beta_ajust = signif(beta_ajust, 3), p_ajust = signif(p_ajust, 3), FDR_ajust = signif(FDR_ajust, 3))])
}), fill = TRUE)
out[!is.na(beta_crudo), signo_nuestro_crudo := fifelse(beta_crudo > 0, "+", "-")]
print(out, row.names = FALSE)
fwrite(out, file.path(RES, "lookup_li2025_species.csv"))
cat("\nguardado en results/lookup_li2025_species.csv\n===== FIN OK =====\n")
