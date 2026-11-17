#!/usr/bin/env bash
# KING-robust kinship among QC'd LOS participants (pairs with kinship >= 0.0442, third degree or closer).
# Output is read by R/02_kinship_firewall.R.
set -euo pipefail
PLINK2=${PLINK2:-plink2}
GENO=data/genotypes/los_qc          # PLINK binary fileset (.bed/.bim/.fam) after genotype QC
OUT=derived/kinship; mkdir -p "$OUT"
$PLINK2 --bfile "$GENO" --autosome --maf 0.05 --geno 0.02 \
        --make-king-table --king-table-filter 0.0442 \
        --threads 8 --out "$OUT/los_king"
