#!/usr/bin/env bash
# Independent LOS GWAS of total hip BMD in the GWAS sample only (participants without metagenomes,
# unrelated to the evaluation set). Covariates: age, sex, weight, PC1-PC10 (R/03_gwas_inputs.R).
set -euo pipefail
PLINK2=${PLINK2:-plink2}
GENO=data/genotypes/los_qc
SETS=derived/analysis_sets
OUT=derived/gwas; mkdir -p "$OUT"
$PLINK2 --bfile "$GENO" \
  --keep "$SETS/gwas_keep.txt" \
  --pheno "$SETS/gwas_pheno.txt" --pheno-name HTOT_BMD \
  --covar "$SETS/gwas_covar.txt" --covar-variance-standardize \
  --glm hide-covar cols=+a1freq,+beta,+nobs \
  --autosome --maf 0.01 --geno 0.02 --hwe 1e-6 midp \
  --threads 8 --out "$OUT/los_gwas"
