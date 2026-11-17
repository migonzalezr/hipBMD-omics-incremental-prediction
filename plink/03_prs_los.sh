#!/usr/bin/env bash
# PRS from the LOS GWAS by clumping and thresholding (LD from the GWAS sample; r2 = 0.1, 250 kb),
# scored at nine P-value thresholds in all LOS participants. Threshold selection: R/05a_prs_threshold_los.R.
set -euo pipefail
PLINK2=${PLINK2:-plink2}
GENO=data/genotypes/los_qc
SETS=derived/analysis_sets
GWAS=derived/gwas/los_gwas.HTOT_BMD.glm.linear     # columns: 3=ID 7=A1 12=BETA 15=P
OUT=derived/prs_los; mkdir -p "$OUT"

{ printf "ID\tA1\tBETA\n"; awk 'NR>1 {print $3"\t"$7"\t"$12}' "$GWAS"; } > "$OUT/score_weights.txt"
{ printf "ID\tP\n";        awk 'NR>1 {print $3"\t"$15}' "$GWAS"; }       > "$OUT/snp_pvalues.txt"
cat > "$OUT/range_list.txt" <<'RANGES'
PT_5e-8 0 5e-8
PT_1e-6 0 1e-6
PT_1e-4 0 1e-4
PT_1e-3 0 1e-3
PT_0.01 0 0.01
PT_0.05 0 0.05
PT_0.1 0 0.1
PT_0.5 0 0.5
PT_1.0 0 1.0
RANGES

$PLINK2 --bfile "$GENO" --keep "$SETS/gwas_keep.txt" --clump "$GWAS" \
        --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --threads 8 --out "$OUT/los_clumped"
CL=$(ls "$OUT"/los_clumped.clump* | head -1)
awk 'NR>1 && $3!="" && $3!="NA" {print $3}' "$CL" | sort -u > "$OUT/clumped_snps.txt"

$PLINK2 --bfile "$GENO" --extract "$OUT/clumped_snps.txt" \
        --score "$OUT/score_weights.txt" 1 2 3 header cols=+scoresums \
        --q-score-range "$OUT/range_list.txt" "$OUT/snp_pvalues.txt" 1 2 header \
        --threads 8 --out "$OUT/los_prs"
