#!/usr/bin/env bash
# PRS from external UK Biobank total femur BMD summary statistics (Qian et al. 2025; NGDC GVM GVP000082).
# Harmonized by rsID; biallelic SNVs only; palindromic (A/T, C/G) and duplicated rsIDs removed.
# Clumped with LD from the LOS GWAS sample and scored at nine thresholds. Selection: R/05b_prs_threshold_external.R.
set -eu
PLINK2=${PLINK2:-plink2}
GENO=data/genotypes/los_qc
SETS=derived/analysis_sets
SUMSTATS=data/external/ukb_total_femur_bmd_sumstats.tsv.gz   # columns: 3=effect allele 4=other allele 5=beta 8=P 9=rsID
OUT=derived/prs_external; mkdir -p "$OUT"
READER=cat; case "$SUMSTATS" in *.gz) READER=zcat;; esac

$READER "$SUMSTATS" | awk -F'\t' 'BEGIN{OFS="\t"; print "ID","A1","BETA","P"} NR==1{next}
  $9 ~ /^rs/ && length($3)==1 && length($4)==1 {
    amb=($3=="A"&&$4=="T")||($3=="T"&&$4=="A")||($3=="C"&&$4=="G")||($3=="G"&&$4=="C");
    if(!amb) print $9,$3,$5,$8 }' > "$OUT/external_all.txt"
awk 'NR==FNR{ if(FNR>1) c[$1]++; next } FNR==1{ print; next } c[$1]==1' \
    "$OUT/external_all.txt" "$OUT/external_all.txt" > "$OUT/external_unique.txt"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$OUT/external_unique.txt" > "$OUT/external_score.txt"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$4}'    "$OUT/external_unique.txt" > "$OUT/external_pvalues.txt"
cp derived/prs_los/range_list.txt "$OUT/range_list.txt"   # same nine thresholds as the LOS PRS

$PLINK2 --bfile "$GENO" --keep "$SETS/gwas_keep.txt" --clump "$OUT/external_pvalues.txt" \
        --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --threads 8 --out "$OUT/external_clumped"
CL=$(ls "$OUT"/external_clumped.clump* | head -1)
awk 'NR>1 && $3!="" && $3!="NA" {print $3}' "$CL" | sort -u > "$OUT/clumped_snps.txt"

$PLINK2 --bfile "$GENO" --extract "$OUT/clumped_snps.txt" \
        --score "$OUT/external_score.txt" 1 2 3 header cols=+scoresums \
        --q-score-range "$OUT/range_list.txt" "$OUT/external_pvalues.txt" 1 2 header \
        --threads 8 --out "$OUT/external_prs"
