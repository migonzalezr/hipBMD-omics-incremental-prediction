#!/usr/bin/env bash
# Scores LOS participants with the two polygenic scores at two levels of restriction, to measure the cost
# of transferring each score to the Guangzhou cohort (Table S12):
#   transferable        palindromic variants and variants without an rsID removed (the score sent to Guangzhou)
#   guangzhou_variants  further restricted to variants present in the Guangzhou genotypes
# Duplicated IDs are excluded in every copy. Requires plink/05_prs_guangzhou.sh to have been run first.
set -eu
PLINK2=${PLINK2:-plink2}
GEN=data/genotypes/guangzhou          # guangzhou_chr{1..22}.bim
LOSBED=data/genotypes/los_qc
CH=derived/prs_guangzhou              # score_los_weights.txt, score_external_weights.txt
OUT=derived/prs_restricted; mkdir -p "$OUT"

cat "$GEN"/guangzhou_chr{1..22}.bim | awk '{print $2}' | sort -u > "$OUT/guangzhou_ids_all.txt"
cat "$GEN"/guangzhou_chr{1..22}.bim | awk '{c[$2]++} END{for(i in c) if(c[i]>1) print i}' | sort -u > "$OUT/guangzhou_ids_dup.txt"
comm -23 "$OUT/guangzhou_ids_all.txt" "$OUT/guangzhou_ids_dup.txt" > "$OUT/guangzhou_ids_ok.txt"

awk '{c[$2]++} END{for(i in c) if(c[i]>1) print i}' "${LOSBED}.bim" | sort -u > "$OUT/los_ids_dup.txt"
EXCL=""; [ -s "$OUT/los_ids_dup.txt" ] && EXCL="--exclude $OUT/los_ids_dup.txt"

for TAG in los_weights external_weights; do
  awk 'NR==FNR{k[$1]=1; next} ($1 in k)' "$OUT/guangzhou_ids_ok.txt" "$CH/score_${TAG}.txt" > "$OUT/score_${TAG}_guangzhou_variants.txt"
  D=$(awk '{c[$1]++} END{n=0; for(i in c) if(c[i]>1) n++; print n}' "$OUT/score_${TAG}_guangzhou_variants.txt")
  [ "$D" -eq 0 ] || { echo "duplicated IDs in restricted score $TAG"; exit 1; }
done

for LEVEL in transferable guangzhou_variants; do
  for TAG in los_weights external_weights; do
    if [ "$LEVEL" = "transferable" ]; then SCF="$CH/score_${TAG}.txt"; else SCF="$OUT/score_${TAG}_guangzhou_variants.txt"; fi
    $PLINK2 --bfile "$LOSBED" $EXCL --score "$SCF" 1 2 3 cols=+scoresums --out "$OUT/${TAG}_${LEVEL}"
  done
done
