#!/usr/bin/env bash
# Applies the two frozen PRS (LOS weights at PT 0.01; external weights at PT 0.1) to the Guangzhou cohort.
# Matching by rsID (independent of genome build); palindromic variants and duplicated rsIDs excluded.
# Per-chromosome scoring with PLINK 1.9, then summed and averaged per allele.
set -eu
PLINK=${PLINK:-plink}
GEN=data/genotypes/guangzhou          # per-chromosome filesets: guangzhou_chr{1..22}.bed/.bim/.fam
LOSBIM=data/genotypes/los_qc.bim
LOS=derived/prs_los; EXT=derived/prs_external
OUT=derived/prs_guangzhou; mkdir -p "$OUT/dups"

awk 'NR==FNR{c[$1]=1; next} ($1 in c) && $2<0.01 {k[$1]=1} END{for(i in k) print i}' \
    "$LOS/clumped_snps.txt" "$LOS/snp_pvalues.txt" > "$OUT/los_ids.txt"
awk 'NR==FNR{k[$1]=1; next} ($2 in k){a1=toupper($5); a2=toupper($6);
     pal=((a1=="A"&&a2=="T")||(a1=="T"&&a2=="A")||(a1=="C"&&a2=="G")||(a1=="G"&&a2=="C"));
     if(!pal && $2 ~ /^rs/) print $2}' "$OUT/los_ids.txt" "$LOSBIM" > "$OUT/los_ids_ok.txt"
awk 'NR==FNR{k[$1]=1; next} ($1 in k) && !seen[$1]++ {print $1"\t"$2"\t"$3}' \
    "$OUT/los_ids_ok.txt" "$LOS/score_weights.txt" > "$OUT/score_los_weights.txt"
awk 'NR==FNR{c[$1]=1; next} ($1 in c) && $2<0.1 {k[$1]=1} END{for(i in k) print i}' \
    "$EXT/clumped_snps.txt" "$EXT/external_pvalues.txt" > "$OUT/external_ids.txt"
awk 'NR==FNR{k[$1]=1; next} ($1 in k){print $1"\t"$2"\t"$3}' \
    "$OUT/external_ids.txt" "$EXT/external_score.txt" > "$OUT/score_external_weights.txt"
for F in "$OUT/score_los_weights.txt" "$OUT/score_external_weights.txt"; do
  awk 'NR==FNR{c[$1]++; next} c[$1]==1' "$F" "$F" > "$F.tmp" && mv "$F.tmp" "$F"
done

for chr in {1..22}; do
  awk '{c[$2]++} END{for(i in c) if(c[i]>1) print i}' "$GEN/guangzhou_chr${chr}.bim" > "$OUT/dups/chr${chr}.txt"
done
for TAG in los_weights external_weights; do
  D="$OUT/perchr_${TAG}"; mkdir -p "$D"
  for chr in {1..22}; do
    EXCL=""; [ -s "$OUT/dups/chr${chr}.txt" ] && EXCL="--exclude $OUT/dups/chr${chr}.txt"
    $PLINK --bfile "$GEN/guangzhou_chr${chr}" --allow-extra-chr $EXCL \
           --score "$OUT/score_${TAG}.txt" 1 2 3 sum --out "$D/chr${chr}"
  done
  awk 'FNR==1{next} { cnt[$2]+=$4; cnt2[$2]+=$5; s[$2]+=$6 }
       END{ print "IID\tALLELE_CT\tDOSAGE_SUM\tSCORE_SUM\tSCORE_AVG";
            for(i in s) printf "%s\t%d\t%d\t%.6f\t%.8f\n", i, cnt[i], cnt2[i], s[i], s[i]/cnt[i] }' \
      "$D"/chr*.profile > "$OUT/guangzhou_prs_${TAG}.txt"
done
