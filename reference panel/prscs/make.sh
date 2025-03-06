#!/bin/bash
#SBATCH --mem=250g
#SBATCH --cpus-per-task=4
#SBATCH --partition=cu,privority,batch01,fat

OUT=/gpfs/chencao/zhenghuili/data/geno/pgsfusion/ldblk_1kg/ldblk_1kg_eur

chr="$1"

i=1
while read -r line; do
echo "$i"
start=$(echo "$line" | awk '{print $2}')
start=$(( start + 1 ))
end=$(echo "$line" | awk '{print $3}')

plink --bfile ~/data/geno/pgsfusion/EUR_w_ld_chr/"$chr" --chr "$chr" --from-bp "$start" --to-bp "$end" --make-bed --out "$OUT"/"$chr"_"$i"

awk '{print $2}' "$OUT"/"$chr"_"$i".bim > "$OUT"/"$chr"_"$i".snplist

plink \
        --bfile ~/data/geno/pgsfusion/EUR_w_ld_chr/"$chr" \
        --r --ld-snp-list "$OUT"/"$chr"_"$i".snplist \
        --out "$OUT"/"$chr"_"$i" \
        --ld-window-r2 0 \
        --ld-window 100000 \
        --ld-window-kb 100000

    ((i+=1))
done < "$OUT"/region/chr"$chr".bed

rows=$(wc -l "$OUT"/region/chr"$chr".bed | cut -d' ' -f1)
echo "$rows"
python -u "$OUT"/script/make_matrix.py "$chr" "$rows"
