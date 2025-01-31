LDAK=/root/biosoft/ldak6.1.linux

# Construct LD matrix
GENO=/disk/reference_pgsfusion/1kg/EAS/hm3_imp/mergeout
OUT=/disk/reference_pgsfusion/megaPRS/EAS/cor/
for chr in `seq 1 22`
do
${LDAK} --calc-cors ${OUT}chr${chr} \
        --bfile ${GENO} --window-kb 1000 --chr ${chr}
done

# Merge LD matrix
cd ${OUT}
rm -rf list.txt; for j in {1..22}; do echo "chr$j" >> list.txt; done
${LDAK} --join-cors cors --corslist list.txt
rm -rf chr*




CORTHIN=/gpfs/chencao/Temporary_Files/pgsfusion/thin
${LDAK} --thin ${CORTHIN} --bfile ${GENO} --window-prune .98 --window-kb 100

cd /gpfs/chencao/Temporary_Files/pgsfusion
awk < thin.in '{print $1, 1}' > weights.thin
TAG=/gpfs/chencao/Temporary_Files/pgsfusion/ldak.thin
${LDAK} --calc-tagging ${TAG} --bfile ${REF} --weights weights.thin --power -.25 --window-kb 1000 \
   --save-matrix YES


LDAK=/root/biosoft/ldak6.1.linux

GENO=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
cd 
${LDAK} --cut-weights sections --bfile ${GENO}
${LDAK} --calc-weights-all EUR_sections --bfile ${GENO}
mv EUR_sections/weights.short bld65



GENO=/disk/reference_pgsfusion/1kg/EAS/hm3_imp/mergeout
${LDAK} --calc-tagging bld.ldak --bfile ${GENO} --power -.25 --annotation-number 65 \
--annotation-prefix bld --save-matrix YES