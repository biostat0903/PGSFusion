# SDPRx: Construct LD matrix for two ancestries
PY=/root/anaconda3/envs/python3/bin/python

#######################
## Notice: Use Large Memory Size
#######################
REF1=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
REF2=/disk/reference_pgsfusion/1kg/EAS/hm3_imp/mergeout
OUT=/disk/reference_pgsfusion/SDPRx/EUR-EAS/cor/
for chr in `seq 1 22`
do
echo ${chr}
${PY} ${CALCLD} --ref_path1 ${REF1} --ref_path2 ${REF2}\
	  --chrom ${chr} --threads 1 --out ${OUT}
done

REF1=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
REF2=/disk/reference_pgsfusion/1kg/AFR/hm3_imp/mergeout
OUT=/disk/reference_pgsfusion/SDPRx/EUR-AFR/cor/
for chr in `seq 1 22`
do
echo ${chr}
${PY} ${CALCLD} --ref_path1 ${REF1} --ref_path2 ${REF2}\
	  --chrom ${chr} --threads 1 --out ${OUT}
done


# CALCLD=/gpfs/chencao/Temporary_Files/pgsfusion/calc_ref.py
# REF1=/disk/reference_pgsfusion/1kg/AFR/hm3_imp/mergeout
# REF2=/disk/reference_pgsfusion/1kg/EAS/hm3_imp/mergeout
# OUT=/disk/reference_pgsfusion/SDPRx/EAS-AFR/cor/
# for chr in `seq 1 22`
# do
# echo ${chr}
# ${PY} ${CALCLD} --ref_path1 ${REF1} --ref_path2 ${REF2}\
	  # --chrom ${chr} --threads 1 --out ${OUT}
# done
CALCLD=/gpfs/chencao/Temporary_Files/pgsfusion/calc_ref.py
REF1=/gpfs/chencao/Temporary_Files/pgsfusion/EAS/hm3_imp/mergeout
REF2=/gpfs/chencao/Temporary_Files/pgsfusion/AFR/hm3_imp/mergeout
OUT=/gpfs/chencao/Temporary_Files/pgsfusion/SDPRx/EAS-AFR/cor
for chr in 5
do
echo ${chr}
python ${CALCLD} --ref_path1 ${REF1} --ref_path2 ${REF2}\
	  --chrom ${chr} --threads 2 --out ${OUT}
done



# PopCorn: Esitmate cross-population score
PY=/root/anaconda3/envs/python3/bin/python
POPCORN=/root/biosoft/PopCorn/Popcorn-master/popcorn/__main__.py
## EUR-EAS
for chr in `seq 1 22`
do
REF1=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/chr${chr}
REF2=/disk/reference_pgsfusion/1kg/EAS/hm3_imp/chr${chr}
SCORE=/disk/reference_pgsfusion/SDPRx/EUR-EAS/score/chr${chr}.txt
${PY} ${POPCORN} compute --bfile1 ${REF1} --bfile2 ${REF2} ${SCORE}
rm -rf ${SCORE}.e.log
rm -rf ${SCORE}.o.log
done
## EUR-AFR
for chr in `seq 1 22`
do
REF1=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/chr${chr}
REF2=/disk/reference_pgsfusion/1kg/AFR/hm3_imp/chr${chr}
SCORE=/disk/reference_pgsfusion/SDPRx/EUR-AFR/score/chr${chr}.txt
${PY} ${POPCORN} compute --bfile1 ${REF1} --bfile2 ${REF2} ${SCORE}
rm -rf ${SCORE}.e.log
rm -rf ${SCORE}.o.log
done
## EAS-AFR
for chr in `seq 1 22`
do
REF1=/disk/reference_pgsfusion/1kg/EAS/hm3_imp/chr${chr}
REF2=/disk/reference_pgsfusion/1kg/AFR/hm3_imp/chr${chr}
SCORE=/disk/reference_pgsfusion/SDPRx/EAS-AFR/score/chr${chr}.txt
${PY} ${POPCORN} compute --bfile1 ${REF1} --bfile2 ${REF2} ${SCORE}
rm -rf ${SCORE}.e.log
rm -rf ${SCORE}.o.log
done
