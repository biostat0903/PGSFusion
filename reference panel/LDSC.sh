LDSC=/root/biosoft/ldsc/ldsc-master/ldsc.py
REF_IN=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/
REF_OUT=/disk/reference_pgsfusion/EUR_w_ld_chr/

source activate ldsc
for chr in `seq 1 22`
do
${LDSC} --bfile ${REF_IN}chr${chr}\
	--l2 --ld-wind-kb 1000\
	--out ${REF_OUT}${chr}
done
conda deactivate