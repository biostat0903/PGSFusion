#!/bin/bash
while getopts ":s:p:" opt; do
  case $opt in
    s) Summary="$OPTARG"
    ;;
    p) parameter="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
printf "\033[33mArgument Summary is %s  \033[0m\n" "$Summary"
printf "\033[33mArgument parameter is %s  \033[0m\n" "$parameter"

# Set parameters
Summary_stat=`sed -n '2p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
ancestry=`sed -n '3p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
dat_type=`sed -n '7p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
phenocode=`sed -n '6p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
sex=`sed -n '5p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
if [ ${dat_type} = "binary" ]
then
    n_case=`sed -n '8p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
    n_control=`sed -n '9p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
    N=${n_case},${n_control}
else
    nobs=`sed -n "2p" ${Summary_stat} | awk '{print $5}'`
    nmis=`sed -n "2p" ${Summary_stat} | awk '{print $4}'`
    N=$(echo "${nobs}+${nmis}" | bc -l)
fi
val_phenotype=`sed -n '2p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
val_genotype=`sed -n '1p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`
BLOCK=/disk/reference_pgsfusion/genome_block/${ancestry}/

# Fit DBSLMM
DB_path=/root/biosoft/DBSLMM/
if [ $ancestry = "EUR" ]
then
	ref_panel=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/
	# ref_panel=/disk/reference_pgsfusion/1kg/EUR/hm3_imp/
	if [ "$sex" = "Mixed" ]
	then
    		sex_label="All"
	else
    		sex_label="$sex"
	fi
	cov=/disk/validationSet/phenotype/${sex_label}/cov.txt
	Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
	${Rscript} ${DB_path}/DBSLMM.R --summary ${Summary_stat} --dbslmm ${DB_path}/dbslmm \
	 							   --type tuning --model DBSLMM --block ${BLOCK} \
	 							   --reference ${ref_panel} --N ${N} --h2f 0.8,1,1.2 \
	 							   --outPath ${outpath}/
	Summary_prefix=$(basename "$Summary_stat" .txt)
	${Rscript} ${DB_path}/TUNE.R --summary ${Summary_stat} --dbslmm_eff "$outpath"/${Summary_prefix}\
								 --h2f 0.8,1,1.2 --cov ${cov}\
								 --validation_g ${val_genotype}/ --validation_p ${val_phenotype}
	cat "$outpath"/"$Summary_prefix"_chr{1..22}_best.dbslmm.txt > ${outpath}/esteff.txt
	# Process results
	PROCEFF=/root/pgsfusion/procEffect.R
	${Rscript} ${PROCEFF} --method DBSLMM --esteff ${outpath}/esteff.txt --summ ${Summary_stat}
	rm -rf ${outpath}/${Summary_prefix}_chr*
else
	echo 'PGSFusion does not support EAS/AFR for DBSLMM model!'
fi
