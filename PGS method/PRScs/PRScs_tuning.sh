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
Summary_info=`echo "$Summary_stat" | awk -F'.txt' '{print $1}'`
ancestry=`sed -n '3p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
sex=`sed -n '5p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`

# Set PRS-CS
CODEDIR=/root/biosoft/PRScs/
PRSCS=/root/biosoft/PRScs/PRSCS_script.sh
val_genotype=`sed -n '1p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
val_phenotype=`sed -n '2p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter' '{print $1}'`
ptype=`sed -n '7p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
if [ "$ancestry" = "EUR" ]; then
	
	ref_panel=/disk/reference_pgsfusion/ldblk_1kg/ldblk_1kg_eur
	if [ "$ptype" = "quantitative" ]; then
		index=r2
	else
		index=auc
	fi
	if [ "$sex" = "Mixed" ]; then
		sex_label="All"
	else
		sex_label=$sex
	fi
	# Fit PRS-CS
	sh ${PRSCS} -C ${CODEDIR} -s ${Summary_info} -L ${ref_panel} -T tuning \
				-G ${val_genotype} -o ${outpath} -P ${val_phenotype} -i ${index} -c ${ptype} -m ${sex_label}
	cat ${outpath}/esteff_bestphi_chr{1..22}.txt > ${outpath}/esteff.txt

	# Process results
	Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
	PROCEFF=/root/pgsfusion/procEffect.R
	${Rscript} ${PROCEFF} --method PRSCS --esteff ${outpath}/esteff.txt --summ ${Summary_stat}

	for chr in `seq 1 22`
	do
		rm -rf ${outpath}/esteff_bestphi_chr${chr}.txt
		rm -rf ${outpath}/val_bestphi_chr${chr}.profile
	done
else
	echo 'PGSFusion does not support EAS/AFR for PRSC-CS model!'
fi
