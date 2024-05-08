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
sex=`sed -n '5p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
val_phenotype=`sed -n '2p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
val_genotype=`sed -n '1p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
pheno_type=`sed -n '3p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
type=`sed -n '7p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
phenocode=`sed -n '6p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`

# Set LDpred2 method
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
LDpred2=/root/pgsfunsion/LDpred2/LDpred2.R

# Run LDpred2
outpath=`echo "$parameter" | awk -F'/parameter' '{print $1}'`
if [ "$ancestry" = "EUR" ]
then
	if [ "$sex" = "Mixed" ]
	then
    		sex_label="All"
	else
    		sex_label="$sex"
	fi
	cov=/disk/validationSet/coef/${sex_label}/${phenocode}.txt
	${Rscript} ${LDpred2}  --summ ${Summary_stat} --val_phenotype ${val_phenotype} \
						   --val_genotype ${val_genotype} --output ${outpath}\
						   --dat ${dat_type} --cov ${cov}
else
	echo 'PGSFusion does not support EAS/AFR for LDpred2 model!'
fi