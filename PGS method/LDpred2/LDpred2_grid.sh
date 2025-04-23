#!/bin/bash
while getopts ":s:p:j:" opt; do
  case $opt in
    s) Summary="$OPTARG"
    ;;
    p) parameter="$OPTARG"
    ;;
    j) jobid="$OPTARG"
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
else
    let n_case=0
    let n_control=0
fi
val_phenotype=`sed -n '2p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
val_genotype=`sed -n '1p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
pheno_type=`sed -n '3p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`

# Set LDpred2 method
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
LDpred2=/root/pgsfusion/LDpred2/LDpred2_grid.R

# Run LDpred2
outpath=`echo "$parameter" | awk -F'/parameter' '{print $1}'`
if [ "$ancestry" = "EUR" ]
then
	ref_panel=/disk/reference_pgsfusion/EUR_UKB_ref/LD_matrix
	if [ "$sex" = "Mixed" ]
	then
    		sex_label="All"
	else
    		sex_label="$sex"
	fi
	cov=/disk/validationSet/phenotype/${sex_label}/cov.txt
	${Rscript} ${LDpred2}  --summ ${Summary_stat} --val_phenotype ${val_phenotype} \
						   --val_genotype ${val_genotype} --reference ${ref_panel} \
						   --dat ${dat_type}  --ncase ${n_case} \
						   --ncontrol ${n_control} --output ${outpath} \
						   --cov ${cov} --jobid ${jobid}
else
	echo 'PGSFusion does not support EAS/AFR for LDpred2 model!'
fi
