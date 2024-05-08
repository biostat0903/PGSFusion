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
window_size=`sed -n '2p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
r2=`sed -n '3p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
val_phenotype=`sed -n '5p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
plen=`sed -n '1p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
val_genotype=`sed -n '4p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
type=`sed -n '7p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
phenocode=`sed -n '6p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`

# Set CT method
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
CT=/root/pgsfunsion/CT/CT.R

# Run CT
if [ $ancestry = "EUR" ]
then
	ref_panel=/disk/reference_pgsfusion/1kg/EUR/hm3_imp/mergeout
	if [ "$sex" = "Mixed" ]
	then
		sex_label="All"
	else
		sex_label=$sex
	fi
	cov=/disk/validationSet/coef/${sex_label}/${phenocode}.txt
	${Rscript} ${CT}  --summ ${Summary_stat} --window ${window_size} \
					  --r2 ${r2} --val_phenotype ${val_phenotype} \
					  --val_genotype ${val_genotype} --plen ${plen} --dat ${dat_type}\
					  --output ${outpath} --ref_genotype ${ref_panel} --cov ${cov}
else
	echo 'PGSFusion does not support EAS/AFR for CT model!'
fi