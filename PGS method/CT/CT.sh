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
sex=`sed -n '5p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
phenocode=`sed -n '6p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
window_size=`sed -n '2p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
r2=`sed -n '3p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
val_phenotype=`sed -n '5p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
plen=`sed -n '1p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
val_genotype=`sed -n '4p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`

# Set CT method
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
CT=/root/pgsfusion/CT/CT.R

# Run CT
if [ $ancestry = "EUR" ]
then
	ref_panel=/disk/reference_pgsfusion/EUR_UKB_ref
	if [ "$sex" = "Mixed" ]
	then
		sex_label="All"
	else
		sex_label=$sex
	fi
	${Rscript} ${CT}  --summ ${Summary_stat} --phenocode ${phenocode} \
					  --dat ${dat_type} --sex_label ${sex_label} --reference ${ref_panel} \
					  --window ${window_size} --r2 ${r2} --plen ${plen}\
					  --output ${outpath} --jobid ${jobid}
else
	echo 'PGSFusion does not support EAS/AFR for CT model!'
fi
