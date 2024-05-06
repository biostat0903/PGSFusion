#!/bin/bash
while getopts ":e:s:m:" opt; do
  case $opt in
    e) JOBPATH="$OPTARG"
    ;;
    s) SUMMINFO="$OPTARG"
    ;;
    m) METHOD="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
printf "\033[33mArgument Jobpath is %s  \033[0m\n" "$JOBPATH"
printf "\033[33mArgument SUMMINFO is %s  \033[0m\n" "$SUMMINFO"
printf "\033[33mArgument METHOD is %s  \033[0m\n" "$METHOD"

# for frontend
cp /root/biosoft/get_picture/fig.txt ${Jobpath}

# Set parameters
attr=`sed -n '1p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
ancestry=`sed -n '3p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
inc_ukbb=`sed -n '4p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
sex=`sed -n '5p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
if [ "$sex" = "Mixed" ]
then
    sex_label=All
else
    sex_label=$sex
fi
phenocode=`sed -n '6p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
plink=/root/biosoft/plink
PRS_app=/root/biosoft/get_picture/PRS_application.R
test_path=/disk/testSet/

# Set path
PATH=${Jobpath}${METHOD}
est_file=${PATH}esteff.txt

# Perform application 
if [ "$sex_label" = "Unknown" ]
then
	echo 'Do not obtain the sex information. The step can not be performed!'
fi
if [ "$inc_ukbb" = "No" ]
then
	if [ "$phenocode" = "Unavailable" ]
	then
		echo 'phenocode is not found! The step can not be performed!'
	else
		### Calculate PGS
		out_prefix=${PATH}pred_hm3_${sex_label}
		for chr in `seq 1 22`
		do
			bfile=${test_path}${ancestry}/genotype/${sex_label}/chr${chr}
			${plink} --bfile ${bfile} \
					 --score ${est_file} 1 2 3 sum \
					 --out ${out_prefix}_chr${chr}
			sudo rm ${out_prefix}_chr${chr}.log
			sudo rm ${out_prefix}_chr${chr}.nopred
		done
		### Visualization
		${Rscript} ${PRS_app} --phenocode ${phenocode} \
							  --sex ${sex_label} \
							  --anc ${ancestry} \
							  --out_path ${out_path}
	fi
else
	echo 'Summary statistics includes ukbb individuals. The step can not be performed!'
fi
