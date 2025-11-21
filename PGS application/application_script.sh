#!/bin/bash

while getopts ":j:m:" opt; do
  case $opt in
    j) JOBPATH="$OPTARG"
    ;;
    m) METHOD="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
printf "\033[33mArgument Jobpath is %s  \033[0m\n" "$JOBPATH"

# Set parameters
APP_PATH=/root/biosoft/get_picture/
PGS_APP=${APP_PATH}PGS_application.R
TEST_PATH=/disk/testSet/
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
PLINK=/root/biosoft/plink

# Perform PGS application
## Set PGS parameters
OUTPATH=${JOBPATH}/${METHOD}/
SUMMINFO=${JOBPATH}/${METHOD}/summary_info.txt
CATEGORY=`sed -n '1p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
if [ $CATEGORY = "cross_ancestry" ]; then
	ANCESTRY=`sed -n '9p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
else
	ANCESTRY=`sed -n '3p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
fi
IN_UKBB=`sed -n '4p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
SEX=`sed -n '5p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
if [ "$SEX" = "Mixed" ]; then
	SEX_LABEL="All"
else
	SEX_LABEL="$SEX"
fi
PHENOCODE=`sed -n '6p' ${SUMMINFO} | sed 's/^[^\t]\+[\t]\+//'`
EFF_FILE=${OUTPATH}esteff.txt
## Determine PGS
if [ "$IN_UKBB" = "Yes" ]; then
	echo 'PGSFusion does not support summary statistics including UKBB samples for PGS application!'
fi
if [ "$PHENOCODE" = "Unavailable" ]; then
	echo 'PGSFusion does not include the trait for PGS application!'
fi
if [ "$SEX_LABEL" = "Unknown" ]; then
	echo 'PGSFusion does not support Unknown for PGS application!'
fi
if [ "$IN_UKBB" = "No" ]; then
	## Calculate PGS
	mkdir ${OUTPATH}/pred_pheno
	PRED_PREFIX=${OUTPATH}/pred_pheno/pred_hm3_${SEX_LABEL}
	for chr in `seq 1 22`
	do
		TEST_BFILE=${TEST_PATH}${ANCESTRY}/genotype/${SEX_LABEL}/chr${chr}
		${PLINK} --bfile ${TEST_BFILE} --score ${EFF_FILE} 1 2 3 sum --out ${PRED_PREFIX}_chr${chr}
		rm -rf ${PRED_PREFIX}_chr${chr}.log
		rm -rf ${PRED_PREFIX}_chr${chr}.nopred
	done
	## Analyze PGS
	mkdir ${OUTPATH}/performance/
	mkdir ${OUTPATH}/jointeffect/
	${Rscript} ${PGS_APP} --phenocode ${PHENOCODE} --sex ${SEX_LABEL} --anc ${ANCESTRY} --out_path ${OUTPATH}
fi
