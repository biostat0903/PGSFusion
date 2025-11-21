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
if [ "$ancestry" = "EUR" ]
then
ref_panel=/disk/reference_pgsfusion/ldblk_ukbb_eur
val_genotype=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
elif [ "$ancestry" = "AFR" ]
then
ref_panel=/disk/reference_pgsfusion/ldblk_1kg/ldblk_1kg_afr
val_genotype=/disk/reference_pgsfusion/1kg/AFR/hm3_imp/mergeout
elif [ "$ancestry" = "EAS" ]
then
ref_panel=/disk/reference_pgsfusion/ldblk_1kg/ldblk_1kg_eas
val_genotype=/disk/reference_pgsfusion/1kg/EAS/hm3_imp/mergeout
fi

# Fit PRS-CS
CODEDIR=/root/pgsfunsion/PRScs/
PRSCS=/root/biosoft/PRScs/PRSCS_script.sh
val_genotype=`sed -n '1p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter' '{print $1}'`
sh ${PRSCS} -C ${CODEDIR} -s ${Summary_info} -L ${ref_panel} -T auto \
            -G ${val_genotype} -o ${outpath}
cat ${outpath}/esteff_auto_chr{1..22}.txt > ${outpath}/esteff.txt
rm -rf ${outpath}/esteff_auto_chr*

# Process results
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
PROCEFF=/root/pgsfusion/procEffect.R
${Rscript} ${PROCEFF} --method PRSCS --esteff ${outpath}/esteff.txt --summ ${Summary_stat}
