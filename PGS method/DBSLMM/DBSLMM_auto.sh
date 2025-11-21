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
if [ $ancestry = 'EUR' ]; then
	# ref_panel=/disk/reference_pgsfusion/EUR_UKB_ref-archive/hm3_imp/
	ref_panel=/disk/reference_pgsfusion/1kg/EUR/hm3_imp/
else 
	ref_panel=/disk/reference_pgsfusion/1kg/${ancestry}/hm3_imp/
fi
BLOCK=/disk/reference_pgsfusion/genome_block/${ancestry}/
nobs=`sed -n "2p" ${Summary_stat} | awk '{print $5}'`
nmis=`sed -n "2p" ${Summary_stat} | awk '{print $4}'`
N=$(echo "${nobs}+${nmis}" | bc -l)
MODEL=DBSLMM

# DBSLMM-auto
DBSLMM=/root/biosoft/DBSLMM/DBSLMM.R
dbslmm=/root/biosoft/DBSLMM/dbslmm
Rscript=/root/anaconda3/envs/dbslmm/bin/Rscript
outpath=`echo "$parameter" | awk -F'/parameter' '{print $1}'`
${Rscript} ${DBSLMM} --summary ${Summary_stat} --dbslmm ${dbslmm} \
                     --type auto --model ${MODEL} \
					 --reference ${ref_panel} --block ${BLOCK} \
					 --N ${N} --outPath ${outpath}/

# Process results
#Summary_prefix=`echo "$Summary_stat" | awk -F'/' '{print $8}' | awk -F'.txt' '{print $1}'`
Summary_prefix=$(basename "$Summary_stat" .txt)
cat ${outpath}/${Summary_prefix}_chr{1..22}.dbslmm.txt > ${outpath}/esteff.txt
PROCEFF=/root/pgsfusion/procEffect.R
${Rscript} ${PROCEFF} --method DBSLMM --esteff ${outpath}/esteff.txt --summ ${Summary_stat}
rm -rf ${outpath}/${Summary_prefix}_chr*
rm -rf ${outpath}/merge_sub-*

