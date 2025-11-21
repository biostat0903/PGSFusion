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
Summary_stat1=`sed -n '2p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
Summary_stat2=`sed -n '8p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
ancestry=`sed -n '3p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
if [ $ancestry = 'EUR' ]; then
	ref_panel=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
	# ref_panel=/disk/reference_pgsfusion/1kg/EUR/hm3_imp/merge_imp
else
	ref_panel=/disk/reference_pgsfusion/1kg/${ancestry}/hm3_imp/merge_imp
fi
r2=`sed -n '1p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
pval=`sed -n '2p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
nsin=`sed -n '3p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`
block=/disk/reference_pgsfusion/genome_block/${ancestry}/${ancestry}_LD_Block.txt

# Process summary
sed 's/ /	/g' ${Summary_stat1} > ${outpath}/sumstat1.txt
sed 's/ /	/g' ${Summary_stat2} > ${outpath}/sumstat2.txt

# Get N
nobs1=`sed -n "2p" ${Summary_stat1} | awk '{print $5}'`
nmis1=`sed -n "2p" ${Summary_stat1} | awk '{print $4}'`
n1=$(echo "${nobs1}+${nmis1}" | bc -l)
nobs2=`sed -n "2p" ${Summary_stat2} | awk '{print $5}'`
nmis2=`sed -n "2p" ${Summary_stat2} | awk '{print $4}'`
n2=$(echo "${nobs2}+${nmis2}" | bc -l)

# Estimate genetic correlation
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
cal_v=/root/pgsfusion/mtPGS/mtPGS_v.R
${Rscript} ${cal_v} --summ1 ${Summary_stat1} --summ2 ${Summary_stat2} \
					--nsin ${nsin} --output ${outpath} --anc ${ancestry}

# Fit mtPGS
mtPGS=/root/biosoft/mtPGS/mtPGS-main/src/mtPGS_int_only
plink=/root/biosoft/plink
PROCEFF=/root/pgsfusion/procEffect.R
${mtPGS} --summstat ${outpath}/sumstat1.txt ${outpath}/sumstat2.txt \
		 --n ${n1} ${n2} --block ${block} --target 0 --ref ${ref_panel} --mafMax 0.8 \
		 --vg ${outpath}/v_g.txt --ve ${outpath}/v_e.txt \
		 --r2 ${r2} --pval ${pval} --plink ${plink} --c_t ${outpath}/test \
		 --output ${outpath}/esteff_mtPGS_target ${outpath}/esteff_mtPGS_relevant
${Rscript} ${PROCEFF} --method mtPGS --esteff ${outpath}/esteff_mtPGS_target.txt --summ ${Summary_stat1}

# Remove files
rm -rf ${outpath}/sumstat1.txt
rm -rf ${outpath}/sumstat2.txt
rm -rf ${outpath}/v_g.txt
rm -rf ${outpath}/v_e.txt
rm -rf ${outpath}/test_clumping_result.log
rm -rf ${outpath}/test_clumping_result.nosex
rm -rf ${outpath}/test_LRT_output.txt
rm -rf ${outpath}/test_l_snp_list.txt
rm -rf ${outpath}/test_clumping_result.clumped
rm -rf ${outpath}/esteff_mtPGS_target.txt
rm -rf ${outpath}/esteff_mtPGS_relevant.txt
