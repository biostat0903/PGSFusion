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

# Set parameters
Summary_stat=`sed -n '2p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
ancestry=`sed -n '3p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "${Summary}" | awk -F'/summary_info.txt' '{print $1}'`
# Summary_stat=/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/trait1.txt
# outpath=/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/
LDAK=/root/biosoft/ldak6.1.linux
ref_panel=/disk/reference_pgsfusion/megaPRS/${ancestry}

# Process summary
DIR=$(dirname ${Summary_stat})
awk '(NR>1){snp=$2;a1=$6;a2=$7;stat=($9/$10);n=($4+$5)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, stat, n}' ${Summary_stat}>${DIR}/quant.summaries 

# Fit MegaPRS
${LDAK} --mega-prs ${outpath}/megabayesr --model bayesr \
		--summary ${DIR}/quant.summaries \
		--cors ${ref_panel}/cor/cors --cv-proportion .1 \
		--power .25 --allow-ambiguous YES --check-sums NO

# Format output
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
PROCEFF=/root/pgsfusion/procEffect.R
${Rscript} ${PROCEFF} --method megaPRS \
					  --esteff ${outpath}/megabayesr.effects \
					  --summ ${Summary_stat}

# Remove files
rm -rf ${outpath}/megabayesr.best
rm -rf ${outpath}/megabayesr.cors
rm -rf ${outpath}/megabayesr.parameters
rm -rf ${outpath}/megabayesr.probs
rm -rf ${outpath}/megabayesr.progress
rm -rf ${outpath}/megabayesr.effects
