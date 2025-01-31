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
# Summary_stat=/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/trait1.txt
ref_panel=/disk/reference_pgsfusion/SDPR/$ancestry

# Get parameters
nobs=`sed -n "2p" ${Summary_stat} | awk '{print $5}'`
nmis=`sed -n "2p" ${Summary_stat} | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
DIR=$(dirname ${Summary_stat})
awk '{print $2, $6, $7, $9, $11}' ${Summary_stat} > ${DIR}/GWAS_SDPR1.txt
awk -v var="$n" '{print $0, var}' ${DIR}/GWAS_SDPR1.txt > ${DIR}/GWAS_SDPR.txt
sed -i '1s/.*/SNP A1 A2 BETA P N/' ${DIR}/GWAS_SDPR.txt
sed -i 's/ /	/g' ${DIR}/GWAS_SDPR.txt

# Fit SDPR
SDPR=/root/biosoft/SDPR/SDPR-main/SDPR
outpath=`echo "$Summary" | awk -F'/summary_info.txt' '{print $1}'`
for chr in `seq 1 22`
do
${SDPR} -mcmc -ref_dir ${ref_panel} -ss ${DIR}/GWAS_SDPR.txt \
		-N ${n} -chr ${chr} -out ${outpath}/esteff_SDPR_chr_${chr}.txt \
		-n_threads 1
    # -iter 10 -burn 2
done

# Format output
for chr in `seq 1 22`
do
sed -i '1d' ${outpath}/esteff_SDPR_chr_${chr}.txt
done
cat ${outpath}/esteff_SDPR_chr_{1..22}.txt > ${outpath}/esteff.txt
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
PROCEFF=/root/pgsfusion/procEffect.R
${Rscript} ${PROCEFF} --method SDPR \
                      --esteff ${outpath}/esteff.txt \
					  --summ ${Summary_stat}

# Remove files
rm -rf ${DIR}/GWAS_SDPR.txt
rm -rf ${DIR}/GWAS_SDPR1.txt
rm -rf ${outpath}/esteff_SDPR_chr_*.txt
