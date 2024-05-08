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
ancestry1=`sed -n '3p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
ancestry2=`sed -n '9p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
PY=/root/anaconda3/envs/python3/bin/python
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`
val_genotype=`sed -n '1p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
ref=/disk/reference_pgsfusion/SDPRx/${ancestry1}_${ancestry2}/SDPRX
score=/disk/reference_pgsfusion/SDPRx/${ancestry1}_${ancestry2}_all_gen_eff.cscore

# Get sample size
nobs1=`sed -n "2p" ${Summary_stat1} | awk '{print $5}'`
nmis1=`sed -n "2p" ${Summary_stat1} | awk '{print $4}'`
n1=$(echo "${nobs1}+${nmis1}" | bc -l)
nobs2=`sed -n "2p" ${Summary_stat2} | awk '{print $5}'`
nmis2=`sed -n "2p" ${Summary_stat2} | awk '{print $4}'`
n2=$(echo "${nobs2}+${nmis2}" | bc -l)

# Process summary statistics for popcorn
awk '{print $2, $6, $7, $9, $10, $8}' ${Summary_stat1} > ${outpath}/summ1.txt
awk -v var="$n1" '{print $0, var}' ${outpath}/summ1.txt > ${outpath}/summ11.txt
sed -i '1s/.*/rsid A1 A2 beta SE af N/' ${outpath}/summ11.txt
sed -i 's/ /	/g' ${outpath}/summ11.txt
awk '{print $2, $6, $7, $9, $10, $8}' ${Summary_stat2} > ${outpath}/summ2.txt
awk -v var="$n2" '{print $0, var}' ${outpath}/summ2.txt > ${outpath}/summ21.txt
sed -i '1s/.*/rsid A1 A2 beta SE af N/' ${outpath}/summ21.txt
sed -i 's/ /	/g' ${outpath}/summ21.txt

# Estimate genetic correlation
popcorn=/root/biosoft/PopCorn/Popcorn-master/popcorn/__main__.py
${PY} ${popcorn} fit -v 1 --cfile ${score} \
						  --sfile1 ${outpath}/summ11.txt \
						  --sfile2 ${outpath}/summ21.txt \
						  ${outpath}/correlation.txt
sed -i 's/\t/_/g' ${outpath}/correlation.txt
fourth_line=$(sed -n '4p' ${outpath}/correlation.txt)
first_number=$(echo "$fourth_line" | grep -oE '[0-9]+(\.[0-9]+)?' | head -n 1)

# Process summary statistics for SDPRX
awk 'NR==1 {print $0} NR>1 {if($10!=0) print $0, $9 / $10; else print $0, "NaN"}' OFS=" " ${Summary_stat1} > ${outpath}/summ1.txt
awk -v var="$n1" '{print $0, var}' ${outpath}/summ1.txt > ${outpath}/summ11.txt
awk '{print $2, $6, $7, $12, $13}' ${outpath}/summ11.txt > ${outpath}/summ1.txt
sed -i '1s/.*/SNP A1 A2 Z N/' ${outpath}/summ1.txt
sed -i 's/ /	/g' ${outpath}/summ1.txt
awk 'NR==1 {print $0} NR>1 {if($10!=0) print $0, $9 / $10; else print $0, "NaN"}' OFS=" " ${Summary_stat2} > ${outpath}/summ2.txt
awk -v var="$n2" '{print $0, var}' ${outpath}/summ2.txt > ${outpath}/summ21.txt
awk '{print $2, $6, $7, $12, $13}' ${outpath}/summ21.txt > ${outpath}/summ2.txt
sed -i '1s/.*/SNP A1 A2 Z N/' ${outpath}/summ2.txt
sed -i 's/ /	/g' ${outpath}/summ2.txt

# Fit SDPRX
SDPRx=/root/biosoft/SDPRx/SDPRX-main/SDPRX.py
for chr in `seq 1 22`
do
	${PY} ${SDPRx} --ss1 ${outpath}/summ1.txt --ss2 ${outpath}/summ2.txt \
	               --N1 ${n1} --N2 ${n2} \
				   --force_shared TRUE --load_ld ${ref}/ \
				   --valid ${val_genotype}.bim --chr ${chr} \
				   --rho ${first_number} --out ${outpath}/esteff_SDPRx_chr${chr}
done

# Formate output
for chr in `seq 1 22`
do
sed -i '1d' ${outpath}/esteff_SDPRx_chr${chr}_1.txt
sed -i '1d' ${outpath}/esteff_SDPRx_chr${chr}_2.txt
done
cat ${outpath}/esteff_SDPRx_chr*_1.txt > ${outpath}/esteff_auxiliary.txt
cat ${outpath}/esteff_SDPRx_chr*_2.txt> ${outpath}/esteff.txt

# Remove files
sudo rm ${outpath}/esteff_SDPRx_chr*
sudo rm ${outpath}/summ1.txt
sudo rm ${outpath}/summ2.txt
sudo rm ${outpath}/summ11.txt
sudo rm ${outpath}/summ21.txt
sudo rm ${outpath}/correlation.txt.e.log
sudo rm ${outpath}/correlation.txt.o.log
sudo rm ${outpath}/correlation.txt



