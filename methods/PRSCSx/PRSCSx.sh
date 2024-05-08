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
val_genotype=`sed -n '1p' ${parameter} | sed 's/^[^\t]\+[\t]\+//'`
ref_panel=/disk/reference_pgsfusion/ldblk_1kg

# Process summary statistics
dir=$(dirname ${Summary_stat1})
awk '{print $2, $6, $7, $9, $11}' ${Summary_stat1} > ${dir}/summ1.txt
sed -i '1s/.*/SNP A1 A2 BETA P/' ${dir}/summ1.txt
sed -i 's/ /	/g' ${dir}/summ1.txt
awk '{print $2, $6, $7, $9, $11}' ${Summary_stat2} > ${dir}/summ2.txt
sed -i '1s/.*/SNP A1 A2 BETA P/' ${dir}/summ2.txt
sed -i 's/ /	/g' ${dir}/summ2.txt

# Get n1,n2
nobs1=`sed -n "2p" ${Summary_stat1} | awk '{print $5}'`
nmis1=`sed -n "2p" ${Summary_stat1} | awk '{print $4}'`
n1=$(echo "${nobs1}+${nmis1}" | bc -l)
nobs2=`sed -n "2p" ${Summary_stat2} | awk '{print $5}'`
nmis2=`sed -n "2p" ${Summary_stat2} | awk '{print $4}'`
n2=$(echo "${nobs2}+${nmis2}" | bc -l)

# Run PRSCSx
PY=/root/anaconda3/envs/python3/bin/python
prscsx=/root/biosoft/PRScsx-master/PRScsx.py
outpath=`echo "$Summary" | awk -F'/summary_info.txt' '{print $1}'`
for chr in `seq 1 22`
do
${PY} ${prscsx} --ref_dir=${ref_panel} --bim_prefix=${val_genotype} \
                 --sst_file=${dir}/summ1.txt,${dir}/summ2.txt \
                 --n_gwas=${n1},${n2} --pop=${ancestry1},${ancestry2} --chrom=${chr} \
                 --out_dir=${outpath} --out_name=test --n_iter=500 --n_burnin=100
done

# Format summary statistics
for chr in `seq 1 22`
do
	mv ${outpath}/test_${ancestry1}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt \
	${outpath}/esteff_prscsx_${ancestry1}_chr${chr}.txt
	mv ${outpath}/test_${ancestry2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt \
	${outpath}/esteff_prscsx_${ancestry2}_chr${chr}.txt
done
cat ${outpath}/esteff_prscsx_${ancestry1}_chr*.txt  > ${outpath}/esteff11.txt
cat ${outpath}/esteff_prscsx_${ancestry2}_chr*.txt > ${outpath}/esteff12.txt
awk '{print $2, $4, $6}' ${outpath}/esteff11.txt > ${outpath}/esteff_auxiliary.txt
awk '{print $2, $4, $6}' ${outpath}/esteff12.txt > ${outpath}/esteff.txt

# Remove files
sudo rm ${outpath}/summ1.txt
sudo rm ${outpath}/summ2.txt
sudo rm ${outpath}/esteff_prscsx_${ancestry1}_chr*.txt
sudo rm ${outpath}/esteff11.txt
sudo rm ${outpath}/esteff_prscsx_${ancestry2}_chr*.txt
sudo rm ${outpath}/esteff12.txt
