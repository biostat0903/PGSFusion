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
AnnoPred=/root/biosoft/AnnoPred-master/AnnoPred.py
PY=/root/anaconda3/envs/python27/bin/python2.7
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
sel=/root/pgsfusion/AnnoPred/annopred_sel.R
Summary_stat=`sed -n '2p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
ancestry=`sed -n '3p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`

# Fit AnnoPred
if [ $ancestry = "EUR" ]
then
	ref_panel=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
	val_genotype=`sed -n '1p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
	annotation=`sed -n '2p' ${parameter}| sed 's/^[^\t]\+[\t]\+//'`
	## get N
	nobs=`sed -n "2p" ${Summary_stat} | awk '{print $5}'`
	nmis=`sed -n "2p" ${Summary_stat} | awk '{print $4}'`
	n=$(echo "${nobs}+${nmis}" | bc -l)
	## coordiante summary statistics 
	dir=$(dirname ${Summary_stat})
	awk '{print $1, $2, $6, $7, $3, $9, $11}' ${Summary_stat} > ${dir}/GWAS_sumstats1.txt
	sed -i 's/^\([^[:blank:]]\+\)/chr\1/' ${dir}/GWAS_sumstats1.txt
	awk 'NR==1 {print $0} NR>1 {printf "%s %s %s %s %d %.10f %.300f\n", $1, $2, $3, $4, $5, exp($6), $7}' ${dir}/GWAS_sumstats1.txt > ${dir}/GWAS_sumstats.txt
	sed -i '1s/.*/hg19chrc snpid a1 a2 bp or p/' ${dir}/GWAS_sumstats.txt
	## 
	mkdir ${outpath}/test_output
	mkdir ${outpath}/tmp_test
	output=${outpath}/test_output
	temp=${outpath}/tmp_test
	for p in 1.0 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 3e-05 1e-05
	do
		#/usr/bin/time -o /root/pgsfusion_test/job-AnnoPred_file/${p}.txt ${python}
		${PY} ${AnnoPred}\
		  --sumstats=${dir}/GWAS_sumstats.txt\
		  --ref_gt=${ref_panel}\
		  --val_gt=${val_genotype}\
		  --coord_out=${output}/coord_out\
		  --N_sample=${n}\
		  --annotation_flag=${annotation}\
		  --P=${p}\
		  --local_ld_prefix=${temp}/local_ld\
		  --out=${output}/test\
		  --temp_dir=${temp}
	done
	## Select best model
	if [ -f ${output}/select_p.txt ]; then
		rm -r ${output}/select_p.txt
	fi
	for p in 1.0 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 3e-05 1e-05
	do
	for model in h2_inf h2_non_inf pT_inf pT_non_inf
	do
	file=${output}/test_${model}_auc_${p}.txt
	if [ ! -f ${file} ]; then
		echo "NA" >> ${output}/select_p.txt
		continue
	fi
	if echo "$model" | grep -q "non_inf"
	then
		
		# # Extract the numeric part from the last column and append to output.txt
		# #awk 'END { match($0, /[0-9]+(\.[0-9]+)?/); if (RSTART) print substr($0, RSTART, RLENGTH) }' "$file" >> ${output}/select_p.txt
		echo "$(tail -n 2 ${file} | head -n 1 | cut -d' ' -f7 | xargs echo)" >> ${output}/select_p.txt
		# #tail -n 2 "$file" | head -n 1 | awk 'END { if(match($0, /[0-9]+(\.[0-9]+)?/)) print substr($0, RSTART, RLENGTH); else print "" }' "$file" >> ${output}/select_p.txt
	else
		echo "$(tail -n 1 ${file} | head -n 1 | cut -d' ' -f7 | xargs echo)" >> ${output}/select_p.txt
		# #awk 'END { if(match($0, /[0-9]+(\.[0-9]+)?/)) print substr($0, RSTART, RLENGTH); else print "" }' "$file" >> ${output}/select_p.txt
	fi
	done
	done
	# Output effect file
	
	${Rscript} ${sel} --output ${outpath}/
	awk '{print $3, $4, $7, $6, $2}' ${outpath}/esteff.txt > ${outpath}/esteff1.txt
	sed '1d' ${outpath}/esteff1.txt > ${outpath}/esteff.txt
	rm -rf ${outpath}/test_output
	rm -rf ${outpath}/tmp_test
else
	echo 'PGSFusion does not support EAS/AFR for AnnoPred model!'
fi
