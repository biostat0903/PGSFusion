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

# 1. Set software
LDAK=/root/biosoft/ldak6.2.beta
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
PROCEFF=/root/pgsfusion/procEffect.R

# 2. Set parameter
analysis=$(sed -n '2p' ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
DIR=$(dirname `sed -n '4p' ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r'`)
outpath=${DIR}/GIGAPRS
cor_list=${DIR}/gigaprs_cor.txt
summary_list=${DIR}/gigaprs_summarystat.txt

# 3. Creat correlation and summary statitics file
touch ${cor_list}
touch ${summary_list}

# 4. Analysis A: Single trait/Meta data of multiple trait & Single ancestry
if [ "$analysis" = "A" ];then
	## Write ancestry and summary statistics
	Summary_stat=$(sed -n '4p' ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
	ancestry=$(sed -n '5p' ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
	outfile="${DIR}/giga_summ.txt"
	awk '(NR>1){snp=$2;a1=$6;a2=$7;sd=$10;stat=($9/$10);n=($4+$5)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && sd!=0 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, stat, n}' ${Summary_stat} > "${outfile}"
	echo "/disk/reference_pgsfusion/gigaPRS/HAPMAP.${ancestry}" > ${cor_list}
	echo "${DIR}/giga_summ.txt	F	1" > ${summary_list}
	## Fit model A
	cd ${outpath}
	echo "DEBUG: about to run LDAK" >&2
	${LDAK} --giga-prs gigaA --corslist ${cor_list} --sumslist ${summary_list}  --check-MCMC NO 
	## Format output
	${Rscript} ${PROCEFF} --method gigaPRS_A --esteff ${outpath}/gigaA.effects\
						  --summ ${Summary_stat}
	## Remove files
	rm -rf ${outpath}/gigaA*
fi

# 5. Analysis B: Multiple trait & Single ancestry
if [ "$analysis" = "B" ];then
	## Write ancestry and summary statistics
	for line in $(seq 4 6 $(wc -l < ${Summary})); do
		Summary_stat=$(sed -n "${line}p" ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
		ancestry=$(sed -n "$((line+1))p" ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
		iter=$(( (line - 4) / 6 + 1 ))
		if [ $iter -eq 1 ]; then
			first_ancestry="${ancestry}"
		elif [ "${ancestry}" != "${first_ancestry}" ]; then
			echo "Error: Ancestry mismatch at line ${line}: expected '${first_ancestry}', got '${ancestry}'" >&2
			exit 1
		fi
		outfile="${DIR}/giga_summ_${iter}.txt"
		awk '(NR>1){snp=$2;a1=$6;a2=$7;sd=$10;stat=($9/$10);n=($4+$5)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && sd!=0 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, stat, n}' ${Summary_stat} > "${outfile}"
		echo -e "/disk/reference_pgsfusion/gigaPRS/HAPMAP.${ancestry}" > ${cor_list}
		if [ $iter -eq 1 ]; then
			echo -e "${outfile}\tF\t1" >> ${summary_list}
		else
			echo -e "${outfile}\tS\t1" >> ${summary_list}
		fi
	done
	## Fit model B
	cd ${outpath}
	${LDAK} --giga-prs gigaB --corslist ${cor_list} --sumslist ${summary_list} --check-MCMC NO
	## Format output
	${Rscript} ${PROCEFF} --method gigaPRS_B --esteff ${outpath}/gigaB.effects\
						  --summ ${Summary_stat}
	## Remove files
	rm -rf ${outpath}/gigaB*
fi

# 6. Analysis C: Single trait & Multiple ancestry
if [ "$analysis" = "C" ];then
	## Write ancestry and summary statistics
	for line in $(seq 4 6 $(wc -l < ${Summary})); do
		Summary_stat=$(sed -n "${line}p" ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
		ancestry=$(sed -n "$((line+1))p" ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
		iter=$(( (line - 4) / 6 + 1 ))
		if [ $iter -eq 1 ]; then
			first_ancestry="${ancestry}"
		elif [ "${ancestry}" == "${first_ancestry}" ]; then
			echo "Error: Ancestry mismatch at line ${line}: expected '${first_ancestry}', got '${ancestry}'" >&2
			exit 1
		fi
		outfile="${DIR}/giga_summ_${iter}.txt"
		awk '(NR>1){snp=$2;a1=$6;a2=$7;sd=$10;stat=($9/$10);n=($4+$5)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && sd!=0 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, stat, n}' ${Summary_stat} > "${outfile}"
		echo -e "${outfile}\tF\t${iter}" >> ${summary_list}
		echo -e "/disk/reference_pgsfusion/gigaPRS/HAPMAP.${ancestry}" >> ${cor_list}
	done
	## Fit model C
	cd ${outpath}
	${LDAK} --giga-prs gigaC --corslist ${cor_list} --sumslist ${summary_list} --check-MCMC NO 
	## Format output
	${Rscript} ${PROCEFF} --method gigaPRS_C --esteff ${outpath}/gigaC.combined.effects\
						  --summ ${Summary_stat}
	## Remove files
	rm -rf ${outpath}/gigaC*
fi

# 7. Analysis D: Multiple trait & Multiple ancestry
if [ "$analysis" = "D" ];then
	## Write ancestry and summary statistics
	let s_count=0
	declare -A ancestry_index
	for line in $(seq 4 6 $(wc -l < ${Summary})); do
		Summary_stat=$(sed -n "${line}p" ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
		ancestry=$(sed -n "$((line+1))p" ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
		iter=$(( (line - 4) / 6 + 1 ))
		outfile="${DIR}/giga_summ_${iter}.txt"
		awk '(NR>1){snp=$2;a1=$6;a2=$7;sd=$10;stat=($9/$10);n=($4+$5)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && sd!=0 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, stat, n}' ${Summary_stat} > "${outfile}"
		if [ -z "${ancestry_index[$ancestry]+x}" ]; then
			s_count=$((s_count + 1))
			ancestry_index[$ancestry]=$s_count
			echo -e "/disk/reference_pgsfusion/gigaPRS/HAPMAP.${ancestry}" >> ${cor_list}
			echo -e "${outfile}\tF\t${s_count}" >> ${summary_list}
		else
			existing_idx=${ancestry_index[$ancestry]}
			echo -e "${outfile}\tS\t${existing_idx}" >> ${summary_list}
		fi
	done
	
	## Fit model D
	cd ${outpath}
	${LDAK} --giga-prs gigaPRS_D --corslist ${cor_list} --sumslist ${summary_list} --check-MCMC NO
		## Format output
	${Rscript} ${PROCEFF} --method gigaPRS_D --esteff ${outpath}/gigaPRS_D.combined.effects\
						  --summ ${Summary_stat}
	## Remove files
	rm -rf ${outpath}/gigaPRS_D*
fi

# 8. Analysis E: Single trait/Meta data of multiple trait & Single ancestry
if [ "$analysis" = "E" ];then
	## Write ancestry and summary statistics
	Summary_stat=$(sed -n '4p' ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
	ancestry_raw=$(sed -n '5p' ${Summary} | sed 's/^[^\t]\+[\t]\+//' | tr -d '\n\r')
	outfile="${DIR}/giga_summ.txt"
	awk '(NR>1){snp=$2;a1=$6;a1Freq=$8;a2=$7;sd=$10;stat=($9/$10);n=($4+$5)}(NR==1){print "Predictor A1 A1Freq A2 Z n"}(NR>1 && sd!=0 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a1Freq, a2, stat, n}' ${Summary_stat} > "${outfile}"
	echo "${DIR}/giga_summ.txt	M	1" > ${summary_list}
	IFS=',' read -ra ancestry_list <<< "${ancestry_raw}"
	for anc in "${ancestry_list[@]}"; do
		anc=$(echo "$anc" | tr -d ' ')
		echo "/disk/reference_pgsfusion/gigaPRS/HAPMAP.${anc}" >> ${cor_list}
	done
	## Fit model E
	cd ${outpath}
	echo "DEBUG: about to run LDAK" >&2
	${LDAK} --giga-prs gigaE --corslist ${cor_list} --sumslist ${summary_list} --check-MCMC NO 
	## Format output
	${Rscript} ${PROCEFF} --method gigaPRS_E --esteff ${outpath}/gigaE.effects\
						  --summ ${Summary_stat}
	## Remove files
	rm -rf ${outpath}/gigaE*
fi

# 9. Remove temp files
rm -rf ${cor_list}
rm -rf ${summary_list}