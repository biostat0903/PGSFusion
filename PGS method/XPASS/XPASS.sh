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
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`

if [ "$ancestry1" = 'EUR' ]; then
	ref1=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
	ref2=/disk/reference_pgsfusion/1kg/${ancestry2}/hm3_imp/merge_imp
	pc1=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/PC5.txt
	pc2=/disk/reference_pgsfusion/1kg/${ancestry2}/hm3_imp/PC5.txt
else
	ref1=/disk/reference_pgsfusion/1kg/${ancestry1}/hm3_imp/merge_imp
 	pc1=/disk/reference_pgsfusion/1kg/${ancestry1}/hm3_imp/PC5.txt
  	if [ "$ancestry2" = EUR ]; then
		ref2=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/merge_imp
    	pc2=/disk/reference_pgsfusion/EUR_UKB_ref/hm3_imp/PC5.txt
  	else
    	ref2=/disk/reference_pgsfusion/1kg/${ancestry2}/hm3_imp/merge_imp
    	pc2=/disk/reference_pgsfusion/1kg/${ancestry2}/hm3_imp/PC5.txt
  	fi
fi

#echo "$ref1"
#echo "$ref2"
#echo "$pc1"
#echo "$pc2"
#exit 1

# Get N
nobs1=`sed -n "2p" ${Summary_stat1} | awk '{print $5}'`
nmis1=`sed -n "2p" ${Summary_stat1} | awk '{print $4}'`
n1=$(echo "${nobs1}+${nmis1}" | bc -l)
nobs2=`sed -n "2p" ${Summary_stat2} | awk '{print $5}'`
nmis2=`sed -n "2p" ${Summary_stat2} | awk '{print $4}'`
n2=$(echo "${nobs2}+${nmis2}" | bc -l)

# Process summary statistics
awk 'NR==1 {print $0} NR>1 {if($10!=0) print $0, $9 / $10; else print $0, "NaN"}' OFS=" " ${Summary_stat1} > ${outpath}/summ1.txt
awk -v var="$n1" '{print $0, var}' ${outpath}/summ1.txt > ${outpath}/summ11.txt
awk '{print $2, $13, $12, $6, $7}' ${outpath}/summ11.txt > ${outpath}/summ1.txt
sed -i '1s/.*/SNP N Z A1 A2/' ${outpath}/summ1.txt
awk 'NR==1 {print $0} NR>1 {if($10!=0) print $0, $9 / $10; else print $0, "NaN"}' OFS=" " ${Summary_stat2} > ${outpath}/summ2.txt
awk -v var="$n2" '{print $0, var}' ${outpath}/summ2.txt > ${outpath}/summ21.txt
awk '{print $2, $13, $12, $6, $7}' ${outpath}/summ21.txt > ${outpath}/summ2.txt
sed -i '1s/.*/SNP N Z A1 A2/' ${outpath}/summ2.txt

# Fit XPASS
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
XPASS=/root/pgsfusion/XPASS/XPASS.R
${Rscript} ${XPASS} --summ1 ${outpath}/summ1.txt --summ2 ${outpath}/summ2.txt \
                    --ref1 ${ref1} --ref2 ${ref2} \
                    --pc1 ${pc1} --pc2 ${pc2} \
                    --output ${outpath}
