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
reference=/disk/reference_pgsfusion/ukb${ancestry}_HM3/
outpath=`echo "$parameter" | awk -F'/parameter.txt' '{print $1}'`

# Summary_stat=/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/trait1.txt
# reference=/disk/reference_pgsfusion/ukbEUR_HM3/
# outpath=/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/SBAYESRC/

Summary_prefix=`basename $Summary_stat`
awk '(NR>1){snp=$2;a1=$6;a2=$7;fq=$8;beta=$9;bse=$10;pval=$11;N=($4+$5)}{print snp, a1, a2, fq, beta, bse, pval, N}(NR==1){print "SNP A1 A2 freq b se p N"}' ${Summary_stat} > ${outpath}/${Summary_prefix}
sed -i '1d' ${outpath}/${Summary_prefix}
sed -i 's/ /	/g' ${outpath}/${Summary_prefix}

# Fit SBayesRC
annot=/disk/reference_pgsfusion/annot_baseline2.2.txt
Rscript1=/root/anaconda3/envs/sbayesrc/bin/Rscript
${Rscript1} -e "SBayesRC::tidy(mafile='${outpath}/${Summary_prefix}', LDdir='$reference', \
               output='${outpath}/summary_tidy.ma', log2file=TRUE)"
${Rscript1} -e "SBayesRC::impute(mafile='${outpath}/summary_tidy.ma', LDdir='$reference', \
			    output='${outpath}/summary_imp.ma', log2file=TRUE)"
${Rscript1} -e "SBayesRC::sbayesrc(mafile='${outpath}/summary_imp.ma', LDdir='$reference', \
               outPrefix='${outpath}/summary_sbrc', annot='$annot', log2file=TRUE)"
			  #  niter=10, burn=5)"

# Format output
PROCEFF=/root/pgsfusion/procEffect.R
Rscript2=/root/anaconda3/envs/pgscalc2/bin/Rscript
${Rscript2} ${PROCEFF} --method SBayesRC \
					  --esteff ${outpath}/summary_sbrc.txt \
					  --summ ${Summary_stat}

# Remove file
 rm -rf ${outpath}/${Summary_prefix}.txt
 rm -rf ${outpath}/summary_tidy.*
 rm -rf ${outpath}/summary_imp.*
 rm -rf ${outpath}/summary_sbrc.*
