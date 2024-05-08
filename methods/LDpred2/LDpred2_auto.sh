#!/bin/bash
while getopts ":s:" opt; do
  case $opt in
    s) Summary="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
printf "\033[33mArgument Summary is %s  \033[0m\n" "$Summary"

# Set parameters
Summary_stat=`sed -n '2p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
ancestry=`sed -n '3p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
ref_panel=/disk/reference_pgsfusion/LD_reference_3cM/$ancestry

# Run LDpred2-auto
LDpred2=/root/pgsfunsion/LDpred2_auto/LDpred2_auto.R
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript
outpath=`echo "$Summary" | awk -F'/summary' '{print $1}'`
${Rscript} ${LDpred2}  --summ ${Summary_stat} \
                       --reference ${ref_panel} \
                       --output ${outpath}
