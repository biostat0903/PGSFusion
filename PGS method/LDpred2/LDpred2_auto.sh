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
dat_type=`sed -n '7p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
if [ $ancestry = 'EUR' ]; then
	ref_panel=/disk/reference_pgsfusion/EUR_UKB_ref/LD_matrix
else 
	ref_panel=/disk/reference_pgsfusion/1kg/${ancestry}/LD_matrix
fi
if [ ${dat_type} = "binary" ]
then
    n_case=`sed -n '8p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
    n_control=`sed -n '9p' ${Summary} | sed 's/^[^\t]\+[\t]\+//'`
else
    let n_case=0
    let n_control=0
fi

# Run LDpred2-auto
LDpred2=/root/pgsfusion/LDpred2/LDpred2_auto.R
Rscript=/root/anaconda3/envs/pgscalc2/bin/Rscript
outpath=`echo "$Summary" | awk -F'/summary' '{print $1}'`
${Rscript} ${LDpred2}  --summ ${Summary_stat} \
                       --reference ${ref_panel} \
                       --ncase ${n_case} \
                       --ncontrol ${n_control} \
                       --output ${outpath}
