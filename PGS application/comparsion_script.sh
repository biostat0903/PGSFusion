#!/bin/bash

while getopts ":j:" opt; do
  case $opt in
    j) JOBPATH="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
printf "\033[33mArgument Jobpath is %s  \033[0m\n" "$JOBPATH"

# Set parameters
APP_PATH=/root/biosoft/get_picture/
PGS_COMP=${APP_PATH}PGS_comparison.R
Rscript=/root/anaconda3/envs/pgscalc/bin/Rscript

# Compare PGS
METHOD_NUM=`cat ${JOBPATH}/method_use.txt | wc -l`
if [ $METHOD_NUM -ne 1 ]; then
	${Rscript} ${PGS_COMP} --job_path ${JOBPATH}/
fi

