#!/bin/bash


# 1.
# for i in $(seq 2 22); do
#     sbatch --job-name=chr"$i" --out=chr"$i".out make.sh "$i"
# done

# 2.
for i in $(seq 2 22); do
{
    /gpfs/chencao/zhenghuili/software/anaconda3/envs/py27/bin/Rscript make_hdf5.R $i
} &
done

echo "ok"

wait