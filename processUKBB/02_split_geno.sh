#
project_path=/public/home/biostat04/Project/19_PGS_fusion/
geno_path=/public/home/Datasets/ukb/geno/

## validation, validation2 and test sets in EUR ancestry
for chr in `seq 1 22`
do

for type in valid valid2 test
do

for sex in All Female Male
do

bfile=${geno_path}EUR/hm3_nosex/chr${chr}
idx=${project_path}02_ukb/${type}/eid_${type}_${sex}.txt
outx=${project_path}02_ukb/${type}/geno/${sex}/chr${chr}
# valid
plink --bfile ${bfile} \
--keep ${idx} \
--make-bed \
--out ${outx}
done
done
done

## test sets in ASA and AFR ancestries
for chr in `seq 1 22`
do

for ancestry in AFR ASA
do

for sex in All Female Male
do

type=test_${ancestry}
bfile=${geno_path}${ancestry}/hm3_nosex/chr${chr}
idx=${project_path}02_ukb/${type}/eid_${type}_${sex}.txt
outx=${project_path}02_ukb/${type}/geno/${sex}/chr${chr}
# valid
plink --bfile ${bfile} \
--keep ${idx} \
--make-bed \
--out ${outx}
done
done
done

