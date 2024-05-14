###### extract for valid and valid2 set ######
eth=EUR
geno_path=/public/home/Datasets/ukb/geno/${eth}/hm3_nosex/
out_path=/public/home/biostat04/Project/19_PGS_fusion/02_ukb/
plink=/public/home/biostat04/biosoft/plink/plink
#
for type in valid valid2
do
for sex in All Female Male
do
#
out_pathx=${out_path}${type}/geno/${sex}/
merge_list=${out_pathx}merge_list.txt
touch ${merge_list}
#
for chr in `seq 1 22`
do
#
id_list=${out_path}${type}/eid_${type}_${sex}.txt
bfile=${geno_path}chr${chr}

${plink} --bfile ${bfile} \
--keep ${id_list} \
--indiv-sort f ${id_list} \
--make-bed \
--out ${out_pathx}chr${chr}

if [[ $chr -ne 1 ]]
then
echo ${out_pathx}chr${chr} >> ${merge_list}
fi

done
#
${plink} --bfile ${out_pathx}chr1 \
--merge-list ${merge_list} \
--make-bed \
--indiv-sort f ${id_list} \
--out ${out_pathx}merge

done
done

###### extract for test set ######
out_path=/public/home/biostat04/Project/19_PGS_fusion/02_ukb/test/
plink=/public/home/biostat04/biosoft/plink/plink

for eth in AFR ASA EUR
do

geno_path=/public/home/Datasets/ukb/geno/${eth}/hm3_nosex/

for sex in All Female Male
do

out_pathx=${out_path}${eth}/geno/${sex}/
merge_list=${out_pathx}merge_list.txt
touch ${merge_list}

# extract for each chr
for chr in `seq 1 22`
do

id_list=${out_path}${eth}/eid_test_${sex}.txt
bfile=${geno_path}chr${chr}

${plink} --bfile ${bfile} \
--keep ${id_list} \
--indiv-sort f ${id_list} \
--make-bed \
--out ${out_pathx}chr${chr}

# format merge_list
if [[ $chr -ne 1 ]]
then
echo ${out_pathx}chr${chr} >> ${merge_list}
fi

done
# memrge
${plink} --bfile ${out_pathx}chr1 \
--merge-list ${merge_list} \
--make-bed \
--indiv-sort f ${id_list} \
--out ${out_pathx}merge

done
done
