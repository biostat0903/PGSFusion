
library(bigreadr)

# copy the AFR and EAS fam file to {PATH}

# Set parameters
PATH <- "/public/home/biostat04/Project/19_PGS_fusion/sample_id/"

# Process EAS test set
eas_fam <- fread2(paste0(PATH, "test_eas_all.fam"))
eas_fam_male <- eas_fam[eas_fam$V5 == 1, ]
eas_fam_female <- eas_fam[eas_fam$V5 == 2, ]
write.table(eas_fam_male, file = paste0(PATH, "test_eas_male.fam"), 
            col.names = F, row.names = F, quote = F)
write.table(eas_fam_female, file = paste0(PATH, "test_eas_female.fam"), 
            col.names = F, row.names = F, quote = F)

# Process AFR test set
afr_fam <- fread2(paste0(PATH, "test_afr_all.fam"))
afr_fam_male <- afr_fam[afr_fam$V5 == 1, ]
afr_fam_female <- afr_fam[afr_fam$V5 == 2, ]
write.table(afr_fam_male, file = paste0(PATH, "test_afr_male.fam"), 
            col.names = F, row.names = F, quote = F)
write.table(afr_fam_female, file = paste0(PATH, "test_afr_female.fam"), 
            col.names = F, row.names = F, quote = F)

# Process EUR validation and test set
eur_fam <- fread2("/public/home/Datasets/ukb/geno/EUR/hm3/chr1.fam")
eur_fam_male <- eur_fam[eur_fam$V5 == 1, ]
eur_fam_female <- eur_fam[eur_fam$V5 == 2, ]
set.seed(20250407)
eur_fam_male_sub <- sample(eur_fam_male$V1, 50000)
set.seed(20250407)
eur_fam_female_sub <- sample(eur_fam_female$V1, 50000)
## validation set
eur_fam_male_val <- eur_fam[eur_fam$V1 %in% eur_fam_male_sub[1: 25000], ]
eur_fam_female_val <- eur_fam[eur_fam$V1 %in% eur_fam_female_sub[1: 25000], ]
eur_fam_all_val <- rbind(eur_fam_male_val, eur_fam_female_val)
## test set
eur_fam_male_test <- eur_fam[eur_fam$V1 %in% eur_fam_male_sub[25001: 50000], ]
eur_fam_female_test <- eur_fam[eur_fam$V1 %in% eur_fam_female_sub[25001: 50000], ]
eur_fam_all_test <- rbind(eur_fam_male_test, eur_fam_female_test)
## output
write.table(eur_fam_male_val, file = paste0(PATH, "val_male.fam"), 
            col.names = F, row.names = F, quote = F)
write.table(eur_fam_female_val, file = paste0(PATH, "val_female.fam"), 
            col.names = F, row.names = F, quote = F)
write.table(eur_fam_all_val, file = paste0(PATH, "val_all.fam"), 
            col.names = F, row.names = F, quote = F)

write.table(eur_fam_male_test, file = paste0(PATH, "test_eur_male.fam"), 
            col.names = F, row.names = F, quote = F)
write.table(eur_fam_female_test, file = paste0(PATH, "test_eur_female.fam"), 
            col.names = F, row.names = F, quote = F)
write.table(eur_fam_all_test, file = paste0(PATH, "test_eur_all.fam"), 
            col.names = F, row.names = F, quote = F)

# Process EUR reference panel 
eur_fam_male_resid <- eur_fam_male$V1[!eur_fam_male$V1%in%eur_fam_male_sub]
eur_fam_female_resid <- eur_fam_female$V1[!eur_fam_female$V1%in%eur_fam_female_sub]
set.seed(1988)
eur_fam_male_ref_id <- sample(eur_fam_male_resid, 1000) 
eur_fam_male_ref <- eur_fam_male[eur_fam_male$V1 %in% eur_fam_male_ref_id, ]
set.seed(1988)
eur_fam_female_ref_id <- sample(eur_fam_female_resid, 1000)
eur_fam_female_ref <- eur_fam_female[eur_fam_female$V1 %in% eur_fam_female_ref_id, ]
eur_fam_ref <- rbind(eur_fam_male_ref, eur_fam_female_ref)
write.table(eur_fam_ref, file = paste0(PATH, "ref_all.fam"), 
            col.names = F, row.names = F, quote = F)


