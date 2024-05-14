#
library(bigreadr)
library(dplyr)

#
ukb_path <- "/public/home/Datasets/ukb/"
out_path <- "/public/home/biostat04/Project/19_PGS_fusion/"

eth <- "EUR"
n_pc <- 20
seed_use <- 20230816
n_each_sex <- 25000

# load sample info
eid_eth <- fread2(paste0(ukb_path, "pheno/eid_", eth, ".txt"),
                      header = F)[, 1, drop = T]
load(paste0(ukb_path, "pheno/sqc_",eth,".RData"))
sqc_eth <- sqc_eth[match(eid_eth, sqc_eth$eid), ]
all_female_idx <- which(sqc_eth$Sex == 0)
all_male_idx <- which(sqc_eth$Sex == 1)

## sample id
set.seed(seed_use + 1)
idx_all_female <- sample(all_female_idx, n_each_sex*2, replace = F) %>% sort()
set.seed(seed_use + 2)
idx_valid_famale <- sample(idx_all_female, n_each_sex, replace = F) %>% sort()
idx_test_famale <- setdiff(idx_all_female, idx_valid_famale) %>% sort()

set.seed(seed_use - 1)
idx_all_male <- sample(all_male_idx, n_each_sex*2, replace = F) %>% sort()
set.seed(seed_use - 2)
idx_valid_male <- sample(idx_all_male, n_each_sex, replace = F) %>% sort()
idx_test_male <- setdiff(idx_all_male, idx_valid_male) %>% sort()

##
eid_valid_all <- eid_eth[c(idx_valid_famale, idx_valid_male) %>% sort()]
eid_valid_female <- eid_eth[idx_valid_famale]
eid_valid_male <- eid_eth[idx_valid_male]
eid_test_all <- eid_eth[c(idx_test_famale, idx_test_male) %>% sort()]
eid_test_female <- eid_eth[idx_test_famale]
eid_test_male <- eid_eth[idx_test_male]

##### extract validation set 2 for coef calc #####
set.seed(seed_use + 6)
idx_valid2_famale <- sample(setdiff(all_female_idx, idx_all_female), 
                            n_each_sex, replace = F) %>% sort()
set.seed(seed_use - 6)
idx_valid2_male <- sample(setdiff(all_male_idx, idx_all_male), 
                          n_each_sex, replace = F) %>% sort()
##
eid_valid2_all <- eid_eth[c(idx_valid2_famale, idx_valid2_male) %>% sort()]
eid_valid2_female <- eid_eth[idx_valid2_famale]
eid_valid2_male <- eid_eth[idx_valid2_male]

##
write.table(data.frame(eid_valid_all, eid_valid_all),
            file = paste0(out_path, "02_ukb/valid/eid_valid_All.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_valid_female, eid_valid_female),
            file = paste0(out_path, "02_ukb/valid/eid_valid_Female.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_valid_male, eid_valid_male),
            file = paste0(out_path, "02_ukb/valid/eid_valid_Male.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_test_all, eid_test_all),
            file = paste0(out_path, "02_ukb/test/EUR/eid_test_All.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_test_female, eid_test_female),
            file = paste0(out_path, "02_ukb/test/EUR/eid_test_Female.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_test_male, eid_test_male),
            file = paste0(out_path, "02_ukb/test/EUR/eid_test_Male.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_valid2_all, eid_valid2_all),
            file = paste0(out_path, "02_ukb/valid2/eid_valid2_All.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_valid2_female, eid_valid2_female),
            file = paste0(out_path, "02_ukb/valid2/eid_valid2_Female.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)
write.table(data.frame(eid_valid2_male, eid_valid2_male),
            file = paste0(out_path, "02_ukb/valid2/eid_valid2_Male.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)

#################
library(bigreadr)
library(dplyr)

#
ukb_path <- "/public/home/Datasets/ukb/"
out_path <- "/public/home/biostat04/Project/19_PGS_fusion/"

for (eth in c("AFR", "ASA")) {
  # load sample info
  load(paste0(ukb_path, "pheno/sqc_",eth,".RData"))
  eid_test_all <- fread2(paste0(ukb_path, "pheno/eid_", eth, ".txt"),
                         header = F)[, 1, drop = T]
  sex_eth <- sqc_eth$Sex[match(eid_test_all, sqc_eth$eid)]
  eid_test_female <- eid_test_all[which(sex_eth == 0)]
  eid_test_male <- eid_test_all[which(sex_eth == 1)]

  write.table(data.frame(eid_test_all, eid_test_all),
              file = paste0(out_path, "02_ukb/test/",eth , "/eid_test_All.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  write.table(data.frame(eid_test_female, eid_test_female),
              file = paste0(out_path, "02_ukb/test/",eth , "/eid_test_Female.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  write.table(data.frame(eid_test_male, eid_test_male),
              file = paste0(out_path, "02_ukb/test/",eth , "/eid_test_Male.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)


}

###################
# load eid for valid and test
rm(list = ls())
gc()
library(bigreadr)
library(dplyr)
library(glue)

setwd("/public/home/biostat04/Project/19_PGS_fusion/")

valid_eid_all <- fread2("02_ukb/valid/eid_valid_All.txt")[, 1, drop = T]
valid_eid_female <- fread2("02_ukb/valid/eid_valid_Female.txt")[, 1, drop = T]
valid_eid_male <- fread2("02_ukb/valid/eid_valid_Male.txt")[, 1, drop = T]
test_eid_eur_all <- fread2("02_ukb/test/EUR/eid_test_All.txt")[, 1, drop = T]
test_eid_eur_female <- fread2("02_ukb/test/EUR/eid_test_Female.txt")[, 1, drop = T]
test_eid_eur_male <- fread2("02_ukb/test/EUR/eid_test_Male.txt")[, 1, drop = T]
test_eid_afr_all <- fread2("02_ukb/test/AFR/eid_test_All.txt")[, 1, drop = T]
test_eid_afr_female <- fread2("02_ukb/test/AFR/eid_test_Female.txt")[, 1, drop = T]
test_eid_afr_male <- fread2("02_ukb/test/AFR/eid_test_Male.txt")[, 1, drop = T]
test_eid_asa_all <- fread2("02_ukb/test/ASA/eid_test_All.txt")[, 1, drop = T]
test_eid_asa_female <- fread2("02_ukb/test/ASA/eid_test_Female.txt")[, 1, drop = T]
test_eid_asa_male <- fread2("02_ukb/test/ASA/eid_test_Male.txt")[, 1, drop = T]
valid2_eid_all <- fread2("02_ukb/valid2/eid_valid2_All.txt")[, 1, drop = T]
valid2_eid_female <- fread2("02_ukb/valid2/eid_valid2_Female.txt")[, 1, drop = T]
valid2_eid_male <- fread2("02_ukb/valid2/eid_valid2_Male.txt")[, 1, drop = T]

id_list <- list("valid_All" = valid_eid_all,
                "valid_Female" = valid_eid_female,
                "valid_Male" = valid_eid_male,
                "test_EUR_All" = test_eid_eur_all,
                "test_EUR_Female" = test_eid_eur_female,
                "test_EUR_Male" = test_eid_eur_male,
                "test_AFR_All" = test_eid_afr_all,
                "test_AFR_Female" = test_eid_afr_female,
                "test_AFR_Male" = test_eid_afr_male,
                "test_ASA_All" = test_eid_asa_all,
                "test_ASA_Female" = test_eid_asa_female,
                "test_ASA_Male" = test_eid_asa_male,
                "valid2_All" = valid2_eid_all,
                "valid2_Female" = valid2_eid_female,
                "valid2_Male" = valid2_eid_male)

saveRDS(id_list, file = "02_ukb/id_list.rds")

pheno_all <- fread2("/public/home/Datasets/ukb/ukb47503.csv.gz")
pheno_all <- pheno_all[match(unique(unlist(id_list)), pheno_all$eid),]
fwrite2(pheno_all, file = "02_ukb/pheno_all.txt", na = "NA", quote = T, sep = "\t")
system("gzip -f 02_ukb/pheno_all.txt")


#############
library(bigreadr)
library(dplyr)
library(stringr)

#
ukb_path <- "/public/home/Datasets/ukb/"
out_path <- "/public/home/biostat04/Project/19_PGS_fusion/"

id_list <- readRDS(paste0(out_path, "02_ukb/id_list.rds"))
# BMI_eid <- paste0("21001-", 0, ".0")
cov_var <- c("Age", "Sex", paste0("PC", 1:20))
sqc_all <- fread2("/public/home/Datasets/ukb/pheno/sqc.txt")

for (type in c("valid", "test_EUR", "test_AFR", "test_ASA", "valid2")) {
  
  for (sex in c("All", "Female", "Male")) {
    
    out_pathx <- paste0(out_path, "02_ukb/", 
                        gsub("\\_", "\\/", type), "/pheno/", sex, "/")
    eidx <- id_list[[paste0(type, "_", sex)]]
    if (sex == "All") {
      cov_var_use <- c("Age_of_attending_assessment_centre",
                       "Sex",
                       paste0("PC", 1:20))
    } else {
      cov_var_use <- c("Age_of_attending_assessment_centre",
                       paste0("PC", 1:20))
    }
    cov_df <- sqc_all[match(eidx, sqc_all$eid), 
                      cov_var_use]
    
    fwrite2(cov_df, file = paste0(out_pathx, "cov.txt"),
            sep = "\t", na = "NA", col.names = F)
    fwrite2(data.frame(cov_var), file = paste0(out_pathx, "cov_var.txt"),
            sep = "\t", col.names = F)
    
    print(paste0("MGS: Cov file for ", sex, " in ", type, " is OK!"))
  }
  
}


