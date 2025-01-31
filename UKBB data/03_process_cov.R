# Load packages
library(bigreadr)
library(dplyr)
library(plyr)

# Set parameters
PATH <- "/public/home/biostat04/Project/19_PGS_fusion/"
TRAIT <- "/public/home/Datasets/ukb/pheno/03_trait/participant.csv.gz"

# Load sample id and field id
sample_id <- list.files(paste0(PATH, "sample_id/")) %>% 
  alply(., 1, function(ff) fread2(ff)[, 1])
names(sample_id) <- list.files(paste0(PATH, "sample_id/")) %>% 
  gsub(".fam", "", .)

# Load covariable
genetic_pc <- paste0("p22009_a", 1:20)
covar <- fread2(TRAIT, select = c("eid", "p21022", "p22001", 
                                   paste0("p22009_a", 1:20)))
covar$p22001 <- ifelse(covar$p22001 == "Female", 1, 0) ## Female = 1
covar_var <- c("Age", "Sex", paste0("PC", 1: 20))
covar_var_sex <- c("Age", paste0("PC", 1: 20))

cov_data <- alply(c(1: length(sample_id)), 1, function (ss){
  
  if (grepl("all", names(sample_id)[ss])){
    
    covar_s <- covar[match(sample_id[[ss]], covar$eid), -1]
    write.table(covar_s, col.names = F, row.names = F, quote = F,
                file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov.txt"))
    write.table(covar_var, col.names = F, row.names = F, quote = F,
                file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov_var.txt"))
  } else {
    
    covar_s <- covar[match(sample_id[[ss]], covar$eid), -c(1, 3)]
    write.table(covar_s, col.names = F, row.names = F, quote = F,
                file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov.txt"))
    write.table(covar_var_sex, col.names = F, row.names = F, quote = F,
                file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov_var.txt"))
    
  }
  return(ss)
}) 

