# Load packages
library(bigreadr)
library(dplyr)
library(plyr)

# Function 1: fit model with only covariables
est.coef <- function(pheno, cov){
  
  na_idx <- ifelse(is.na(pheno), T, F)
  pheno_na <- pheno[!na_idx]
  model_type <- ifelse(length(unique(pheno_na)) == 2, "binomial", "gaussian")
  if (model_type == "gaussian")
    coef <- lm(pheno_na ~ as.matrix(cov)[!na_idx, ])$coef
  else 
    coef <- glm(pheno_na ~ as.matrix(cov)[!na_idx, ], family = "binomial")$coef
  return(coef)
}

# Set parameters
PATH <- "/public/home/biostat04/Project/19_PGS_fusion/"
TRAIT <- "/public/home/Datasets/ukb/pheno/03_trait/participant.csv.gz"

# Load sample id and field id
sample_id <- list.files(paste0(PATH, "sample_id/")) %>% 
  alply(., 1, function(ff) fread2(paste0(PATH, "sample_id/", ff))[, 1])
names(sample_id) <- list.files(paste0(PATH, "sample_id/")) %>% 
  gsub(".fam", "", .)

# Get phenotype ID
pheno_idx <- list.files(paste0(PATH, "out_pheno/test_afr_all")) %>% grep("PFID", .)
pheno_files <- list.files(paste0(PATH, "out_pheno/test_afr_all"))[pheno_idx]

# Output covariables and coefficients
genetic_pc <- paste0("p22009_a", 1:20)
covar <- fread2(TRAIT, select = c("eid", "p21022", "p22001", 
                                   paste0("p22009_a", 1:20)))
covar$p22001 <- ifelse(covar$p22001 == "Female", 1, 0) ## Female = 1
covar_var <- c("Age", "Sex", paste0("PC", 1: 20))
covar_var_sex <- c("Age", paste0("PC", 1: 20))
cov_data <- alply(c(1: length(sample_id)), 1, function (ss){
  
  if (grepl("all", names(sample_id)[ss])){
    
    covar_s <- covar[match(sample_id[[ss]], covar$eid), -1]
    # write.table(covar_s, col.names = F, row.names = F, quote = F,
    #             file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov.txt"))
    if (grepl("val", names(sample_id)[ss])){
     
      covar_coef <- alply(c(1: length(pheno_files)), 1, function (pp){
        
        pheno <- fread2(paste0(PATH, "out_pheno/", names(sample_id)[ss], "/", pheno_files[pp]))[, 1]
        coef_s <- try(est.coef(pheno, covar_s), silent = T)
        if (!inherits(coef_s, "try-error")){
          
          write.table(coef_s, col.names = F, row.names = F, quote = F,
                      file = paste0(PATH, "out_coef/", names(sample_id)[ss], "/", pheno_files[pp]))
        }
        return(pp)
      })
    }
    # write.table(covar_var, col.names = F, row.names = F, quote = F,
    #             file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov_var.txt"))
  } else {
    
    covar_s <- covar[match(sample_id[[ss]], covar$eid), -c(1, 3)]
    # write.table(covar_s, col.names = F, row.names = F, quote = F,
    #             file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov.txt"))
    if (grepl("val", names(sample_id)[ss])){
      
      covar_coef <- alply(c(1: length(pheno_files)), 1, function (pp){
        
        pheno <- fread2(paste0(PATH, "out_pheno/", names(sample_id)[ss], "/", pheno_files[pp]))[, 1]
        coef_s <- try(est.coef(pheno, covar_s), silent = T)
        if (!inherits(coef_s, "try-error")){
          
          write.table(coef_s, col.names = F, row.names = F, quote = F,
                      file = paste0(PATH, "out_coef/", names(sample_id)[ss], "/", pheno_files[pp]))
        }
        return(pp)
      })
    }
    # write.table(covar_var_sex, col.names = F, row.names = F, quote = F,
    #             file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/cov_var.txt"))
  }
  return(ss)
}) 
