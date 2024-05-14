#
library(bigreadr)
library(dplyr)
library(stringr)

project_path <- "/public/home/biostat04/Project/19_PGS_fusion/"
setwd(project_path)

## format code use
cancer_diag_code_b <- c("2453-")
diab_diag_code_b <- c("2443-")
ICD10_code_b <- c("41202-", "41204-", "41270-")
self_cancer_code_b <- c("20001-")
self_noncancer_code_b <- c("20002-")

## format code list for definition
code_df_b <- fread2("02_ukb/pheno_code_b.txt")
format_list_b <- lapply(1:nrow(code_df_b), function(x){
  #
  value_icd10x <- value_cancer <- value_noncancer <- 
    value_other_T <- value_other_F <- NA
  code_bx <- strsplit(code_df_b$UKBB_code[x], ", ") %>% 
    unlist() %>% 
    as.integer() %>% 
    paste0(., "-")
  
  # code for ICD10
  if (any(ICD10_code_b %in% code_bx)) {
    value_icd10x <- code_df_b$ICD10[x] %>% 
      as.character() %>%
      strsplit(. ,", ") %>% 
      unlist()
  }
  # code for cancer OR diabetes
  cancer <- any(cancer_diag_code_b %in% code_bx) 
  diabetes <- any(diab_diag_code_b %in% code_bx)

    # code for self-cancer
  if (any(self_cancer_code_b %in% code_bx)) {
    value_cancer <- code_df_b$Self_reported_code[x] %>% 
      as.character() %>%
      strsplit(. ,", ") %>% 
      unlist()
  }
  # code for self-noncancer
  if (any(self_noncancer_code_b %in% code_bx)) {
    value_noncancer <- code_df_b$Self_reported_code[x] %>% 
      as.character() %>%
      strsplit(. ,", ") %>% 
      unlist()
  }
  # code for other (more specific)
  other_code_bx <- setdiff(code_bx, 
                           c(ICD10_code_b, 
                             self_cancer_code_b, self_noncancer_code_b,
                             cancer_diag_code_b, diab_diag_code_b))
  if (length(other_code_bx) > 0) {
    value_other_T <- code_df_b$Yes[x] %>% 
      as.character() %>%
      strsplit(. ,", ") %>% 
      unlist()
    value_other_F <- code_df_b$No[x] %>% 
      as.character() %>%
      strsplit(. ,", ") %>% 
      unlist()
  }
  
  return(list("value_ICD10" = value_icd10x,
              "value_cancer" = value_cancer,
              "value_noncancer" = value_noncancer,
              "cancer" = cancer,
              "diabetes" = diabetes,
              "other_code_b" = other_code_bx,
              "value_other_T" = value_other_T,
              "value_other_F" = value_other_F))
})
#
names(format_list_b) <- code_df_b$PFID
all_other_code_b <- lapply(format_list_b, function(x){
  x[["other_code_b"]]
}) %>% unlist %>% unique()

code_b_list_all <- list("cancer_diag_code_b" = cancer_diag_code_b,
                        "diab_diag_code_b" = diab_diag_code_b,
                        "ICD10_code_b" = ICD10_code_b,
                        "self_cancer_code_b" = self_cancer_code_b,
                        "self_noncancer_code_b" = self_noncancer_code_b,
                        "all_other_code_b" = all_other_code_b)
## set corresponding variables
ukb_header <- fread2("/public/home/Datasets/ukb/ukb47503_header.txt", header = F) %>% 
  t() %>% as.vector()
select_var_list_b <- lapply(code_b_list_all, function(code_listx){
  lapply(code_listx, function(codex){
    ukb_header[grep(paste0("^", codex), ukb_header)]
  }) %>% unlist
})
#
names(select_var_list_b) <- c("select_var_cancer_diag",
                              "select_var_diab_diag",
                              "select_var_icd10",
                              "select_var_cancer",
                              "select_var_noncancer",
                              "select_var_other")

## load data with selected variables
id_list <- readRDS("02_ukb/id_list.rds")
pheno_sel_b <- fread2("02_ukb/pheno_all.txt.gz", 
                      select = c("eid", unlist(select_var_list_b) %>%
                                   unique() %>% as.vector()))
pheno_sel_b[pheno_sel_b < 0] <- NA # set all positive values as NA
# set cancer diagnostic status
diag_cancer <- rep(0, nrow(pheno_sel_b))
diag_cancer[rowSums(pheno_sel_b[, select_var_list_b$select_var_cancer_diag] == 1,
                    na.rm = T) > 0] <- 1 
diag_cancer[rowMeans(is.na(pheno_sel_b[, select_var_list_b$select_var_cancer_diag])) == 1] <- NA

# set cancer diabetes status
diag_diab <- rep(0, nrow(pheno_sel_b))
diag_diab[rowSums(pheno_sel_b[, select_var_list_b$select_var_diab_diag] == 1,
                    na.rm = T) > 0] <- 1 
diag_diab[rowMeans(is.na(pheno_sel_b[, select_var_list_b$select_var_diab_diag])) == 1] <- NA

## define base phenotype 
# (Set Positive for individuals with at least ONE Positive in multiple records)
source("code/pheno_def_fun.R")
pheno_df <- lapply(names(format_list_b), function(PFIDx){
  
  format_list_bx <- format_list_b[[PFIDx]]
  ##
  diag_ICD10x <-   deter.b(code_list = format_list_bx$value_ICD10,
                           equal = F,
                           target_df = pheno_sel_b[, select_var_list_b[["select_var_icd10"]]])
  
  diag_cancerx <-   deter.b(code_list = format_list_bx$value_cancer,
                            equal = T,
                            target_df = pheno_sel_b[, select_var_list_b[["select_var_cancer"]]])
  
  diag_noncancerx <-   deter.b(code_list = format_list_bx$value_noncancer,
                               equal = T,
                               target_df = pheno_sel_b[, select_var_list_b[["select_var_noncancer"]]])
    
  diag_otherx <- rep(NA, nrow(pheno_sel_b))
  if (length(format_list_bx$other_code_b) != 0) {
    
    other_var_bx <- ukb_header[grep(paste0("^", format_list_bx$other_code_b), 
                                    ukb_header)]
    
    F_idx <- lapply(other_var_bx, function(other_var_bxx){
      which(as.character(pheno_sel_b[, other_var_bxx, drop = T]) %in% 
              as.character(format_list_bx$value_other_F))
    }) %>% unlist() %>% unique()
    
    T_idx <- lapply(other_var_bx, function(other_var_bxx){
      which(as.character(pheno_sel_b[, other_var_bxx, drop = T]) %in% 
              as.character(format_list_bx$value_other_T))
    }) %>% unlist() %>% unique()
    
    diag_otherx[F_idx] <- 0
    diag_otherx[T_idx] <- 1
    
    
  } 
  
  ## determine the phenotype of PFIDx
  diag_matx <- cbind(diag_ICD10x,
                     diag_cancerx,
                     diag_noncancerx,
                     diag_otherx) 
  
  diagx <- rep(0, nrow(pheno_sel_b))
  diagx[rowSums(diag_matx == T, na.rm = T) > 0] <- 1 
  diagx[rowMeans(is.na(diag_matx)) == 1] <- NA
  
  if (format_list_bx$cancer) {
    diagx <- diagx & diag_cancer
  }

  if (format_list_bx$diabetes) {
    diagx <- diagx & diag_diab
  }
  
  # saveRDS(diagx, file = paste0("diag/", PFIDx, ".rds"))
  print(PFIDx)
  return(diagx)
  
}) %>% Reduce("cbind", .) %>% as.data.frame()
colnames(pheno_df) <- names(format_list_b)
pheno_df$eid <- pheno_sel_b$eid
saveRDS(pheno_df, file = "02_ukb/pheno_sel_b_df.rds")

##
out_all <- lapply(names(format_list_b), function(PFIDx){
  #
  outx <- ukb_pheno.out(col_list = PFIDx,
                        type_c = NA,
                        PFID = PFIDx,
                        do_scale = F,
                        pheno_df = pheno_df,
                        eid_list = id_list)
  outx$PFID = PFIDx
  print(paste0("MGS: Phenotype file for ", PFIDx, " is OK!"))
  return(outx)
}) %>% Reduce("rbind", .) %>% as.data.frame()
colnames(out_all) <- c("Type", "Sex", "Sample", "PFID")
out_all$Sample <- as.integer(out_all$Sample )
out_all_df <- reshape2::dcast(out_all, 
                              PFID + Sex ~Type, 
                              value.var = "Sample")
fwrite2(out_all_df, file = "02_ukb/out_all_b_df.txt", sep = "\t")

