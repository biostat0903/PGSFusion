# Load packages
library(bigreadr)
library(dplyr)
library(plyr)

# Set parameters
PATH <- "/public/home/biostat04/Project/19_PGS_fusion/"
TRAIT <- "/public/home/Datasets/ukb/pheno/03_trait/participant.csv.gz"

# Load sample id
sample_id <- list.files(paste0(PATH, "sample_id/")) %>% 
  alply(., 1, function(ff) fread2(paste0(PATH, "sample_id/", ff))[, 1])
names(sample_id) <- list.files(paste0(PATH, "sample_id/")) %>% 
  gsub(".fam", "", .)

# Load disease
## cancer flied id
cancer_self_report <- c(paste0("p20001_i0_a", c(0:5)), paste0("p20001_i1_a", c(0:5)),
                        paste0("p20001_i2_a", c(0:5)), paste0("p20001_i3_a", c(0:5)))
cancer_type <- paste0("p40006_i", c(0: 21))
## diagnose
main_diagnose <- "p41202"
second_diagnose <- "p41204"
## non-cancer id
noncancer_self_report <- c(paste0("p20002_i0_a", c(0:33)), paste0("p20002_i1_a", c(0:33)),
                            paste0("p20002_i2_a", c(0:33)), paste0("p20002_i3_a", c(0:33)))
disease_id <- c(cancer_self_report,
                cancer_type, 
                main_diagnose, 
                second_diagnose, 
                noncancer_self_report)  
trait_eid <- fread2(TRAIT, select = "eid")
trait_binary_dat <- fread2(TRAIT, select = disease_id)
field_id <- fread2(paste0(PATH, "fieldID_binary.txt"))
field_id$ICD10 <- gsub(",", "\\|", field_id$ICD10)

# Identify cancer
trait_binary_dat_cancer <- trait_binary_dat[, c(cancer_self_report, cancer_type, 
                                                main_diagnose, second_diagnose)]
field_id_cancer <- field_id[field_id$Category == "Cancer", ]
trait_binary_sample_cancer <- alply(c(1: length(sample_id)), 1, function (ss){
  
  trait_binary_dat_cancer_s <- trait_binary_dat_cancer[match(sample_id[[ss]], trait_eid[, 1]), ]
  trait_binary_out <- aaply(c(1: nrow(field_id_cancer)), 1, function(tt) {

    ## select cases
    idx_1 <- alply(trait_binary_dat_cancer_s[, cancer_self_report], 2, function(ss){
      which(ss == field_id_cancer$`Self-reported Trait`[tt])
    }) %>% unlist
    idx_2 <- alply(trait_binary_dat_cancer_s[, cancer_type], 2, function(ss){
      grep(field_id_cancer$`ICD10`[tt], ss[,1])
    }) %>% unlist
    idx_3 <- grep(field_id_cancer$`ICD10`[tt], trait_binary_dat_cancer_s[, main_diagnose])
    idx_4 <- grep(field_id_cancer$`ICD10`[tt], trait_binary_dat_cancer_s[, second_diagnose])
    idx <- c(idx_1, idx_2, idx_3, idx_4) %>% unique
    ## define cases
    label_s <- rep(0, nrow(trait_binary_dat_cancer_s))
    label_s[idx] <- 1
    ## output
    write.table(label_s,
                file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/",
                              field_id_cancer$PFID[tt], ".txt"),
                col.names = F, row.names = F, quote = F)
    trait_binary_summ <- paste0(length(label_s), "(", sum(label_s), ":", sum(label_s==0), ")")
    return(trait_binary_summ)
  })
  return(trait_binary_out)
}) %>% do.call("cbind", .) 

# Identify noncancer
trait_binary_dat_noncancer <- trait_binary_dat[, c(noncancer_self_report, main_diagnose, 
                                                   second_diagnose)]
field_id_noncancer <- field_id[field_id$Category == "Non cancer", ]
trait_binary_sample_noncancer <- alply(c(1: length(sample_id)), 1, function (ss){
  
  trait_binary_dat_noncancer_s <- trait_binary_dat_noncancer[match(sample_id[[ss]], trait_eid[, 1]), ]
  trait_binary_out <- aaply(c(1: nrow(field_id_noncancer)), 1, function(tt) {
    
    ## select cases
    idx_1 <- alply(trait_binary_dat_noncancer_s[, noncancer_self_report], 2, function(ss){
      which(ss == field_id_noncancer$`Self-reported Trait`[tt])
    }) %>% unlist
    idx_2 <- grep(field_id_noncancer$`ICD10`[tt], trait_binary_dat_noncancer_s[, main_diagnose])
    idx_3 <- grep(field_id_noncancer$`ICD10`[tt], trait_binary_dat_noncancer_s[, second_diagnose])
    idx <- c(idx_1, idx_2, idx_3) %>% unique
    ## define cases
    label_s <- rep(0, nrow(trait_binary_dat_noncancer_s))
    label_s[idx] <- 1
    ## output
    write.table(label_s,
                file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/",
                              field_id_noncancer$PFID[tt], ".txt"),
                col.names = F, row.names = F, quote = F)
    trait_binary_summ <- paste0(length(label_s), "(", sum(label_s), ":", sum(label_s==0), ")")
    return(trait_binary_summ)
  })
  return(trait_binary_out)
}) %>% do.call("cbind", .) 


# Identify other traits
## Fracture
trait_binary_dat_nondis1 <- fread2(TRAIT, select = paste0("p2463_i", 0: 3))  %>%
  do.call(paste0, .)
trait_binary_sample_nondis1 <- aaply(c(1: length(sample_id)), 1, function (ss){
  
  trait_binary_dat_s <- trait_binary_dat_nondis1[match(sample_id[[ss]], trait_eid[, 1])]
  label_s <- rep(NA, length(trait_binary_dat_s))
  label_s[grep("Yes", trait_binary_dat_s)] <- 1
  label_s[grepl("No", trait_binary_dat_s) & !grepl("Yes", trait_binary_dat_s)] <- 0
  write.table(label_s, 
              file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/",
                           "PFIDB3001.txt"),
              col.names = F, row.names = F, quote = F)
  label_s_na <- label_s[!is.na(label_s)]
  trait_binary_summ <- paste0(length(label_s_na), "(", sum(label_s_na), ":", sum(label_s_na==0), ")")
  return(trait_binary_summ)
}) 

## tanning
trait_binary_dat_nondis2 <- fread2(TRAIT, select = paste0("p1727_i", 0: 2)) %>% 
  do.call(paste0, .)
trait_binary_sample_nondis2 <- aaply(c(1: length(sample_id)), 1, function (ss){
  
  trait_binary_dat_s <- trait_binary_dat_nondis2[match(sample_id[[ss]], trait_eid[, 1])]
  label_s <- rep(NA, length(trait_binary_dat_s))
  label_s[grep("Get moderately tanned|Get mildly or occasionally tanned|Never tan, only burn",
               trait_binary_dat_s) ] <- 0
  label_s[grep("Get very tanned", trait_binary_dat_s)] <- 1
  write.table(label_s, 
              file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/",
                            "PFIDB3002.txt"),
              col.names = F, row.names = F, quote = F)
   label_s_na <- label_s[!is.na(label_s)]
   trait_binary_summ <- paste0(length(label_s_na), "(", sum(label_s_na), ":", sum(label_s_na==0), ")")
   return(trait_binary_summ)
}) 

## type I balding
trait_binary_dat_nondis3 <- fread2(TRAIT, select = paste0("p2395_i", 0: 3)) %>% 
  do.call(paste0, .)
trait_binary_sample_nondis3 <- aaply(c(1: length(sample_id)), 1, function (ss){
  
  trait_binary_dat_s <- trait_binary_dat_nondis3[match(sample_id[[ss]], trait_eid[, 1])]
  label_s <- rep(NA, length(trait_binary_dat_s))
  label_s[grep("Pattern 1", trait_binary_dat_s)] <- 0
  label_s[grep("2|3|4", trait_binary_dat_s)] <- 1
  write.table(label_s, 
              file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/",
                            "PFIDB3003.txt"),
              col.names = F, row.names = F, quote = F)
  label_s_na <- label_s[!is.na(label_s)]
  trait_binary_summ <- paste0(length(label_s_na), "(", sum(label_s_na), ":", sum(label_s_na==0), ")")
  return(trait_binary_summ)
}) 

## Ever smoking
trait_binary_dat_nondis4 <- fread2(TRAIT, select = paste0("p20160_i", 0: 3)) %>% 
  do.call(paste0, .)
trait_binary_sample_nondis4 <- aaply(c(1: length(sample_id)), 1, function (ss){
  
  trait_binary_dat_s <- trait_binary_dat_nondis4[match(sample_id[[ss]], trait_eid[, 1])]
  label_s <- rep(NA, length(trait_binary_dat_s))
  label_s[grep("Yes", trait_binary_dat_s)] <- 1
  label_s[grepl("No", trait_binary_dat_s) & !grepl("Yes", trait_binary_dat_s)] <- 0
  write.table(label_s, 
              file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/",
                            "PFIDB3004.txt"),
              col.names = F, row.names = F, quote = F)
  label_s_na <- label_s[!is.na(label_s)]
  trait_binary_summ <- paste0(length(label_s_na), "(", sum(label_s_na), ":", sum(label_s_na==0), ")")
  return(trait_binary_summ)
}) 

## Snoring
trait_binary_dat_nondis5 <- fread2(TRAIT, select = paste0("p1210_i", 0: 3)) %>% 
  do.call(paste0, .)
trait_binary_sample_nondis5 <- aaply(c(1: length(sample_id)), 1, function (ss){
  
  trait_binary_dat_s <- trait_binary_dat_nondis5[match(sample_id[[ss]], trait_eid[, 1])]
  label_s <- rep(NA, length(trait_binary_dat_s))
  label_s[grep("Yes", trait_binary_dat_s)] <- 1
  label_s[grepl("No", trait_binary_dat_s) & !grepl("Yes", trait_binary_dat_s)] <- 0
  write.table(label_s, 
              file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/",
                            "PFIDB3005.txt"),
              col.names = F, row.names = F, quote = F)
  label_s_na <- label_s[!is.na(label_s)]
  trait_binary_summ <- paste0(length(label_s_na), "(", sum(label_s_na), ":", sum(label_s_na==0), ")")
  return(trait_binary_summ)
}) 

# Summarize data
trait_summary <- list(trait_binary_sample_cancer, trait_binary_sample_noncancer, 
                              trait_binary_sample_nondis1, trait_binary_sample_nondis2, 
                              trait_binary_sample_nondis3, trait_binary_sample_nondis4, 
                              trait_binary_sample_nondis5) %>% do.call("rbind" ,.)
colnames(trait_summary) <- names(sample_id)
write.table(trait_summary, file = paste0(PATH, "out_pheno/bianry_summary.txt"),
            row.names = F, quote = F, sep = "\t")
