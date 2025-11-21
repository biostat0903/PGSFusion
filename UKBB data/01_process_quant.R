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
field_id <- fread2(paste0(PATH, "fieldID_quant.txt"))

# Load traits
trait_eid <- fread2(TRAIT, select = "eid")
trait_quant_dat <- field_id$Code[!grepl("\\/", field_id$Code)] %>%
  fread2(TRAIT, select = .)
trait_quant_dat_ext <- field_id$Code[grepl("\\/", field_id$Code)] %>% 
  strsplit(., "\\/") %>%
  llply(., function(pp) {
    
    dat_tmp <- trait_quant_dat[, match(pp, colnames(trait_quant_dat))]
    return(dat_tmp[, 1] / dat_tmp[, 2])
  }) %>% do.call("cbind", .) 
colnames(trait_quant_dat_ext) <- field_id$Code[grepl("\\/", field_id$Code)] 
trait_quant_dat <- cbind.data.frame(trait_quant_dat, trait_quant_dat_ext)
trait_quant_dat_num <- apply(trait_quant_dat, 2, as.numeric) 

# Output traits
trait_quant_sample <- alply(c(1: length(sample_id)), 1, function (ss){

  trait_quant_dat_scale <- trait_quant_dat_num[match(sample_id[[ss]], trait_eid[, 1]), ] %>%
              apply(., 2, scale)
  trait_quant_summ <- apply(trait_quant_dat_scale, 2, function(tt) sum(!is.na(tt))) 
  trait_quant_out <- aaply(c(1: ncol(trait_quant_dat_scale)), 1, function(tt){
  
    code_s <- colnames(trait_quant_dat_scale)[tt]
    write.table(trait_quant_dat_scale[, tt], 
                file = paste0(PATH, "out_pheno/", names(sample_id)[ss], "/", 
                              field_id$PFID[field_id$Code == code_s], ".txt"),
                col.names = F, row.names = F, quote = F)
    return(tt)
  })
  return(trait_quant_summ)
}) %>% do.call("cbind", .) 

# Summarize data
trait_quant_sample_s <- trait_quant_sample[match(field_id$Code, row.names(trait_quant_sample)), ]
colnames(trait_quant_sample_s) <- names(sample_id)
trait_summary <- cbind(field_id, trait_quant_sample_s)

write.table(trait_summary, file = paste0(PATH, "out_pheno/quant_summary.txt"),
            row.names = F, quote = F, sep = "\t")


