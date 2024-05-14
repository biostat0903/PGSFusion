#
library(bigreadr)
library(dplyr)
library(stringr)

project_path <- "/public/home/biostat04/Project/19_PGS_fusion/"
setwd(project_path)
use_cov_all <- c("Age", "Sex", "BMI", paste0("PC", 1:20))
## format code list
code_df_c <- fread2("02_ukb/pheno_code_c.txt")
code_list_c <- lapply(1:nrow(code_df_c), function(x){
  #
  code_cx <- strsplit(code_df_c$UKBB_code[x], "\\/") %>% 
    unlist() %>% 
    as.integer() %>% 
    paste0(., "-")
  #
  if (code_df_c$BMI[x] == 0) {
    cov_usex <- setdiff(use_cov_all, "BMI")
  } else {
    cov_usex <- use_cov_all
  }
  #
  if (length(code_cx) == 2) {
    typex <- "1/2"
  } else {
    typex <- NA
  }
  return(list("code_c" = code_cx,
              "cov_use" = cov_usex,
              "type" = typex))
})
names(code_list_c) <- code_df_c$PFID
all_code_c <- lapply(code_list_c, function(code_list_cx){
  code_list_cx[["code_c"]]
}) %>% unlist %>% unique

## load ukb data
ukb_header <- fread2("/public/home/Datasets/ukb/ukb47503_header.txt", header = F) %>% 
  t() %>% as.vector()
select_var_list <- lapply(all_code_c, function(codex){
  ukb_header[grep(paste0("^", codex), ukb_header)]
})

names(select_var_list) <- all_code_c

id_list <- readRDS("02_ukb/id_list.rds")
pheno_sel_c <- fread2("02_ukb/pheno_all.txt.gz", 
                      select = c("eid", unlist(select_var_list) %>% as.vector))

## define base phenotype 
# (use the first record for individuals with multiple records)
pheno_df <- lapply(select_var_list, function(varx){
  pheno_sel_cx <- pheno_sel_c[, varx, drop = F]
  phenox <- apply(pheno_sel_cx, 1, function(x){
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    } else {
      return(x[1])
    }
  }) %>% unlist()
  return(phenox)
}) %>% Reduce("cbind", .) %>% as.data.frame()
colnames(pheno_df) <- paste0("code_", all_code_c)
pheno_df$eid <- pheno_sel_c$eid
saveRDS(pheno_df, file = "02_ukb/pheno_sel_c_df.rds")

## save for each PFID
source("code/pheno_def_fun.R")
out_all <- lapply(names(code_list_c), function(PFIDx){
  #
  code_list_cx <- code_list_c[[PFIDx]]
  code_listx <- paste0("code_", code_list_cx[["code_c"]])
  type_c <- code_list_cx[["type"]]
  outx <- ukb_pheno.out(col_list = code_listx,
                        type_c = type_c,
                        PFID = PFIDx,
                        do_scale = T,
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
fwrite2(out_all_df, file = "02_ukb/out_all_c_df.txt", sep = "\t")




