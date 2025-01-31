# Load packages
library(plyr)
library(stringr)
library(bigstatsr)
library(bigsnpr)
library(bigreadr)
library(tidyverse)
library(optparse)

# Input parameters
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--val_phenotype", type="character", default=NULL,
              help="INPUT: phenotype qqnorm", metavar="character"), 
  make_option("--val_genotype", type="character", default=NULL,
              help="INPUT: validation genotype", metavar="character"),
  make_option("--reference", type="character", default=NULL,
              help="INPUT: reference panel", metavar="character"),
  make_option("--window", type="character", default=NULL,
              help="INPUT: window size", metavar="character"),
  make_option("--r2", type="character", default=NULL,
              help="INPUT: r2 value", metavar="character"),
  make_option("--plen", type="integer", default=NULL,
              help="INPUT: p value", metavar="character"),
  make_option("--dat", type="character", default=NULL,
              help="INPUT: dat type", metavar="character"),
  make_option("--cov", type="character", default=NULL,
              help="INPUT: covariates", metavar="character"),
  make_option("--output", type="character", default=NULL,
              help="INPUT: output path", metavar="character"),
  make_option("--jobid", type="character", default=NULL,
              help="jobid", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/trait1.txt",
#             val_phenotype = "/disk/validationSet/phenotype/All/PFIDB024.txt",
#             val_genotype = "/disk/validationSet/genotype/All/merge_imp",
#             reference = "/disk/reference_pgsfusion/EUR_UKB_ref",
#             window = 50,
#             r2 = 0.1,
#             plen = 3,
#             dat = "binary",
#             cov = "/disk/validationSet/coef/All/PFIDB024.txt",
#             output = "/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85",
#             jobid = "14bb978f52a1468f9c9740c3e5bc8b85")

# Set parameters
TEMP_DIR <- paste0("/root/reference/intermediate_file/", opt$jobid)

# Format parameters
if (str_detect(opt$window,',')){
  
  dist <- unlist(strsplit(opt$window, ","))
}else{
  
  dist <- as.numeric(opt$window)
}
if (str_detect(opt$r2,',')){
  
  r2 <- unlist(strsplit(opt$r2, ","))
}else{
  
  r2 <- as.numeric(opt$r2)
}

# QC for summary statistics and hm3 file
map <- transmute(readRDS(paste0(opt$reference,"/map_hm3_plus.rds")),
                 chr = as.integer(chr), pos = pos, rsid, af_val = af_ref,
                 a0 = a0, a1 = a1, ldsc = ldsc)
sumstats <- fread2(opt$summ, header = T)
sumstats <- sumstats[, c(1, 3, 2, 6, 7, 8, 9, 10, 5)]
t <- sumstats$beta/sumstats$se
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
sumstats$pval <- ifelse(p_val == 0, min(p_val[-which(p_val==0)]), p_val)
colnames(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "freq", 
                        "beta", "beta_se", "n_eff", "pval")
df_beta <- as_tibble(snp_match(sumstats, map,
                               return_flip_and_rev = TRUE))
lp_val <- -log10(df_beta$pval)

# Perform CT model
if (!file.exists(paste0(opt$reference, "/hm3_imp/merge_imp.rds")))
  
  snp_readBed(paste0(opt$reference, "/hm3_imp/merge_imp.bed"))
ref_bed <- snp_attach(paste0(opt$reference, "/hm3_imp/merge_imp.rds"))
ref_sub_str <- paste0(TEMP_DIR, "/ref_sub-", as.numeric(as.POSIXlt(Sys.time())))
ref_sub_bed <- snp_attach(snp_subset(ref_bed, 
                                     ind.col = df_beta$`_NUM_ID_`, 
                                     backingfile = ref_sub_str))
# ref_G <- snp_fastImputeSimple(ref_sub_bed$genotypes)
ref_G <- ref_sub_bed$genotypes
all_keep <- snp_grid_clumping(ref_G,
                              grid.thr.r2 = as.numeric(r2),
                              grid.base.size = as.numeric(dist),
                              infos.chr = df_beta$chr,
                              infos.pos = df_beta$pos,
                              lpS = lp_val,
                              ncores = 1)

# Estimate PGS in validation set
## Load validation phenotype
y <- fread2(opt$val_phenotype, header = F)[, 1]
covar <- fread2(opt$cov, header = F)[, 1]
if (!file.exists(paste0(opt$val_genotype, "merge_imp.rds")))
  snp_readBed(paste0(opt$val_genotype, "merge_imp.bed"))
val_bed <- snp_attach(paste0(opt$val_genotype, "merge_imp.rds"))
val_sub_str <- paste0(TEMP_DIR, "/val_sub-", as.numeric(as.POSIXlt(Sys.time())))
if (any(is.na(y))){
  
  row_idx <- which(!is.na(y))
  val_sub_bed <- snp_attach(snp_subset(val_bed, 
                                       ind.row = row_idx, 
                                       ind.col = df_beta$`_NUM_ID_`,
                                       backingfile = val_sub_str))
  y <- y[row_idx]
  if (exists('covar'))
    
    covar <- covar[row_idx]
} else {
  
  val_sub_bed <- snp_attach(snp_subset(val_bed, 
                                       ind.col = df_beta$`_NUM_ID_`,
                                       backingfile = val_sub_str))
}

## Estimate PGS
ct_bk_str <- paste0(TEMP_DIR, "/ct_sub-", as.numeric(as.POSIXlt(Sys.time())))
# val_G <- snp_fastImputeSimple(val_sub_bed$genotypes)
val_G <- val_sub_bed$genotypes
multi_PRS <- snp_grid_PRS(val_G,
                          all_keep,
                          betas = df_beta$beta,
                          lpS = lp_val,
                          n_thr_lpS = as.numeric(opt$plen),
                          backingfile = ct_bk_str,
                          ncores = 1)
nn <- nrow(attr(all_keep, "grid"))
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = c(1:nn)) %>%
  unnest(cols = "thr.lp")
s <- nrow(grid2)

# Select best parameter combination
if (opt$dat == "quantitative"){
  
  grid2$valIdx <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train, covar.train) {
    single_PRS <- rowSums(X[, ind + s * (0:21)] + covar.train)
    return(cor(single_PRS, y.train)^2)
  },
  ind = 1:s,
  s = s,
  y.train = y,
  covar.train = as.vector(covar),
  a.combine = 'c',
  block.size = 1,
  ncores = 1
  )
} else{
  
  grid2$valIdx <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train, covar.train) {
    single_PRS <- rowSums(X[, ind + s * (0:21)] + covar.train)
    return(bigstatsr::AUC(single_PRS, y.train))
  },
  ind = 1:s,
  s = s,
  y.train = y,
  covar.train = as.vector(covar),
  a.combine = 'c',
  block.size = 1,
  ncores = 1
  )
}
max_prs <- grid2 %>% arrange(desc(valIdx)) %>% dplyr::slice(1)
c_idx <- c(1: nrow(df_beta)) %in% unlist(map(all_keep, max_prs$id))
t_idx <- c(1: nrow(df_beta)) %in% which(lp_val >= 0.999999*max_prs$thr.lp)
idx <- ifelse(c_idx == T & t_idx == T, T, F)
snp_sig_CT <- data.frame(rsid = df_beta$rsid[idx],
                         a1 = df_beta$a1[idx],
                         beta = df_beta$beta[idx],
                         pos = df_beta$pos[idx])

# output
write.table(snp_sig_CT, file = paste0(opt$output,"/esteff.txt"),
            col.names = F, row.names = F, quote = F)

# remove file
# system(paste0("rm ", ref_sub_str, ".bk"))
# system(paste0("rm ", ref_sub_str, ".rds"))
# system(paste0("rm ", val_sub_str, ".bk"))
# system(paste0("rm ", val_sub_str, ".rds"))
# system(paste0("rm ", ct_bk_str, ".bk"))
# system(paste0("rm ", ct_bk_str, ".rds"))
system(paste0("rm -r ", TEMP_DIR))