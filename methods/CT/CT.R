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
  make_option("--ref_genotype", type="character", default=NULL,
              help="INPUT: reference panel", metavar="character"),
  make_option("--window", type="character", default=NULL,
              help="INPUT: window size", metavar="character"),
  make_option("--r2", type="character", default=NULL,
              help="INPUT: r2 value", metavar="character"),
  make_option("--plen", type="integer", default=NULL,
              help="INPUT: p value", metavar="character"),
  make_option("--dat", type="character", default=NULL,
              help="INPUT: dat type", metavar="character"),
  make_option("--val_phenotype", type="character", default=NULL,
              help="INPUT: phenotype qqnorm", metavar="character"), 
  make_option("--val_genotype", type="character", default=NULL,
              help="INPUT: validation genotype", metavar="character"),
  make_option("--cov", type="character", default=NULL,
              help="INPUT: covariates", metavar="character"),
  make_option("--output", type="character", default=NULL,
              help="INPUT: output path", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# Set global parameters
TEMP_DIR <- "/root/reference/intermediate_file"
                   
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
ref_str <- opt$ref_genotype
val_str <- opt$val_genotype
ref_sub_str <- paste0(TEMP_DIR, "/ref_panel", "_sub-", 
                      as.numeric(as.POSIXlt(Sys.time())))
val_sub_str <- paste0(TEMP_DIR, "/val_genotype", "_sub-",
                      as.numeric(as.POSIXlt(Sys.time())))

# Coordinate validation set, reference panel and summary statistics
## Coordinate validation set and reference panel
if(file.exists(paste0(val_str, ".rds")) == F | file.exists(paste0(val_str, ".bk")) == F){
  if(file.exists(paste0(val_str, ".bk"))){
    system(paste0("rm ", val_str, ".bk"))
  }
  val_bed <- snp_readBed(paste0(val_str, ".bed"))
}
ref_bed <- snp_attach(paste0(ref_str, ".rds"))
val_bed <- snp_attach(paste0(val_str, ".rds"))
val_snp <- fread2(paste0(val_str, ".bim"))
if(all(ref_bed$map$marker.ID == val_bed$map$marker.ID) == F){
  
  snp_inter <- intersect(ref_bed$map$marker.ID, val_bed$map$marker.ID)
  ref_bed <- snp_attach(snp_subset(ref_bed, 
                                   ind.col = which(ref_bed$map$marker.ID%in%snp_inter), 
                                   backingfile = ref_sub_str))
  val_bed <- snp_attach(snp_subset(val_bed,
                                   ind.col = which(val_bed$map$marker.ID%in%snp_inter),
                                   backingfile = val_sub_str))
}
val_n_snp <- dim(val_bed$genotypes)[2]
## Process reference panel
ref_G <- ref_bed$genotypes
ref_CHR <- ref_bed$map$chromosome
ref_POS <- ref_bed$map$physical.pos
ref_n_snp <- dim(ref_G)[2]
ref_map <- ref_bed$map[, -3]
names(ref_map) <- c("chr", "rsid", "pos", "a1", "a0")
# Process summary statistics
summstats <- fread2(opt$summ, select =  c(1, 2, 3, 7, 6, 9, 10),header=T)
colnames(summstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "se") 
t <- summstats$beta/summstats$se
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
summstats$pval <- ifelse(p_val == 0,
                         min(p_val[-which(p_val==0)]),
                         p_val)
# Coordinate reference panel and summary statistics
info_snp <- snp_match(summstats, ref_map)
beta <- rep(0, ref_n_snp)
lp_val <- rep(0, ref_n_snp)
beta[ref_map[, 2]%in%info_snp[, 5]] <- info_snp$beta
lp_val[ref_map[, 2]%in%info_snp[, 5]] <- -log10(info_snp$pval)

# Perform CT model
all_keep <- snp_grid_clumping(ref_G,
                              grid.thr.r2 = as.numeric(r2),
                              grid.base.size = as.numeric(dist),
                              infos.chr = ref_CHR,
                              infos.pos = ref_POS,
                              lpS = lp_val,
                              ncores = 1)

# Estimate PGS in validation set
## Load validation phenotype
y <- fread2(opt$val_phenotype, header=F)[,1]
covar <- fread2(opt$cov,header=F)[, 1]

val_bk_str2 <- paste0(TEMP_DIR, "/val_sub-", as.numeric(as.POSIXlt(Sys.time())))
if (any(is.na(y)) == T){
  
  idx <- which(!is.na(y))
  val_bed <- snp_attach(snp_subset(val_bed, ind.row = idx, 
                                   backingfile = val_bk_str2))
  y <- y[which(!is.na(y))]
  if (exists('covar')){
    covar <- covar[which(!is.na(y))]
  }
}

## Estimate PGS
ct_bk_str <- paste0(TEMP_DIR, "/ct2_sub-",
                     as.numeric(as.POSIXlt(Sys.time())))
multi_PRS <- snp_grid_PRS(val_bed$genotypes,
                          all_keep,
                          betas = beta,
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
c_idx <- c(1: ref_n_snp) %in% unlist(map(all_keep, max_prs$id))
t_idx <- c(1: ref_n_snp) %in% which(lp_val >= 0.999999*max_prs$thr.lp)
idx <- ifelse(c_idx == T & t_idx == T, T, F)
snp_sig_CT <- data.frame(ref_map$rsid[idx],
                         ref_map$a1[idx],
                         beta[idx])

#output
write.table(snp_sig_CT, file = paste0(opt$output,"/esteff.txt"),
            col.names = F, row.names = F, quote = F ,append = F)

# remove file
system(paste0("rm ", ref_sub_str, ".bk"))
system(paste0("rm ", ref_sub_str, ".rds"))
system(paste0("rm ", val_sub_str, ".bk"))
system(paste0("rm ", val_sub_str, ".rds"))
if(file.exists(paste0(val_bk_str2, ".bk"))){
  system(paste0("rm ", val_bk_str2, ".bk"))
  system(paste0("rm ", val_bk_str2, ".rds"))
}
if(file.exists(paste0(ct_bk_str, ".bk"))){
  system(paste0("rm ", ct_bk_str, ".bk"))
  system(paste0("rm ", ct_bk_str, ".rds"))
}
system(paste0('rm ',val_str,'.bk'))
system(paste0('rm ',val_str,'.rds'))
