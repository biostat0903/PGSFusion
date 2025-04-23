# Load packages
library(bigsnpr)
library(bigreadr)
library(dplyr)
library(optparse)
library(plyr)

# Input parameters
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--reference", type="character", default=NULL,
              help="INPUT: reference panel", metavar="character"), 
  make_option("--ncase", type="integer", default=NULL,
              help="INPUT: case numeber", metavar="character"), 
  make_option("--ncontrol", type="integer", default=NULL,
              help="INPUT: control number", metavar="character"), 
  make_option("--output", type="character", default=NULL,
              help="INPUT: output path", metavar="character"),
  make_option("--jobid", type="character", default=NULL,
              help="jobid", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85//trait1.txt",
#             reference = "/disk/reference_pgsfusion/EUR_UKB_ref/LD_matrix",
#             ncase = 21982,
#             ncontrol = 41944,
#             output = "/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/LDPRED_AUTO/")

# Set parameter
NCORES <- 1
TEMP_DIR <- paste0("/root/reference/intermediate_file/", opt$jobid)

# QC for summary statistics and hm3 file
map <- transmute(readRDS(paste0(opt$reference,"/map_hm3_plus.rds")),
                 chr = as.integer(chr), pos = pos, af_val = af_ref,
                 a0 = a0, a1 = a1, ldsc = ldsc)
sumstats <- fread2(opt$summ, header = T)
if (opt$ncase != 0){
  
  neff <- 4 / (1 / opt$ncase + 1 / opt$ncontrol)
  sumstats <- sumstats[, c(1, 3, 2, 6, 7, 8, 9, 10)]
  sumstats$n_eff <- 4 / (1 / opt$ncase + 1 / opt$ncontrol)
}else{
  
  neff <- sumstats[1, 4] + sumstats[1, 5]
  sumstats <- sumstats[, c(1, 3, 2, 6, 7, 8, 9, 10, 5)]
}
colnames(sumstats) <- c("chr", "pos", "rsid", "a1", "a0", "freq", "beta", "beta_se", "n_eff")

# QC for SNPs
info_snp <- snp_match(sumstats, map)
sd_val <- with(info_snp, sqrt(2 * af_val * (1 - af_val)))
if(opt$ncase == 0 & opt$ncontrol == 0){
  
  sd_y_est = median(sd_val * info_snp$beta_se * sqrt(info_snp$n_eff))
  sd_ss = with(info_snp, sd_y_est / sqrt(n_eff * beta_se^2))
} else {
  
  sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
}
is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
cat("Whole genome has", sum(is_bad), "bed SNPs.\n")

# Estimate heritability
sumstats <- info_snp[!is_bad, ]
ldsc <- with(sumstats, snp_ldsc(ldsc, length(ldsc), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL))
if(ldsc[["h2"]] < 0.05) ldsc[["h2"]] <- 0.05
cat("Whole genome hertiability is ", ldsc[["h2"]] , "\n")

# Create genome-wide sparse LD matrix
ref_sub_str <- paste0(TEMP_DIR,"/ref_sub-", as.numeric(as.POSIXlt(Sys.time())))
for (chr in c(1: 22)) {

  ind.chr <- which(sumstats$chr == chr)
  ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
  ind.chr3 <- match(ind.chr2, which(map$chr == chr))
  corr0 <- readRDS(paste0(opt$reference, '/LD_with_block_chr', chr, '.rds'))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    
    corr <- as_SFBM(corr0, ref_sub_str, compact = TRUE)
  } else {
    
    corr$add_columns(corr0, nrow(corr))
  }
}

# Fit model for whole genome
coef_shrink <- 0.95
multi_auto <- snp_ldpred2_auto(corr, sumstats, h2_init = ldsc[["h2"]],
                               vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                               allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                               ncores = NCORES)

# Output
cnd <- sapply(multi_auto, function(auto) sum(is.na(auto$beta_est)) == length(auto$beta_est)) %>%
  any(. == rep(TRUE, 30))
if (cnd){
  
  ##
  cat("LDpred2-auto fails to construct PGS!")
  esteff <- data.frame(rsid = sumstats$rsid, 
                       a1 = sumstats$a1, 
                       beta_est = 0,
                       beta_raw = sumstats$beta, 
                       pos = sumstats$pos)
} else {
  
  ## Filter bed chain 
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  esteff <- data.frame(rsid = sumstats$rsid, 
                       a1 = sumstats$a1, 
                       beta_est = beta_auto,
                       beta_raw = sumstats$beta, 
                       pos = sumstats$pos)
}
write.table(esteff, file=paste0(opt$output,'/esteff.txt'),
            col.names=F, row.names=F, quote=F)

# Remove temp files
system(paste0("rm -rf ", ref_sub_str, ".sbk"))

