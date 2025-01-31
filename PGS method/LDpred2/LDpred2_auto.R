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

# opt <- list(summ = "/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/trait1.txt",
#             reference = "/disk/reference_pgsfusion/EUR_UKB_ref/LD_matrix",
#             ncase = 10000,
#             ncontrol = 90000,
#             output = "/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/LDPRED_AUTO/")

# Set parameter
NCORES <- 1
TEMP_DIR <- paste0("/root/reference/intermediate_file/", opt$jobid)

# QC for summary statistics and hm3 file
map <- transmute(readRDS(paste0(opt$reference,"/map_hm3_plus.rds")),
                 chr = as.integer(chr), pos = pos, rsid, af_val = af_ref,
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
colnames(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "freq", "beta", "beta_se", "n_eff")

# Fit model for each chromosome
esteff <- alply(c(1: 22), 1, function(CHR){

  # CHR=22
  # cat("chr", CHR, "\n")
  sumstats_chr <- sumstats[sumstats$chr == CHR, ]
  map_chr <- map[map$chr == CHR, ]
  df_beta <- as_tibble(snp_match(sumstats_chr, map_chr, return_flip_and_rev = TRUE))

  ## Estimate heritability
  ldsc <- with(df_beta, snp_ldsc(ldsc, length(ldsc), chi2 = (beta / beta_se)^2,
                                 sample_size = n_eff, blocks = NULL))

  ## Estimate PGS for different parameters
  corr <- readRDS(paste0(opt$reference, 
                         "/LD_with_block_chr", CHR, ".rds"))
  corr_sub <- corr[df_beta$`_NUM_ID_`, df_beta$`_NUM_ID_`]
  ref_sub_str <- paste0(TEMP_DIR,"/ref_sub_chr", CHR, "-", 
                        as.numeric(as.POSIXlt(Sys.time())))
  corr_sub <- as_SFBM(corr_sub, ref_sub_str, compact = T)
  coef_shrink <- 0.95
  multi_auto <- snp_ldpred2_auto(corr_sub, df_beta, h2_init = pmax(ldsc[["h2"]], 0.0001),
                                 vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                 # vec_p_init = seq_log(1e-4, 0.2, length.out = 5),
                                 # burn_in = 50, num_iter = 20, report_step = 10,
                                 allow_jump_sign = TRUE, shrink_corr = coef_shrink,
                                 ncores = NCORES)
  ## Output effect 
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  
  if (all(is.na(range)) == F){
    
    keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3)
    beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    esteff_chr <- data.frame(rsid = df_beta$rsid.ss, 
                             a1 = df_beta$a1, 
                             beta_est = beta_auto,
                             beta_raw = df_beta$beta, 
                             pos = df_beta$pos)
  } else {
    
    esteff_chr <- data.frame(rsid = df_beta$rsid.ss, 
                             a1 = df_beta$a1, 
                             beta_est = 0,
                             beta_raw = df_beta$beta, 
                             pos = df_beta$pos)
  }
  ## Remove files
  system(paste0("rm -rf ", ref_sub_str, ".sbk"))
  return(esteff_chr)
  
}) %>% do.call("rbind", .)

# Output whole genome
write.table(esteff, file=paste0(opt$output,'/esteff.txt'),
            col.names=F, row.names=F, quote=F)

