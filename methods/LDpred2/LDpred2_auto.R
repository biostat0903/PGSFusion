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
   make_option("--output", type="character", default=NULL,
               help="INPUT: output path", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# Set parameter
NCORES <- 1

# Matching and QC
map <- transmute(readRDS(paste0(opt$reference,"/map_hm3_plus.rds")),
                 chr = as.integer(chr), pos = pos, rsid, af_val = af_UKBB,
                 a0 = a0, a1 = a1)
sumstats <- fread2(opt$summ, header = T)
sumstats <- sumstats[, c(1, 3, 2, 6, 7, 8, 9, 10, 5)]
colnames(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "freq", "beta", "beta_se", "n_eff")
df_beta <- as_tibble(snp_match(sumstats, map,
                               return_flip_and_rev = TRUE))

# Estimate heritability of LDSC
corr <- runonce::save_run({

  for (chr in 1:22) {

    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map$chr == chr))
    corr_chr <- readRDS(paste0(opt$reference,"/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
    if (chr == 1) {
      
      corr <- as_SFBM(corr_chr, paste0(opt$output,"/corr"), compact = TRUE)
    } else {
      
      corr$add_columns(corr_chr, nrow(corr))
    }
  }
  corr
}, file = paste0(opt$output,"/corr.rds"))
ldsc <- snp_ldsc2(corr, df_beta, blocks = 200,
                  intercept = NULL, ncores = NCORES)

# Fit LDpred2-auto model
coef_shrink <- 0.95
multi_auto <- runonce::save_run(
  snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                   vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                   burn_in = 500, num_iter = 500, report_step = 20,
                   allow_jump_sign = TRUE, shrink_corr = coef_shrink,
                   ncores = NCORES),
  file = paste0(opt$output,"/mod.rds")
)

# Filter bad results
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3)
postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp))
esteff <- data.frame(df_beta$rsid.ss, df_beta$a0, postp)
write.table(esteff, file=paste0(opt$output,'/esteff.txt'),
            col.names=F, row.names=F, quote=F)

# Remove files
system(paste0('rm ',opt$output,'/corr*'))
system(paste0('rm ',opt$output,'/mod*'))
