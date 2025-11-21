# Load packages
library(plyr)
library(dplyr)
library(bigsnpr)
library(bigreadr)
library(optparse)

# Input parameters
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--val_phenotype", type="character", default=NULL,
              help="INPUT: phenotype value", metavar="character"),
  make_option("--val_genotype", type="character", default=NULL,
              help="INPUT: validation genotype", metavar="character"), 
  make_option("--reference", type="character", default=NULL,
              help="INPUT: reference panel", metavar="character"), 
  make_option("--dat", type="character", default=NULL,
              help="INPUT: dat type", metavar="character"),
  make_option("--ncase", type="integer", default=NULL,
              help="INPUT: case numeber", metavar="character"), 
  make_option("--ncontrol", type="integer", default=NULL,
              help="INPUT: control number", metavar="character"), 
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
#             reference = "/disk/reference_pgsfusion/EUR_UKB_ref/LD_matrix",
#             val_phenotype = "/disk/validationSet/phenotype/All/PFIDB024.txt",
#             val_genotype = "/disk/validationSet/genotype/All",
#             dat = "binary",
#             ncase = 10000,
#             ncontrol = 90000,
#             cov = "/disk/validationSet/coef/All/PFIDB024.txt",
#             output = "/public/home/biostat03/project/pgsfusionProject/tmp"
# )

# Set parameters
TEMP_DIR <- paste0("/root/reference/intermediate_file/", opt$jobid)
# TEMP_DIR="/public/home/biostat03/project/pgsfusionProject/tmp"
p_len <- 21
# p_len <- 3

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
colnames(sumstats) <- c("chr", "pos", "rsid", "a1", "a0", "freq", "beta", "beta_se", "n_eff")

# Fit model for each chromosome
esteff <- alply(c(1: 22), 1, function(CHR) {
  
  sumstats_chr <- sumstats[sumstats$chr == CHR, ]
  map_chr <- map[map$chr == CHR, ]
  df_beta <- as_tibble(snp_match(sumstats_chr, map_chr,
                                 return_flip_and_rev = TRUE))
  
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
  p_seq <- signif(seq_log(1e-5, 1, length.out = p_len), 2)
  # h_seq <- pmax(ldsc[2], 0.0001)
  # params <- expand.grid(p = p_seq, h2 = h_seq,
  #                       sparse = FALSE)
  h_seq <- pmax(ldsc[2], 0.0001) * c(0.7, 1, 1.4)
  params <- expand.grid(p = p_seq, h2 = h_seq,
                        sparse = c(FALSE, TRUE))
  beta_grid <- snp_ldpred2_grid(corr_sub, 
                                df_beta,
                                params, 
                                ncores = 1)
  ## validation
  y <- fread2(opt$val_phenotype, header = F)[, 1]
  if (!is.null(opt$cov))
    covar <- fread2(opt$cov, header = F)
  geno_chr <- paste0(opt$val_genotype, "/chr", CHR)
  if (!file.exists(paste0(geno_chr, ".rds")))
    
    snp_readBed(paste0(geno_chr, ".bed"))
  val_bed <- snp_attach(paste0(geno_chr, ".rds"))
  val_sub_str <- paste0(TEMP_DIR, "/val_sub_chr", CHR, "-", 
                        as.numeric(as.POSIXlt(Sys.time())))
  if (any(is.na(y))){
    
    row_idx <- which(!is.na(y))
    val_sub_bed <- snp_attach(snp_subset(val_bed, 
                                         ind.row = row_idx, 
                                         ind.col = df_beta$`_NUM_ID_`,
                                         backingfile = val_sub_str))
    y <- y[row_idx]
    # y=y-1
    if (exists('covar'))
      
      covar <- covar[row_idx, ]
  } else {
    
    val_sub_bed <- snp_attach(snp_subset(val_bed, 
                                         ind.col = df_beta$`_NUM_ID_`,
                                         backingfile = val_sub_str))
  }
  ## Select the best parameter
  # val_G <- snp_fastImputeSimple(val_sub_bed$genotypes)
  val_G <- val_sub_bed$genotypes
  pred_grid <- big_prodMat(val_G, beta_grid)
  idx_na <- apply(pred_grid, 2, function(a) all(is.na(a))) | 
    apply(beta_grid, 2, function(a) all(is.na(a)))
  beta_grid_na <- beta_grid[, !idx_na]
  pred_grid_na <- pred_grid[, !idx_na]
  params_na <- params[!idx_na, ]
  
  # Output effect 
  if (all(idx_na)){
    
    esteff_chr <- data.frame(rsid = df_beta$rsid,
                             a0 = df_beta$a0,
                             beta_est =0, 
                             beta_raw = df_beta$beta,
                             pos =  df_beta$pos)
  } else {
    
    if (opt$dat == "quantitative"){
      params_na[c("coef", "score")] <- big_univLinReg(big_copy(pred_grid_na), 
                                                      y, 
                                                      covar.train = as.matrix(covar))[c("estim", "score")]
      params_na$idx_val <- apply(pred_grid_na, 2, function(a) cor(a, y)^2)
    } else {
      params_na[c("coef", "score")] <- big_univLinReg(big_copy(pred_grid_na), 
                                                      y, 
                                                      covar.train = as.matrix(covar))[c("estim", "score")]
      params_na$idx_val  <- apply(pred_grid_na, 2, AUC, target = y)
    }
    ## Select parameters
    params_na <- params_na %>%
      mutate(sparsity = colMeans(beta_grid_na == 0, na.rm = T),
             id = c(1: nrow(params_na))) %>%
      arrange(desc(score)) %>%
      mutate_at(4:8, signif, digits = 3)
    best_grid <- params_na %>%
      mutate(id = c(1: nrow(params_na))) %>%
      filter(!sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      { beta_grid_na[, .$id] * .$coef }
    esteff_chr <- data.frame(rsid = df_beta$rsid, 
                             a1 = df_beta$a1,
                             beta_est = best_grid, 
                             beta_raw = df_beta$beta,
                             pos =  df_beta$pos)
    
    # Remove files
    system(paste0("rm -rf ", ref_sub_str, ".sbk"))
    system(paste0("rm -rf ", val_sub_str, ".bk"))
    system(paste0("rm -rf ", val_sub_str, ".rds"))
  }
  
  return(esteff_chr)
}) %>% do.call("rbind", .)

# Output whole genome
write.table(esteff, file = paste0(opt$output,'/esteff.txt'), 
            row.names = F, col.names = F, quote = F)
