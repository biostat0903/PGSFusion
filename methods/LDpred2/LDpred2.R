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
  make_option("--output", type="character", default=NULL,
              help="INPUT: output path", metavar="character"),
  make_option("--dat", type="character", default=NULL,
              help="INPUT: dat type", metavar="character"),
  make_option("--cov", type="character", default=NULL,
              help="INPUT: covariates", metavar="character"),
  make_option("--val_genotype", type="character", default=NULL,
              help="INPUT: validation genotype", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# Set parameters
TEMP_DIR <- '/root/reference/intermediate_file'
REF_PATH <- "/disk/reference_pgsfusion/LD_reference_3cM/EUR"
p_len <- 21

# Fit LDpred2 in each chromosome
lastout <- alply(c(1: 22), 1, function(chr){

  ref_str <- paste0(REF_PATH, '/LD_with_blocks_chr', chr)
  val_str <- opt$val_genotype
  if(file.exists(paste0(val_str, ".rds")) == F | file.exists(paste0(val_str, ".bk")) == F){
    if(file.exists(paste0(val_str, ".bk"))){
      system(paste0("rm ", val_str, ".bk"))
    }
    val_bed <- snp_readBed(paste0(val_str, ".bed"))
  }
  ## Process reference panel
  ref_corr <- readRDS(paste0(ref_str, ".rds"))
  ref_info <- readRDS(paste0(REF_PATH,'/map_hm3_plus.rds'))
  ref_info <- ref_info[which(ref_info$chr == chr),]
  ref_map <- data.frame(ref_info$chr,ref_info$rsid,ref_info$pos,
                        ref_info$a1,ref_info$a0)
  names(ref_map) <- c("chr", "marker.ID", "pos", "a1", "a0")
  ## Process summary statistics
  summstats <- fread2(opt$summ,
                      select =  c(1, 2, 3, 5, 6, 7, 8, 9, 10),header=T)
  colnames(summstats) <- c("chr", "rsid", "pos", "n_obs", 
                           "a1", "a0", "MAF", "beta", "beta_se")
  summstats$n_eff <- summstats$n_obs
  summstats$n_obs <- NULL
  summstats <- summstats[summstats[, 1] == chr, ]
  val_sub_str <- paste0(TEMP_DIR,"/val_sub-", as.numeric(as.POSIXlt(Sys.time())))
  snp_ref_summ_inter <- snp_match(summstats, ref_map, 
                                  match.min.prop = 0.05)
  snp_ref_summ_inter <- snp_ref_summ_inter[snp_ref_summ_inter$rsid == snp_ref_summ_inter$marker.ID, ]
  val_bed <- snp_attach(paste0(val_str, ".rds"))
  val_map <- val_bed$map[-3]
  names(val_map) <- c("chr", "marker.ID", "pos", "a1", "a0")
  snp_val_summ_inter <- snp_match(summstats, val_map, 
                                  match.min.prop = 0.05)
  snp_val_summ_inter <- snp_val_summ_inter[snp_val_summ_inter$rsid == snp_val_summ_inter$marker.ID, ]
  snp_inter <- intersect(snp_val_summ_inter$rsid, snp_ref_summ_inter$rsid)
  val_sub_bed <- snp_attach(snp_subset(val_bed,
                                       ind.col = which(val_bed$map$marker.ID%in%snp_inter),
                                       backingfile = val_sub_str))
  df_beta_h2 <- snp_ref_summ_inter[, c("rsid", "a1", "beta", "beta_se", "n_eff")]
  df_beta_h2 <- df_beta_h2[df_beta_h2$rsid %in% snp_inter, ]
  ref_pos <- which(ref_info$rsid %in% snp_inter)
  corr <- ref_corr[ref_pos,ref_pos]
  h2 <- snp_ldsc2(corr, df_beta_h2)
  cat("Heritability: ", h2[2], "\n")
  if (h2[2] < 0){
    
    beta_LDpred2 <- data.frame(df_beta_h2$rsid[df_beta_h2$rsid %in% snp_inter],
                               df_beta_h2$a1[df_beta_h2$rsid %in% snp_inter],
                               0)
  } else {
    
    y <- fread2(opt$val_phenotype,header=F)[,1]
    covar <- fread2(opt$cov,header=F)[, 1]
    p_seq <- signif(seq_log(1e-5, 1, length.out = p_len), 2)
    h_seq <- h2[2] * c(0.7, 1, 1.4)
    params <- expand.grid(p = p_seq, h2 = h_seq, 
                          sparse = c(FALSE, TRUE))
    beta_grid <- snp_ldpred2_grid(as_SFBM(corr), df_beta_h2,
                                  params, ncores = 1)
    val_sub_str2 <- paste0(TEMP_DIR,"/val2_sub-", as.numeric(as.POSIXlt(Sys.time())))
    if (any(is.na(y))){
      
      row_idx <- which(!is.na(y))
      val_sub_bed <- snp_attach(snp_subset(val_sub_bed, 
                                           ind.row = row_idx, 
                                           backingfile = val_sub_str2))
      y <- y[row_idx]
      if (exists('covar')){
        covar <- covar[row_idx]
      }
    }
    val_G <- val_sub_bed$genotypes
    beta_grid_sub <- beta_grid[df_beta_h2$rsid %in% snp_inter, ]
    pred_grid <- big_prodMat(val_G, beta_grid_sub)
    idx_na <- apply(pred_grid, 2, function(a) all(is.na(a))) | 
      apply(beta_grid, 2, function(a) all(is.na(a)))
    beta_grid_na <- beta_grid[, !idx_na]
    pred_grid_na <- pred_grid[, !idx_na]
    params_na <- params[!idx_na, ]
    
    if (all(idx_na)){
      
      beta_LDpred2 <- data.frame(df_beta_h2$rsid[df_beta_h2$rsid %in% snp_inter],
                                 df_beta_h2$a1[df_beta_h2$rsid %in% snp_inter],
                                 0)
    } else {
      if (opt$dat == "quantitative"){
        params_na[c("coef", "score")] <-
        big_univLinReg(big_copy(pred_grid_na), 
                       y)[c("estim", "score")]
	params_na$idx_val <- apply(pred_grid_na, 2, function(a) cor(a, y)^2)
      } else {
        params_na[c("coef", "score")] <-
        big_univLinReg(big_copy(pred_grid_na), 
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
      best_grid_nosp <- params_na %>%
        mutate(id = c(1: nrow(params_na))) %>%
        filter(!sparse) %>%
        arrange(desc(score)) %>%
        slice(1) %>%
        { beta_grid_na[, .$id] * .$coef }
      beta_LDpred2 <- data.frame(df_beta_h2$rsid[df_beta_h2$rsid %in% snp_inter],
                                 df_beta_h2$a1[df_beta_h2$rsid %in% snp_inter],
                                 best_grid_nosp[df_beta_h2$rsid %in% snp_inter])
    }
  }
  colnames(beta_LDpred2) <- c('rsid','a1','effect')
  system(paste0("rm ", val_sub_str, ".bk"))
  system(paste0("rm ", val_sub_str, ".rds"))
  if (file.exists(paste0(val_sub_str2, ".bk"))){
    system(paste0("rm ", val_sub_str2, ".bk"))
    system(paste0("rm ", val_sub_str2, ".rds"))
  }
  return(beta_LDpred2)
})

# Output
write.table(lastout, file = paste0(opt$output,'/esteff.txt'), 
            row.names = F, col.names = F, quote = F)
system(paste0('rm ',val_str,'.rds'))
system(paste0('rm ',val_str,'.bk'))

