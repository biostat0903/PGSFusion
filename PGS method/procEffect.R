#! /usr/bin/env Rscript

# load packages
library(dplyr)
library(plyr)
library(bigreadr)
library(optparse)

# input parameters
args_list = list(
  make_option("--method", type="character", default=NULL,
              help="INPUT: method", metavar="character"), 
  make_option("--esteff", type="character", default=NULL,
              help="INPUT: effet size", metavar="character"), 
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary statistic", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)
# opt=list(esteff = "/home/chencao_pgs/website/pgsfusion-server/job/d449824760b74d55a1cf29aa02fe7d71/MTPGS/esteff.txt",
#          summ = "/home/chencao_pgs/website/pgsfusion-server/job/14bb978f52a1468f9c9740c3e5bc8b85/trait1.txt")

esteff <- fread2(opt$esteff)
summ <- fread2(opt$summ)
colnames(summ) <- c("chr", "rs", "ps", "n_mis", "n_obs", "allele1", 
                    "allele0", "af", "beta", "se", "p_wald")


if (opt$method %in% c("DBSLMM", "mtPGS")){
  
  esteff_j <- inner_join(esteff, summ, by = join_by("V1" == "rs"))
  esteff_out <- data.frame(rs = esteff_j$V1,
                           a1 = esteff_j$V2,
                           beta = esteff_j$V4,
                           beta_raw = esteff_j$beta,
                           pos = esteff_j$ps)
}

if (opt$method == "megaPRS"){
  
  esteff_j <- inner_join(esteff, summ, by = join_by("Predictor" == "rs"))
  esteff_out <- data.frame(rs = esteff_j$Predictor,
                           a1 = esteff_j$A1,
                           beta = esteff_j[, 5],
                           beta_raw = esteff_j$beta,
                           pos = esteff_j$ps)
}

if (opt$method %in% c("SDPR", "SDPRX")){

  esteff_j <- inner_join(esteff, summ, by = join_by("V1" == "rs"))
  esteff_out <- data.frame(rs = esteff_j$V1,
                           a1 = esteff_j$V2,
                           beta = esteff_j$V3,
                           beta_raw = esteff_j$beta,
                           pos = esteff_j$ps)
}

if (opt$method == "PRSCS"){

  esteff_j <- inner_join(esteff, summ, by = join_by("V2" == "rs"))
  esteff_out <- data.frame(rs = esteff_j$V2,
                           a1 = esteff_j$V4,
                           beta = esteff_j$V6,
                           beta_raw = esteff_j$beta,
                           pos = esteff_j$ps)
}

if (opt$method == "SBayesRC"){
  
  esteff_j <- inner_join(esteff, summ, by = join_by("SNP" == "rs"))
  esteff_out <- data.frame(rs = esteff_j$SNP,
                           a1 = esteff_j$A1,
                           beta = esteff_j$BETA,
                           beta_raw = esteff_j$beta,
                           pos = esteff_j$ps)
}



outpath <- dirname(opt$esteff)
write.table(esteff_out, file = paste0(outpath,'/esteff.txt'),
            col.names = F,row.names = F, quote = F)
