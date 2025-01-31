#! Load packages
library(plyr)
library(dplyr)
library(bigsnpr)
library(bigreadr)
library(optparse)
library(XPASS)
library(data.table)

## Input parameters
args_list = list(
  make_option("--summ1", type="character", default=NULL,
              help="INPUT: gemma file1", metavar="character"), 
  make_option("--summ2", type="character", default=NULL,
              help="INPUT: gemma file2", metavar="character"),
  make_option("--val_geno", type="character", default=NULL,
              help="INPUT: val genotype file", metavar="character"),
  make_option("--ref1", type="character", default=NULL,
              help="INPUT: ref file1", metavar="character"), 
  make_option("--ref2", type="character", default=NULL,
              help="INPUT: ref file2", metavar="character"), 
  make_option("--output", type="character", default=1,
              help="INPUT: output path", metavar="character"),
  make_option("--pc1", type="character", default=1,
              help="INPUT: pc1 file", metavar="character"), 
  make_option("--pc2", type="character", default=NULL,
              help="INPUT: pc2 file", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# Fit XPASS
xpass <- XPASS( file_z1 = opt$summ1 , file_z2 = opt$summ2 ,
                file_ref1 = opt$ref1 , file_ref2 = opt$ref2 ,
                file_cov1 = opt$pc1 , file_cov2 = opt$pc2 ,
                sd_method = "LD_block" , compPosMean = T,
                file_out = paste0(opt$output , '/esteff_XPASS'))

# Format output
esteff <- xpass$mu
esteff_auxiliary <- esteff[, c(2, 4, 8, 3)]
esteff_target <- esteff[, c(2, 4, 9, 3)]
write.table(esteff_auxiliary, col.names=F, row.names=F, quote=F, 
            file=paste0(opt$output,'/esteff_auxiliary.txt'))
write.table(esteff_target, col.names=F, row.names=F, quote=F, 
            file=paste0(opt$output,'/esteff.txt'))

# Remove files
system(paste0('rm ',opt$output,'/esteff_XPASS.log'))
system(paste0('rm ',opt$output,'/esteff_XPASS_PosteriorMean.txt'))
system(paste0('rm ',opt$output,'/summ1.txt'))
system(paste0('rm ',opt$output,'/summ2.txt'))
system(paste0('rm ',opt$output,'/summ11.txt'))
system(paste0('rm ',opt$output,'/summ21.txt'))