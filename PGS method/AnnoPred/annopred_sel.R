#!/bin/bash
library(bigreadr)
library(optparse)

# Parameter setting
args_list <- list(
  make_option("--output", type = "character", default = NULL,
              help = "INPUT: a output path", metavar = "character")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

# opt <- list(output='/home/chencao_pgs/website/pgsfusion-server/job/f4de6a6da78a4c09a727d3647de9e52b/ANNOPRED/')

# Set index: auc or cor
index <- fread2(paste0(opt$output, "/test_output/select_p.txt"))[, 1]

# Select best model
model <- c('h2_inf','h2_non_inf','pT_inf','pT_non_inf')
pval <- c("1.0", "0.3", "0.1", "0.03", "0.01", "0.003", "0.001", 
           "0.0003", "0.0001", "3e-05", "1e-05")
model_label <- rep(model, length(pval))
pval_label <- rep(pval, each = length(model))
slc_p <- data.frame(model_label, pval_label, index)

if (sum(is.na(slc_p)) == nrow(slc_p)){
  
  stop("All AUC/R2 are na!")
} else {
  
  ind_max <- slc_p[which.max(slc_p$index), ]
}

# Output best beta
cp_cmd <- paste0('cp ', opt$output, '/test_output/test_', ind_max[,1], 
                 '_betas_', ind_max[,2], '.txt ', opt$output, '/')
system(cp_cmd)
mv_cmd <- paste0('mv ', opt$output, '/test_', ind_max[,1], '_betas_',
                 ind_max[,2],'.txt ', opt$output, '/esteff.txt')
system(mv_cmd)

