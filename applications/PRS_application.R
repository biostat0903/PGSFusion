# Load packages
library(bigreadr)
library(dplyr)
library(ggplot2)
library(optparse)
library(trend)

# Input parameters
args_list = list(
  make_option("--phenocode", type="character", default=NULL,
              help="INPUT: phenotype code", metavar="character"),
  make_option("--sex", type="character", default=NULL,
              help="INPUT: sex type of phenotype", metavar="character"),
  make_option("--anc", type="character", default=NULL,
              help="INPUT: ancestry type of phenotype", metavar="character"),
  make_option("--out_path", type="character", default=NULL,
              help="INPUT: out_path file prefix", metavar="character")
)
opt_parser = OptionParser(option_list = args_list)
opt = parse_args(opt_parser)
anc <- opt$anc
phenocode <- opt$phenocode
out_path <- opt$out_path
sex <- opt$sex

## Set global variables
test_path <- paste0("/disk/testSet/", anc, '/')
code_path <- "/root/biosoft/get_picture/"
COV_USE <- c("Age", "Sex", paste0("PC", 1:20))
PGS_GROUP_LIST <- c(3:10)
PGS_SUB_LIST <- c(3:5)

##### analysis #####
### format PGS
## 0.1 set pheno, strata, and covar
pheno_path <- paste0(test_path, "phenotype/", sex, "/")
pheno_test <- fread2(paste0(pheno_path, phenocode, ".txt"))[, 1, drop = T]
type <- ifelse(length(unique(pheno_test[!is.na(pheno_test)])) == 2, "b", "c")
#
if (type == "c") {
  cov_test <- coef_cov_test <- NULL
  cov_sum_test <- rep(0, length(pheno_test))
} else {
  cov_test <- fread2(paste0(pheno_path, "cov.txt"))
  colnames(cov_test) <- read.table(paste0(pheno_path, "cov_var.txt"), header = F)[, 1, drop = T]
  
  coef_cov_test <- fread2(paste0("/disk/testSet/cov_beta/",sex, "/coef_cov_", phenocode, ".txt"))[, 1, drop = T]
  names(coef_cov_test) <- c("Intercept", colnames(cov_test))
  
  cov_test <- cov_test[, intersect(colnames(cov_test), COV_USE)]
  coef_cov_test <- coef_cov_test[c("Intercept", colnames(cov_test))]
  cov_sum_test <- as.matrix(cbind(1, cov_test)) %*% as.matrix(coef_cov_test)
}
#
strata_test <- readRDS(paste0(pheno_path, "strata_df.rds"))
strata_var <- colnames(strata_test)

## 0.2 sum PGS
all_pred_file <- list.files(out_path, 
                            pattern = "^pred_hm3.*\\.profile$", 
                            full.names = T)
PGS_test_df <- data.frame(eid = fread2(all_pred_file[1])[,1])
PGS_test_mat <- lapply(all_pred_file, function(all_pred_filex){
  predx <- fread2(all_pred_filex)
  if (all(!is.na(predx))){
    scorex <- predx[match(PGS_test_df$eid, predx$IID), "SCORESUM"]
  } else {
    scorex <- data.frame(SCORESUM=rep(0,nrow(predx)))[,1]
  }
  return(scorex)
}) %>% Reduce("cbind", .)

## 0.3 set PGS groups
PGS_test_df$PGS <- rowSums(PGS_test_mat)
PGS_group <- lapply(PGS_GROUP_LIST, function(n){
  cut(PGS_test_df$PGS, 
      breaks = quantile(PGS_test_df$PGS, probs = seq(0, 1, 1/n)), 
      labels = seq(n))
}) %>% Reduce("cbind", .) %>% as.data.frame()
colnames(PGS_group) <- paste0("G", PGS_GROUP_LIST)
PGS_test_df <- cbind(PGS_test_df,
                     PGS_group)

##### plot #####
source(paste0(code_path, "viz_fun.R"))
### 1. plot the whole performance of PGS
## 1.1 performance test
if (type == "c") {
  
  perform_all_plt <- dens.plt(pheno = pheno_test,
                              pred = PGS_test_df$PGS)
  
} else {
  
  perform_all_plt <- roc.plt(pheno = pheno_test,
                             pred = PGS_test_df$PGS + cov_sum_test)
  
}
## 1.2 plot
out_filex <- paste0(out_path, "Performance.png")
png(out_filex, height = 6, width = 7, units = "in",
    res = 300)
print(perform_all_plt)
dev.off()

### 2. plot eff across PGS group
PGS_trend_plt <- lapply(PGS_GROUP_LIST, function(nx){
  ## 2.1 format
  datt_plt <- data.frame(y = pheno_test,
                         PGS_g = PGS_test_df[[paste0("G", nx)]])
  if (!is.null(cov_test)) {
    datt_plt <- cbind(datt_plt,
                      cov_test)
  }
  ## 2.2 plot
  PGS_trend_pltx <- PGS.trend.plt(datt = datt_plt,
                                  n = nx,
                                  type = type)
  # save
  out_filex <- paste0(out_path, "Trend_PGS_group", nx, ".png")
  png(out_filex, height = 6, width = 7, units = "in",
      res = 300)
  print(PGS_trend_pltx)
  dev.off()
  
  message(paste0("MSG: Trend plot of PGS across ", nx, " groups is saved!"))
  return(out_filex)
})


### 3. plot {PGS} and performance (R2/AUC) by strata
strata_plt <- lapply(strata_var, function(varx){
  
  ## 3.1 format base data frame
  datt_strata <- data.frame(y = pheno_test,
                            PGS = PGS_test_df[["PGS"]])
  ## 3.2 add strata var
  datt_strata$Group <- strata_test[[varx]]
  
  ## 3.3 plot PGS in subgroups
  PGS_strata_pltx <- PGS.strata.plt(datt = datt_strata,
                                    strata_var = varx)
  # save plot
  out_PGS_filex <- paste0(out_path, "PGS_subgroup_by_", varx, ".png")
  plt_widx <- 3 + length(unique(datt_strata$Group)) * 1
  png(out_PGS_filex, height = 6, width = 6, units = "in",
      res = 300)
  print(PGS_strata_pltx)
  dev.off()
  message(paste0("MSG: PGS plot by ", varx, " is saved!"))
  
  ## 3.4 plot effect in subgroups
  datt_perform_strata <- perform.strata.test(pheno = datt_strata$y,
                                             pred = datt_strata$PGS + cov_sum_test,
                                             subgroup = datt_strata$Group,
                                             type = type)
  perform_strata_plt <- perform.strata.plt(datt = datt_perform_strata$perform_strata,
                                           test = datt_perform_strata$test,
                                           type = type,
                                           strata_var = varx)
  # save plot
  out_perform_filex <- paste0(out_path, "Performance_subgroup_by_", varx, ".png")
  plt_widx <- 3 + length(unique(datt_strata$Group)) * 0.8
  png(out_perform_filex, height = 5, width = plt_widx, units = "in",
      res = 300)
  print(perform_strata_plt)
  dev.off()
  message(paste0("MSG: Performance plot grouped by ", varx, " is saved!"))
  
  return(c(out_PGS_filex, out_perform_filex))
})


### 4. PGS effect (Beta/OR) by strata (with interaction test)
interact_plt <- lapply(strata_var, function(varx){
  
  ## 4.1 format base data frame
  datt_strata_plt <- data.frame(y = pheno_test,
                                PGS = PGS_test_df[["PGS"]])
  if (!is.null(cov_test)) {
    strata_cov_rmx <- ifelse(grepl("^Age", varx), "Age", varx)
    strata_covx <- setdiff(colnames(cov_test), strata_cov_rmx)
    datt_strata_plt <- cbind(datt_strata_plt, 
                             cov_test[, strata_covx])
  }
  ## 4.2 add strata var
  datt_strata_plt$Group <- strata_test[[varx]]
  
  ## 4.3 plot effect in subgroups
  effect_strata_pltx <- effect.strata.plt(datt = datt_strata_plt,
                                          strata_var = varx,
                                          type = type)
  # save
  out_effect_filex <- paste0(out_path, "Interaction_with_", varx, ".png")
  png(out_effect_filex, height = 6, width = 6, units = "in",
      res = 300)
  print(effect_strata_pltx)
  dev.off()
  message(paste0("MSG: Interaction plot with ", varx, " is saved!"))
  
  return(out_effect_filex)
})

### 4. Subgroup analysis of PGS group and strata
sub_pgsg_plt <- lapply(strata_var, function(varx){
  
  ## 4.1 format base data frame
  datt_sub_plt <- data.frame(y = pheno_test,
                             Group = strata_test[[varx]])
  if (!is.null(cov_test)) {
    
    strata_cov_rmx <- ifelse(grepl("^Age", varx), "Age", varx)
    strata_covx <- setdiff(colnames(cov_test), strata_cov_rmx)
    datt_sub_plt <- cbind(datt_sub_plt, 
                          cov_test[, strata_covx])
    
  }
  ## 4.2 plot for each PGS group setting
  sub_pltx <- lapply(PGS_SUB_LIST, function(nx){
    
    # add PGS group var
    datt_sub_plt$PGS_g <- PGS_test_df[[paste0("G", nx)]] %>% 
      paste0("PGS_group", .) %>%
      factor(., levels = paste0("PGS_group", 1:nx))
    
    ## subgroup analysis
    datt_sub_pltx <- subset(datt_sub_plt, rowSums(is.na(datt_sub_plt)) == 0)
    datt_sub_pltx[[varx]] <- datt_sub_pltx$Group
    ## sub Group by PGS_g
    sub_plt1 <- subgroup.plt(datt = datt_sub_pltx,
                              type = type,
                              sub_var = varx,
                              eff_var = "PGS_g")
    # save
    out_subgroup_prefix1 <- paste0(out_path, "Subgroup_of_PGS_g", nx, "_by_", varx)
    plt_widx <- 4 + max(nlevels(datt_sub_plt$Group), nx) * 1.3
    png(paste0(out_subgroup_prefix1, ".png"), 
        height = 5, width = plt_widx, units = "in", res = 300)
    print(sub_plt1$plt)
    dev.off()
    write.table(sub_plt1$tab, file = paste0(out_subgroup_prefix1, ".txt"),
                sep = "\t", col.names = T, row.names = F, quote = F)
    message(paste0("MSG: Subgroup analysis plot of PGS in group ", nx, 
                   " by ", varx, " is saved!"))
    
    ## sub PGS_g by Group
    sub_plt2 <- subgroup.plt(datt = datt_sub_pltx,
                             type = type,
                             sub_var = "PGS_g",
                             eff_var = varx)
    # save
    out_subgroup_prefix2 <- paste0(out_path, "Subgroup_of_", varx, "_by_PGS_g", nx)
    plt_widx <- 4 + max(nlevels(datt_sub_plt$Group), nx) * 1.3
    png(paste0(out_subgroup_prefix2, ".png"), 
        height = 5, width = plt_widx, units = "in", res = 300)
    print(sub_plt2$plt)
    dev.off()
    
    write.table(sub_plt2$tab, file = paste0(out_subgroup_prefix2, ".txt"),
                sep = "\t", col.names = T, row.names = F, quote = F)
    message(paste0("MSG: Subgroup analysis plot of ", varx, 
                   " by PGS in group ", nx, " is saved!"))
    
    return(c(out_subgroup_prefix1, out_subgroup_prefix2))
    
  }) %>% unlist
  
  return(sub_pltx)
})
