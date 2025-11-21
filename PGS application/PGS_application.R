# Load packages
library(bigreadr)
library(dplyr)
library(ggplot2)
library(pROC)
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

# opt <- list(phenocode = "PFIDB2006",
#            sex = "All",
#            anc = "EUR",
#            out_path = "/home/chencao_pgs/website/pgsfusion-server/job/f4de6a6da78a4c09a727d3647de9e52b/PRSCS_AUTO")

## Set global variables
TEST_PATH <- paste0("/disk/testSet/", opt$anc, '/')
VAL_PATH <- "/disk/validationSet/coef/"
CODE_PATH <- "/root/biosoft/get_picture/"
# PGS_GROUP_LIST <- c(3: 5)
# PGS_SUB_LIST <- 3
PGS_GROUP_LIST <- c(3: 10)
PGS_SUB_LIST <- c(3: 6)

#############################################
##### Application 0: Estimate pheno_hat #####
#############################################
##### Get pheno_hat #####
## Load data
pheno_test <- fread2(paste0(TEST_PATH, "phenotype/", opt$sex, 
                            "/", opt$phenocode, ".txt"))[, 1, drop = T]
cov_test <- fread2(paste0(TEST_PATH, "phenotype/", opt$sex, 
                          "/cov.txt"))  %>%
              cbind(1, .) %>% as.matrix
cov_coef_val <- fread2(paste0(VAL_PATH, opt$sex, "/", opt$phenocode, ".txt")) %>% 
  as.matrix
p_type <- ifelse(grepl("Q", opt$phenocode), "gaussian", "binomial")
## Calculate pheno_hat
all_pred_file <- list.files(paste0(opt$out_path, "/pred_pheno"), 
                            pattern = "^pred_hm3.*\\.profile.gz$", 
                            full.names = T)
eid <- fread2(all_pred_file[1])[, 1]
PGS <- lapply(all_pred_file, function(all_pred_filex){
  
  predx <- fread2(all_pred_filex)
  if (all(!is.na(predx))){
    
    scorex <- predx[match(eid, predx$IID), "SCORESUM"]
  } else {
    
    scorex <- data.frame(SCORESUM = rep(0, nrow(predx)))[, 1]
  }
  return(scorex)
}) %>% Reduce("cbind", .) %>% rowSums()
p_type <- ifelse(grepl("Q", opt$phenocode), "gaussian", "binomial")
## quantitative traits
if(p_type == "gaussian"){
  
  pheno_hat <- PGS + cov_test %*% cov_coef_val
  PGS_test_res <- data.frame(pheno = pheno_test, 
                             pheno_hat = pheno_hat[, 1], 
                             PGS = PGS) %>%
    filter(!is.na(pheno))
  trait_label <- fread2("/disk/field_id/fieldID_quantitative.txt") %>% 
                  filter(., PFID == opt$phenocode) %>%
                  select(., "Reported Trait") %>% as.character()
  cat("R2:", round(cor(PGS_test_res$pheno, PGS_test_res$pheno_hat)^2, 4), "\n")
## binary traits
} else {
  
  pheno_hat <- 1 / (1 + exp(- PGS - cov_test %*% cov_coef_val))
  PGS_test_res <- data.frame(pheno = pheno_test, 
                             pheno_hat = pheno_hat[, 1], 
                             PGS = PGS) %>%
    filter(!is.na(pheno))
  trait_label <- fread2("/disk/field_id/fieldID_binary.txt") %>% 
    filter(., PFID == opt$phenocode) %>%
    select(., "Reported Trait") %>% as.character()
  cat("AUC:", round(auc(PGS_test_res$pheno, PGS_test_res$pheno_hat), 4), "\n")
}
saveRDS(PGS_test_res, paste0(opt$out_path, "/PGS_test_res.rds"))



######################################
##### Application 1: Performance #####
######################################
source(paste0(CODE_PATH, "viz_fun.R"))
## 1.1 Overall consistency
if (p_type == "gaussian") {

  perform_all_plt_annot <- dens.plt(pheno = PGS_test_res$pheno,
                                    pred = PGS_test_res$pheno_hat, 
                                    trait_label = trait_label, 
                                    watermark = T)
} else {

  perform_all_plt_annot <- roc.plt(pheno = PGS_test_res$pheno,
                                   pred = PGS_test_res$pheno_hat, 
                                   trait_label = trait_label, 
                                   watermark = T)
}
ggsave(file = paste0(opt$out_path, "/performance/Performance.png"), 
       perform_all_plt_annot, height = 6, width = 7,
       units = "in", dpi = 300)
# pdf(paste0(opt$out_path, "/performance/Performance.pdf"),
#     height = 6, width = 7)
#     # , units = "in", res = 300)
# print(perform_all_plt_annot)
# dev.off()

## 1.2 Gradient across grades of genetic risk
cov_test_na <- cov_test[!is.na(pheno_test), -1]
PGS_group <- lapply(PGS_GROUP_LIST, function(n){
  cut(PGS_test_res$PGS,
      breaks = quantile(PGS_test_res$PGS, probs = seq(0, 1, 1/n)),
      labels = seq(n))
}) %>% Reduce("cbind", .) %>% as.data.frame()
colnames(PGS_group) <- paste0("G", PGS_GROUP_LIST)
PGS_trend_plt <- lapply(PGS_GROUP_LIST, function(nx){

  ## format
  datt_plt <- data.frame(y = PGS_test_res$pheno,
                         PGS_g = PGS_group[[paste0("G", nx)]])
  if (!is.null(cov_test_na))
    datt_plt <- cbind(datt_plt, cov_test_na)
  ## plot
  PGS_trend_eff_pltx <- PGS.trend.plt(datt = datt_plt, n = nx, 
                                      type = p_type, watermark = TRUE) 
  ## save
  out_eff_filex <- paste0(opt$out_path, "/performance/Trend_PGS_group", nx, ".png")
  ggsave(file = out_eff_filex, PGS_trend_eff_pltx, height = 6, width = 7,
         units = "in", dpi = 300)
  message(paste0("MSG: Trend plot of PGS across ", nx, " groups is saved!"))
  return(nx)
})

## 1.3 plot {PGS} and performance (R2/AUC) by strata
strata_test <- readRDS(paste0(TEST_PATH, "phenotype/", opt$sex, "/strata_df.rds")) %>%
  filter(!is.na(pheno_test))
if (opt$sex != "All")
  strata_test$Sex <- NULL
strata_var <- colnames(strata_test)
strata_res <- plyr::alply(strata_var, 1, function(varx){

  datt_strata <- PGS_test_res
  datt_strata$Group <- strata_test[[varx]]
  datt_strata <- datt_strata[which(datt_strata$Group != "NA"), ]

  ## 3.3 plot PGS in subgroups
  PGS_strata_eff_pltx <- PGS.strata.plt(datt = datt_strata, 
                                        strata_var = varx, 
                                        watermark = TRUE)
  out_PGS_eff_filex <- paste0(opt$out_path, "/performance/PGS_subgroup_by_", varx, ".png")
  ggsave(file = out_PGS_eff_filex, PGS_strata_eff_pltx, height = 6, width = 6,
         units = "in", dpi = 300)
  message(paste0("MSG: PGS plot by ", varx, " is saved!"))
  ## 3.4 plot effect in subgroups
  datt_perform_strata <- perform.strata.test(pheno = datt_strata$pheno,
                                             pred = datt_strata$pheno_hat,
                                             subgroup = datt_strata$Group,
                                             type = p_type)
  perform_strata_eff_plt <- perform.strata.plt(datt = datt_perform_strata$perform_strata,
                                               test = datt_perform_strata$test,
                                               type = p_type,
                                               strata_var = varx,
                                               watermark = TRUE)
  plt_widx <- ifelse(length(levels(datt_strata$Group)) <= 3, 5, length(levels(datt_strata$Group)) + 2)
  out_perform_eff_filex <- paste0(opt$out_path, "/performance/Performance_subgroup_by_", varx, ".png")
  ggsave(file = out_perform_eff_filex, perform_strata_eff_plt,
         height = 5, width = plt_widx, units = "in", dpi = 300)
  message(paste0("MSG: Performance plot grouped by ", varx, " is saved!"))

  return(cbind(varx, datt_perform_strata$perform_strata))
}) %>% do.call("rbind", .)
# strata_out <- data.frame(Variable = strata_res$varx,
#                          Group = strata_res$Group,
#                          PGS_effect = round(strata_res$perform, 4),
#                          Low_CI = round(strata_res$lowCI, 4),
#                          High_CI = round(strata_res$highCI, 4))
# write.table(strata_out, file = paste0(opt$out_path, "/performance/Performance_subgroup.txt"),
#             row.names = F, quote = F, sep = "\t")



#######################################
##### Application 2: Joint effect #####
#######################################
## 2.1 PGS effect (Beta/OR) by strata (with interaction test)
interaction_res <- plyr::alply(strata_var, 1, function(varx){
 
   # cat(varx, "\n")
   ### 1 format base data frame
   datt_strata_plt <- PGS_test_res[, c("pheno", "PGS")]
 
   if (!is.null(cov_test_na)) {
 
     strata_cov_rmx <- ifelse(grepl("^Age", varx), "Age", varx)
     strata_covx <- setdiff(colnames(cov_test_na), strata_cov_rmx)
     datt_strata_plt <- cbind(datt_strata_plt,
                              cov_test_na[, strata_covx])
   }
   datt_strata_plt$Group <- strata_test[[varx]]
   ### 2 viz
   effect_strata_resx  <- effect.strata.plt(datt = datt_strata_plt,
                                            strata_var = varx,
                                            model = p_type,
                                            watermark = TRUE)
   ### 3 output
   plt_widx <- nlevels(datt_strata_plt$Group) * 2
   out_interaction_eff_filex <- paste0(opt$out_path, "/jointeffect/Interaction_with_", varx, ".png")
   ggsave(file = out_interaction_eff_filex, effect_strata_resx[[1]],
          height = 6, width = plt_widx, units = "in", dpi = 300)
   message(paste0("MSG: Interaction plot with ", varx, " is saved!"))
   return(effect_strata_resx[[2]])
}) %>% do.call("rbind", .)
# write.table(interaction_res, file = paste0(opt$out_path, "/jointeffect/Interaction_result.txt"),
#             row.names = F, quote = F, sep = "\t")

## 2.2 Subgroup analysis of PGS group and strata
sub_pgsg_plt <- lapply(strata_var, function(varx){

  ## 1 format base data frame
  datt_sub_plt <- data.frame(y = PGS_test_res$pheno,
                             Group = strata_test[[varx]])
  if (!is.null(cov_test_na)) {

    strata_cov_rmx <- ifelse(grepl("^Age", varx), "Age", varx)
    strata_covx <- setdiff(colnames(cov_test_na), strata_cov_rmx)
    datt_sub_plt <- cbind(datt_sub_plt,
                          cov_test_na[, strata_covx])

  }
  ## 2 plot for each PGS group setting
  # nx=5
  sub_pltx <- lapply(PGS_SUB_LIST, function(nx){

    # add PGS group var
    datt_sub_plt$PGS_g <- PGS_group[[paste0("G", nx)]] %>%
      paste0("PGS_group", .) %>%
      factor(., levels = paste0("PGS_group", 1:nx))

    # ## subgroup analysis
    datt_sub_pltx <- subset(datt_sub_plt, rowSums(is.na(datt_sub_plt)) == 0)
    datt_sub_pltx[[varx]] <- datt_sub_pltx$Group
    # ## sub Group by PGS_g
    # sub_plt1 <- subgroup.plt(datt = datt_sub_pltx,
    #                           type = p_type,
    #                           sub_var = varx,
    #                           eff_var = "PGS_g")
    # out_subgroup_prefix1 <- paste0(opt$out_path, "/jointeffect/Subgroup_of_PGS_g", nx, "_by_", varx)
    # plt_widx <- 4 + max(nlevels(datt_sub_plt$Group), nx) * 1.3
    # # png(paste0(out_subgroup_prefix1, ".png"),
    # #     height = 5, width = plt_widx, units = "in", res = 300)
    # pdf(paste0(out_subgroup_prefix1, ".pdf"),
    #     height = 5, width = plt_widx)
    # print(sub_plt1$plt)
    # dev.off()
    # write.table(sub_plt1$tab, file = paste0(out_subgroup_prefix1, ".txt"),
    #             sep = "\t", col.names = T, row.names = F, quote = F)
    # message(paste0("MSG: Subgroup analysis plot of PGS in group ", nx,
    #                " by ", varx, " is saved!"))
    # source(paste0(CODE_PATH, "viz_fun.R"))
    ## sub PGS_g by Group
    sub_eff_plt2 <- subgroup.plt(datt = datt_sub_pltx, type = p_type,
                                 sub_var = "PGS_g", eff_var = varx, 
                                 watermark = TRUE)
    # output
    plt_widx <- 4 + max(nlevels(datt_sub_plt$Group), nx) * 1.3
    out_subgroup_eff_prefix2 <- paste0(opt$out_path, "/jointeffect/Subgroup_of_", varx, "_by_PGS_g", nx)
    ggsave(paste0(out_subgroup_eff_prefix2, ".png"), 
                  sub_eff_plt2$plt, 
           height = 5, width = plt_widx, units = "in", dpi = 300)
    # write.table(sub_eff_plt2$tab, file = paste0(out_subgroup_eff_prefix2, ".txt"),
    #             sep = "\t", col.names = T, row.names = F, quote = F)
    message(paste0("MSG: Subgroup analysis plot of ", varx,
                   " by PGS in group ", nx, " is saved!"))
    return(nx)
  }) %>% unlist
  return(varx)
})
