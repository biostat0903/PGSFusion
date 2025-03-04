
library(bigreadr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

PATH <- "/home/chencao_pgs/website/pgsfusion-server/job/f4de6a6da78a4c09a727d3647de9e52b"

gwasResults <- fread2(paste0(PATH, "/trait1.txt"))
ref_EUR <- fread2("/disk/reference_pgsfusion/SNP_BP_GRCh37.txt.gz")

## calculate the length of each chromosome
chr_len <- ref_EUR %>% 
  group_by(CHR) %>% 
  summarise(chr_len = max(BP))
## set the start x of each chr
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 
## calculate the cumulative position of each SNP
SNP_info <- chr_pos %>%
  left_join(ref_EUR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + total)
# set x axis
X_axis <-  SNP_info %>% 
  group_by(CHR) %>% 
  summarize(center=(max(BPcum) + min(BPcum))/2)
SNP_info$PVAL <- gwasResults$p_wald[match(SNP_info$SNP, gwasResults$rs)] %>% as.numeric()
SNP_info$PVAL[SNP_info$PVAL == 0] <- min(SNP_info$PVAL[SNP_info$PVAL != 0])
SNP_info$CHR <- as.factor(SNP_info$CHR)
SNP_info_sub <- SNP_info[!is.na(SNP_info$PVAL), ]

## Manhattan plot
Manh_plt <- ggplot(SNP_info_sub) +
  geom_point(aes(x = BPcum, y = -log10(PVAL), color = CHR),
             alpha = 0.8, size = 0.8)+
  scale_x_continuous(label = X_axis$CHR, 
                     breaks= X_axis$center, 
                     limits = range(SNP_info$BPcum)) +
  scale_y_continuous(limits = c(0, 160), 
                     expand = c(0, 0)) +
  scale_color_manual(values = rep(c("#DE423D", "#3D5587"), 11)) +
  xlab("Chromosome") + ylab(bquote("-"~log[10]~"(P-value)"))+
  geom_hline(yintercept = -log10(5E-8), 
             color = 'red', linewidth = 0.8, linetype = 2) + 
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold", color = "black"),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

## QQ plot
qq_df <- data.frame(obs = -log10(sort(SNP_info_sub$PVAL, decreasing = FALSE)),
                    exp = -log10(ppoints(length(SNP_info_sub$PVAL))))
qq_plt <- ggplot(data = qq_df, aes(exp, obs))+
  geom_point(alpha = 0.8, color = "grey60")+
  geom_abline(color = "red", linewidth = 0.8, linetype = 2)+
  xlab(bquote("Expected -"~log[10]~"(P-value)"))+
  ylab(bquote("Observed -"~log[10]~"(P-value)"))+
  theme_bw()+
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank())

## output
tiff(paste0(PATH, "/gwas_summary.tiff"), 
     height = 6, width = 14, units = "in", 
     res = 300, compression = "lzw")
(Manh_plt | qq_plt) + 
  plot_layout(widths = c(10, 4))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()

CODE_PATH <- "/root/biosoft/get_picture/"
source(paste0(CODE_PATH, "viz_fun.R"))


trait_label = "Dementia"


roc.plt <- function(pheno = NULL,
                    pred = NULL,
                    trait_label = NULL, 
                    type = "b"
){
  
  # AUC
  perform_test <- perform.test(pheno = pheno, 
                               pred = pred, 
                               type = type)
  auc_vec <- perform_test[["perform"]] %>% round(3)
  auc <- paste0("AUC of ",
                trait_label, ": ",
                auc_vec[1] , "(",
                auc_vec[2], "~", 
                auc_vec[3], ")")
  roc <- perform_test[["test"]]
  datt <- data.frame(Model = "ROC",
                     rev_specificity = 1 - roc$specificities,
                     sensitivity = roc$sensitivities)
  # plot
  plt <- ggplot(datt) +
    geom_line(aes(x = rev_specificity, y = sensitivity), 
              color = "#4DAF4A", linetype = 2, size = 1.3) +
    geom_abline(slope = 1, intercept = 0, 
                color = "grey10", linetype = 2) + 
    annotate("text", 
             x = 0.7, y = 0.2, 
             size = 5, colour="grey10",
             label = auc) + 
    xlab("1-Specificity") + ylab("Sensitivity") +
    theme_bw() + 
    theme(axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 12, color = "black"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"))
  
  return(plt)
}


perform_plt <- plyr::alply(c("DBSLMM_AUTO", "MEGAPRS"),  1, function(method){
                        
  PGS_test_res <- readRDS(paste0(PATH, "/", method, "/PGS_test_res.rds"))
  perform_plt_m <- roc.plt(pheno = PGS_test_res$pheno,
                           pred = PGS_test_res$pheno_hat, 
                           trait_label = trait_label)
  return(perform_plt_m)
})

tiff(paste0(PATH, "/ROC.tiff"), 
     height = 6, width = 14, units = "in", 
     res = 300, compression = "lzw")
(perform_plt[[1]] | perform_plt[[2]])  +
  plot_layout(widths = c(7, 7))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()


TEST_PATH <- "/disk/testSet/EUR/"
PGS_test_res <- readRDS(paste0(PATH, "/", "DBSLMM_AUTO", "/PGS_test_res.rds"))
pheno_test <- fread2(paste0(TEST_PATH, "phenotype/All/PFIDB2006.txt"))[, 1, drop = T]
cov_test <- fread2(paste0(TEST_PATH, "phenotype/All/cov.txt"))  %>%
  cbind(1, .) %>% as.matrix
cov_test_na <- cov_test[!is.na(pheno_test), -1]
p_type = "binomial"
PGS_GROUP_LIST = c(3, 4, 5, 7, 8, 9)
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
  PGS_trend_pltx <- PGS.trend.plt(datt = datt_plt,
                                  n = nx,
                                  type = p_type)

  
  message(paste0("MSG: Trend plot of PGS across ", nx, " groups is saved!"))
  return(PGS_trend_pltx)
})

tiff(paste0(PATH, "/group.tiff"), 
     height = 16, width = 14, units = "in", 
     res = 300, compression = "lzw")
((PGS_trend_plt[[1]] | PGS_trend_plt[[2]]) /
    (PGS_trend_plt[[3]] | PGS_trend_plt[[4]]) /
    (PGS_trend_plt[[5]] | PGS_trend_plt[[6]]) ) +
  plot_layout(widths = c(7, 7))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()


## 2.2 Subgroup analysis of PGS group and strata
varx = "SES"
strata_test <- readRDS(paste0(TEST_PATH, "phenotype/All/strata_df.rds")) %>%
  filter(!is.na(pheno_test))
PGS_SUB_LIST <- c(3, 4, 5)
  
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
  sub_pltx <- lapply(PGS_SUB_LIST, function(nx){
    
    # add PGS group var
    datt_sub_plt$PGS_g <- PGS_group[[paste0("G", nx)]] %>%
      paste0("PGS_group", .) %>%
      factor(., levels = paste0("PGS_group", 1:nx))
    
    # ## subgroup analysis
    datt_sub_pltx <- subset(datt_sub_plt, rowSums(is.na(datt_sub_plt)) == 0)
    datt_sub_pltx[[varx]] <- datt_sub_pltx$Group
    ## sub PGS_g by Group
    sub_plt2 <- subgroup.plt(datt = datt_sub_pltx,
                             type = p_type,
                             sub_var = "PGS_g",
                             eff_var = varx)
    # # save
    # out_subgroup_prefix2 <- paste0(opt$out_path, "/jointeffect/Subgroup_of_", varx, "_by_PGS_g", nx)
    # plt_widx <- 4 + max(nlevels(datt_sub_plt$Group), nx) * 1.3
    # png(paste0(out_subgroup_prefix2, ".png"),
    #     height = 5, width = plt_widx, units = "in", res = 300)
    # # pdf(paste0(out_subgroup_prefix2, ".pdf"),
    # #     height = 5, width = plt_widx)
    # print(sub_plt2$plt)
    # dev.off()

    message(paste0("MSG: Subgroup analysis plot of ", varx,
                   " by PGS in group ", nx, " is saved!"))
    
    return(sub_plt2)
  })

  tiff(paste0(PATH, "/group_SES.tiff"), 
       height = 16, width = 12, units = "in", 
       res = 300, compression = "lzw")
  (sub_pltx[[1]]$plt /
   sub_pltx[[2]]$plt /
   sub_pltx[[3]]$plt) +
    # plot_layout(widths = c(7, 7))+
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 15, face = "bold"))
  dev.off()


