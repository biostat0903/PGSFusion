
library(bigreadr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

PATH <- "/home/chencao_pgs/website/pgsfusion-server/job/d449824760b74d55a1cf29aa02fe7d71"

# reference
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

# GWAS1
gwasResults <- fread2(paste0(PATH, "/trait1.txt"))
SNP_info$PVAL <- gwasResults$p_wald[match(SNP_info$SNP, gwasResults$rs)] %>% as.numeric()
SNP_info$PVAL[SNP_info$PVAL == 0] <- min(SNP_info$PVAL[SNP_info$PVAL != 0])
SNP_info$CHR <- as.factor(SNP_info$CHR)
SNP_info_sub <- SNP_info[!is.na(SNP_info$PVAL), ]

## Manhattan plot1
Manh_plt1 <- ggplot(SNP_info_sub) +
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

## QQ plot1
qq_df1 <- data.frame(obs = -log10(sort(SNP_info_sub$PVAL, decreasing = FALSE)),
                    exp = -log10(ppoints(length(SNP_info_sub$PVAL))))
qq_plt1 <- ggplot(data = qq_df1, aes(exp, obs))+
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


# GWAS2
## calculate the cumulative position of each SNP
SNP_info <- chr_pos %>%
  left_join(ref_EUR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + total)
# set x axis
X_axis <-  SNP_info %>% 
  group_by(CHR) %>% 
  summarize(center=(max(BPcum) + min(BPcum))/2)

gwasResults <- fread2(paste0(PATH, "/trait2.txt"))
SNP_info$PVAL <- gwasResults$p_wald[match(SNP_info$SNP, gwasResults$rs)] %>% as.numeric()
SNP_info$PVAL[SNP_info$PVAL == 0] <- min(SNP_info$PVAL[SNP_info$PVAL != 0])
SNP_info$CHR <- as.factor(SNP_info$CHR)
SNP_info_sub <- SNP_info[!is.na(SNP_info$PVAL), ]

## Manhattan plot2
Manh_plt2 <- ggplot(SNP_info_sub) +
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

## QQ plot2
qq_df2 <- data.frame(obs = -log10(sort(SNP_info_sub$PVAL, decreasing = FALSE)),
                     exp = -log10(ppoints(length(SNP_info_sub$PVAL))))
qq_plt2 <- ggplot(data = qq_df2, aes(exp, obs))+
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



## output1
tiff(paste0(PATH, "/gwas_summary.tiff"), 
     height = 12, width = 14, units = "in", 
     res = 300, compression = "lzw")
# ((Manh_plt1 | qq_plt1) / 
#   (Manh_plt2 | qq_plt2)) + 
Manh_plt1 + qq_plt1 + Manh_plt2 + qq_plt2 +
  plot_layout(widths = c(10, 4))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()


















CODE_PATH <- "/root/biosoft/get_picture/"
source(paste0(CODE_PATH, "viz_fun.R"))


trait_label = "Standing height"


perform_plt <- plyr::alply(c("DBSLMM_LMM", "XPASS", "PRSCS_AUTO", "PRSCSX", 
                             "SDPR", "SDPRX"),  1, function(method){
                        
  PGS_test_res <- readRDS(paste0(PATH, "/", method, "/PGS_test_res.rds"))
  perform_plt_m <- dens.plt(pheno = PGS_test_res$pheno,
                           pred = PGS_test_res$pheno_hat, 
                           trait_label = trait_label)
  return(perform_plt_m)
})

tiff(paste0(PATH, "/R2.tiff"), 
     height = 18, width = 14, units = "in", 
     res = 300, compression = "lzw")
  perform_plt[[1]] + perform_plt[[2]] +
  perform_plt[[3]] + perform_plt[[4]] +
  perform_plt[[5]] + perform_plt[[6]] +
  plot_layout(widths = c(7, 7))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()

