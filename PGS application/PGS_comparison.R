# Load packages
library(boot)
library(plyr)
library(dplyr)
library(pROC)
library(ggplot2)
library(optparse)

# Input parameters
args_list = list(
  make_option("--job_path", type="character", default=NULL,
              help="INPUT: phenotype code", metavar="character")
)
opt_parser = OptionParser(option_list = args_list)
opt = parse_args(opt_parser)

# opt = list(job_path = "/home/chencao_pgs/website/pgsfusion-server/job/f4de6a6da78a4c09a727d3647de9e52b/")

# Function: correlation
corr.fun <- function(df, idx){
  
  df_sub <- df[idx, ]
  c(cor(df_sub[, 1], df_sub[, 2], method = "pearson")^2)
}

# Set global variables
BOOT_NUM <- 1000
method_match <- data.frame(file = c("CT", "DBSLMM_TUNING", "DBSLMM_AUTO", "DBSLMM_LMM", "LASSOSUM", 
                                    "LDPRED_NOSP", "LDPRED_AUTO", "MEGAPRS", "PRSCS_TUNING", "PRSCS_AUTO", 
                                    "SDPR", "MTPGS", "ANNOPRED", "SBAYESRC", "PRSCSX", "SDPRX", "XPASS"), 
                           method_name = c("CT", "DBSLMM", "DBSLMM-auto", "DBSLMM-lmm", "lassosum2", 
                                           "LDpred2", "LDpred2-auto", "MegaPRS-BayesR", "PRS-CS", "PRS-CS-auto", 
                                           "SDPR", "mtPGS", "AnnoPred", "SBayseRC", "PRS-CSx", "SDPRX", "XPASS"), 
                           color =  c("#800000", "#2A9D8F", "#46F0F0", "#AAFFC3", "#6A3569",
                                      "#4363D8", "#3D5A80", "#0077b6", "#FFE119", "#E9C46A",
                                      "#A4277C", "#FFD8B1", "#BCF60C", "#3CB44B", "#E6194B",
                                      "#F58231", "#911EB4"))

# Load data
method_use <- read.table(paste0(opt$job_path, "/method_use.txt"))[, 1, drop = T]
# method_use <- c("DBSLMM_LMM", "PRSCS_AUTO",
#                "SDPR", "PRSCSX", "SDPRX", "XPASS")
method_match_sub <- method_match[method_match[, 1] %in% method_use, ]
comp_df <- alply(method_match_sub[, 1], 1, function(m){
  
  test_df <- readRDS(paste0(opt$job_path, m, "/PGS_test_res.rds"))
  if (length(unique(test_df$pheno)) == 2){
    
    roc_fit <- roc(test_df$pheno, test_df$pheno_hat)
    ci_fit <- ci(roc_fit, of = "auc")
    return(c(ci_fit[2], ci_fit[1], ci_fit[3]))
  } else {
    
    set.seed(20250207)
    r2 <- cor(test_df[, 1], test_df[, 2], method = "pearson")^2
    boot_fit <- boot(test_df, corr.fun, R = BOOT_NUM)
    boot_ci <- boot.ci(boot.out = boot_fit, type = "basic")$basic[4:5]
    return(c(r2, boot_ci))
  }
}) %>% do.call('rbind' ,.) %>% as.data.frame()

colnames(comp_df) <- c("Index", "low", "up")
comp_df$method <- method_match_sub[, 2]
comp_df$method <- factor(method_match_sub[, 2],
levels = c("CT", "DBSLMM", "DBSLMM-auto", "DBSLMM-lmm", "lassosum2",
           "LDpred2", "LDpred2-auto", "MegaPRS-BayesR", "PRS-CS", "PRS-CS-auto",
           "SDPR", "mtPGS", "AnnoPred", "SBayseRC", "PRS-CSx", "SDPRX", "XPASS"))
# comp_df$method <- factor(method_match_sub[, 2],
#                         levels =  c("DBSLMM-lmm", "XPASS" ,"PRS-CS-auto", "PRS-CSx",
#                                      "SDPR", "SDPRX"))
# comp_df
# Plot 
plt <- ggplot(comp_df, aes(x = method, y = Index, fill = method)) + 
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = low, ymax = up), 
                size = 0.8, width = 0.2) +
  # geom_text(aes(label = round(Index, 4)), vjust = -0.5, 
  #           size = 4) +
  scale_fill_manual(values = method_match_sub$color) +
  xlab("Compared Methods") +
  theme_bw() +
  geom_text(aes(label = round(Index, 4)), vjust = -0.5, 
            size = 4) + 
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.title.x = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        panel.grid = element_blank(), 
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# y lab
test_df <- readRDS(paste0(opt$job_path, method_match_sub[1, 1], "/PGS_test_res.rds"))
if (length(unique(test_df$pheno)) == 2){
  
  plt1 <- plt + ylab("Prediction AUC")
} else {
  
  plt1 <- plt + ylab(expression("Prediction " ~ italic(R^2))) 
}

# Output
width_set <- ifelse(nrow(comp_df)>=4, 5.5, nrow(comp_df))
plt1 <- plt1 + annotate(geom="text",
                        x = length(comp_df$method)/2,
                        y = mean(c(0, max(comp_df$up))),
                        label = 'Created by PGSFusion\nUKBB application 144904', 
                        color = 'grey', angle = 45, 
                        size = 6, alpha = 0.5)
ggsave(paste0(opt$job_path, "/Performance_comparsion.png"), plt1, 
       width = width_set, height = 6)
# write.table(comp_df, file = paste0(opt$job_path, "/in_Performance_comparsion.txt"), 
#             quote = F, row.names = F)
