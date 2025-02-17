# Load packages
library(boot)
library(plyr)
library(dplyr)
library(ggplot2)
library(optparse)

# Input parameters
args_list = list(
  make_option("--job_path", type="character", default=NULL,
              help="INPUT: phenotype code", metavar="character")
)
opt_parser = OptionParser(option_list = args_list)
opt = parse_args(opt_parser)

# opt = list(job_path = "/home/chencao_pgs/website/pgsfusion-server/job/94e4a83aeb6f4b8498ff3a023e1b8fe5/")

# Function: correlation
corr.fun <- function(df, idx){
  
  df_sub <- df[idx, ]
  c(cor(df_sub[, 1], df_sub[, 2], method = "pearson")^2)
}

# Set global variables
BOOT_NUM <- 1000
method_match <- data.frame(file = c("CT", "DBSLMM_TUNING", "DBSLMM_AUTO", "DBSLMM_LMM", "LASSOSUM", 
                                    "LDPRED_NOSP", "LDPRED_AUTO", "MEGAPRS", "PRSCS_TUNING", "PRSCS_AUTO", 
                                    "SDPR", "MTPGS", "ANNOPRED", "SBayesRC", "PRSCSX", "SDPRX", "XPASS"), 
                           method_name = c("CT", "DBSLMM", "DBSLMM-auto", "DBSLMM-lmm", "lassosum2", 
                                           "LDpred2", "LDpred2-auto", "MegaPRS-BayesR", "PRS-CS", "PRS-CS-auto", 
                                           "SDPR", "mtPGS", "AnnoPred", "SbayseRC", "PRS-CSx", "SDPRX", "XPASS"), 
                           color =  c("#800000", "#2A9D8F", "#46F0F0", "#AAFFC3", "#6A3569",
                                      "#4363D8", "#3D5A80", "#0077b6", "#FFE119", "#E9C46A",
                                      "#A4277C", "#FFD8B1", "#BCF60C", "#3CB44B", "#E6194B",
                                      "#F58231", "#911EB4"))


# Load data
method_use <- read.table(paste0(opt$job_path, "/method_use.txt"))[, 1, drop = T]
method_match_sub <- method_match[method_match[, 1] %in% method_use, ]
comp_df <- alply(method_match_sub[, 1], 1, function(m){
  
  test_df <- readRDS(paste0(opt$job_path, m, "/PGS_test_res.rds"))
  set.seed(20250207)
  r2 <- cor(test_df[, 1], test_df[, 2], method = "pearson")^2
  boot_fit <- boot(test_df, corr.fun, R = BOOT_NUM)
  boot_ci <- boot.ci(boot.out = boot_fit, type = "basic")$basic[4:5]
  return(c(r2, boot_ci))
}) %>% do.call('rbind' ,.) %>% as.data.frame()
colnames(comp_df) <- c("R2", "low", "up")
comp_df$method <- method_match_sub[, 2]

# Plot 
plt <- ggplot(comp_df, aes(x = method, y = R2, fill = method)) + 
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = low, ymax = up), 
                size = 0.8, width = 0.2) +
  geom_text(aes(label = round(R2, 4)), vjust = -0.5, 
            size = 4) +
  scale_fill_manual(values = method_match_sub$color) +
  xlab("Compared Methods") +
  ylab(expression("Prediction " ~ italic(R^2))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        panel.grid = element_blank(), 
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# Output
write.table(comp_df, file = paste0(opt$job_path, "/Performance_comparsion.txt"), 
            quote = F, row.names = F)
ggsave(paste0(opt$job_path, "/Performance_comparsion.png"), plt, 
       width = nrow(comp_df)*1.8, height = 6)

