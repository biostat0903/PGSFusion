
# Load pacages
library(bigreadr)
library(plyr)
library(tidyverse)

JOB_PATH="/home/chencao_pgs/website/pgsfusion-server/job/f4de6a6da78a4c09a727d3647de9e52b"

# Function 1: transfer to data frame of memory size and time
toDF <- function(data = NULL, 
                 logs = T){
  
  time <- unlist(strsplit(data[5], " "))[8] %>% 
    strsplit(., ":") %>% unlist %>% as.numeric
  memory <- unlist(strsplit(data[10], " "))[6] %>% 
    as.numeric
  if (logs){
    
    time_min <- ifelse(length(time) == 2, time[1]+time[2]/60, 
                       time[1]*60+time[2]+time[3]/60) %>%
      log10 %>% round(., 3)
    memory_gb <- round(log10(memory/1024/1024), 3)
  } else {
    
    time_min <- ifelse(length(time) == 2, time[1]+time[2]/60, 
                       time[1]*60+time[2]+time[3]/60) %>%
      round(., 3)
    memory_gb <- round(memory/1024/1024, 3)
  }

  return(c(time_min, memory_gb))
}

method_match <- data.frame(file = c("CT", "DBSLMM_TUNING", "DBSLMM_AUTO", "DBSLMM_LMM", "LASSOSUM", 
                                    "LDPRED_NOSP", "LDPRED_AUTO", "MEGAPRS", "PRSCS_TUNING", "PRSCS_AUTO", 
                                    "SDPR", "MTPGS", "ANNOPRED", "SBayesRC", "PRSCSX", "SDPRX", "XPASS"), 
                           method_name = c("CT", "DBSLMM", "DBSLMM-auto", "DBSLMM-lmm", "lassosum2", 
                                           "LDpred2", "LDpred2-auto", "MegaPRS-BayesR", "PRS-CS", "PRS-CS-auto", 
                                           "SDPR", "mtPGS", "AnnoPred", "SbayseRC", "PRS-CSx", "SDPRX", "XPASS"))
method_use <- fread2(paste0(JOB_PATH, "/method_use.txt.tmp"), header = F)[, 1]
method_match_sub <- method_match[match(method_use, method_match$file), ]

summ_res <- aaply(method_use, 1, function(mm){
  
  res_m <- fread2(paste0(JOB_PATH, "/output/", mm, ".txt"))[, 2] %>% 
    toDF(., logs = T)
  return(res_m)
})
df <- cbind.data.frame(Methods = method_match_sub$method_name, 
                       Time = as.numeric(summ_res[, 1]), 
                       Memory = as.numeric(summ_res[, 2]))
df
sf <- max(df$Time)/max(df$Memory)
DF_long <- df %>%
  mutate(Memory = Memory*sf) %>%
  pivot_longer(names_to = "y_new", values_to = "val", Time:Memory)


#Plot
DF_long$Methods_fac <- factor(DF_long$Methods, 
                              # levels = c("DBSLMM-lmm", "XPASS", 
                              #            "PRS-CS-auto", "PRS-CSx",
                              #            "SDPR", "SDPRX"))
                              levels =  c("CT", "DBSLMM", "DBSLMM-auto", "DBSLMM-lmm", "lassosum2", 
                                          "LDpred2", "LDpred2-auto", "MegaPRS-BayesR", "PRS-CS", "PRS-CS-auto", 
                                          "SDPR", "AnnoPred", "SbayseRC"))
plt <- ggplot(DF_long, aes(x = Methods_fac)) +
  geom_bar(aes(y = val, fill = y_new, group = y_new),
           stat ="identity", position = position_dodge(),
           color = "grey", 
           alpha = 0.6, width = 0.9)  +
  scale_fill_manual(values = c('#88B6C3','#F0CF7B')) +
  scale_y_continuous(name = "log10(CPU time (min))", 
                     limits = c(0, 4),
                     breaks = seq(0, 4, 1), 
                     expand=c(0, 0.3),
                     sec.axis = sec_axis(~./1, name = "log10(Memory (GB))",
                                         breaks = seq(0, 4, 1)))+
  labs(fill='')+
  xlab("Compared Methods") +
  theme_bw()+
  theme(legend.position = c(0.18, 0.9),
        axis.text.x = element_text(color = 'black',size = 12,
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', face = 'bold', size = 15),
        legend.text = element_text(color = 'black', face = 'bold', size = 12),
        legend.title = element_text(color = 'black', face = 'bold', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggsave(paste0(JOB_PATH, '/memory_time.pdf'), plt, 
       width = 6, height = 6)

