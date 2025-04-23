# Load packages
library(bigreadr)
library(dplyr)
library(plyr)

# Set parameters
PATH <- "/public/home/biostat04/Project/19_PGS_fusion/"
TRAIT <- "/public/home/Datasets/ukb/pheno/03_trait/participant.csv.gz"

# Load data
anc_dat <- fread2(TRAIT, 
                  select = c("eid", 
                             "p31",             ## Sex
                             "p21000_i0",       ## Ethnic background | Instance 0
                             "p21000_i1",       ## Ethnic background | Instance 1
                             "p21000_i2",       ## Ethnic background | Instance 2
                             "p21000_i3",       ## Ethnic background | Instance 3
                             "p22001",          ## Genetic sex
                             "p22010",          ## Recommended genomic analysis exclusions
                             "p22018",          ## Genetic relatedness exclusions
                             "p22020",          ## Used in genetic principal components
                             "p22028"           ## Used in phasing Chromosomes 1-22
                  ))
fam_file <- fread2("/public/home/Datasets/ukb/geno/plink/chr1.fam")[, 1]

# Define samples
eth_EUR <- c("White", "British", 
             "Any other white background") 
eth_ESA <- c("Chinese") 
eth_AFR <- c("Caribbean", "African",
             "Any other Black background")
cnd_i_EUR <- anc_dat$p21000_i0 %in% eth_EUR | anc_dat$p21000_i1 %in% eth_EUR |
  anc_dat$p21000_i2 %in% eth_EUR | anc_dat$p21000_i3 %in% eth_EUR 
cnd_i_EAS <- anc_dat$p21000_i0 %in% eth_ESA | anc_dat$p21000_i1 %in% eth_ESA |
  anc_dat$p21000_i2 %in% eth_ESA | anc_dat$p21000_i3 %in% eth_ESA 
cnd_i_AFR <- anc_dat$p21000_i0 %in% eth_AFR | anc_dat$p21000_i1 %in% eth_AFR |
  anc_dat$p21000_i2 %in% eth_AFR | anc_dat$p21000_i3 %in% eth_AFR 
cnd_ii <- anc_dat$p22018 == ""
cnd_iii <- anc_dat$p31 == anc_dat$p22001
cnd_iv <- anc_dat$p22010 == ""
cnd_v <- anc_dat$p22028 == "Yes" 
cnd_vi <- anc_dat$p22020 == "Yes" 
cnd_vii <- anc_dat$eid %in% fam_file
cnd <- data.frame(EUR_cnd = cnd_i_EUR & cnd_ii & cnd_iii & cnd_iv & cnd_v & cnd_vi & cnd_vii,
                  EAS_cnd = cnd_i_EAS & cnd_ii & cnd_iii & cnd_iv & cnd_v & cnd_vi & cnd_vii,
                  AFR_cnd = cnd_i_AFR & cnd_ii & cnd_iii & cnd_iv & cnd_v & cnd_vi & cnd_vii)
eid_eur <- anc_dat$eid[cnd[, 1]]
eid_eas <- anc_dat$eid[cnd[, 2]]
eid_afr <- anc_dat$eid[cnd[, 3]]
eid_eur_eas <- intersect(eid_eur, eid_eas) ## eid: 3624367, instance 2: Chinese
## intersect(eid_eur, eid_afr) integer(0)
## intersect(eid_afr, eid_eas) integer(0)
eid_eas <- eid_eas[!eid_eas%in%eid_eur_eas]
cat("Samples Remaining:\n\tEUR: ", 
    length(eid_eur), "\n\tEAS: ",    ## 370524
    length(eid_eas), "\n\tAFR: ",    ## 1418
    length(eid_afr), "\n")           ## 6944

# Output
write.table(cbind(eid_eas, eid_eas), 
            file = paste0("/public/home/Datasets/ukb/pheno/01_eid/EAS.txt"), 
            row.names = F, col.names = F, quote = F)
write.table(cbind(eid_eur, eid_eur), 
            file = paste0("/public/home/Datasets/ukb/pheno/01_eid/EUR.txt"), 
            row.names = F, col.names = F, quote = F)
write.table(cbind(eid_afr, eid_afr), 
            file = paste0("/public/home/Datasets/ukb/pheno/01_eid/AFR.txt"), 
            row.names = F, col.names = F, quote = F)
