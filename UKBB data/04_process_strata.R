# Load packages
library(bigreadr)
library(dplyr)
library(stringr)
library(poLCA)

# Set parameters
PATH <- "/public/home/biostat04/Project/19_PGS_fusion/"
TRAIT <- "/public/home/Datasets/ukb/pheno/03_trait/participant.csv.gz"

# Load sample id and field id
sample_id <- list.files(paste0(PATH, "sample_id/")) %>% 
  alply(., 1, function(ff) fread2(paste0(PATH, "sample_id/", ff))[, 1])
names(sample_id) <- list.files(paste0(PATH, "sample_id/")) %>% 
  gsub(".fam", "", .)

# Load strata variables
habit_field <- c("p21003_i0",        ## 01 age when attended assessment centre
                 "p54_i0",           ## 02 UK Biobank assessment centre
                 "p21001_i0",        ## 03 BMI
                 "p20116_i0",        ## 04 Smoking status
                 "p2897_i0",         ## 05 Age stopped smoking
                 "p884_i0",          ## 06 Number of days/week of moderate physical activity 10+ minutes
                 "p894_i0",          ## 07 Duration of moderate activity
                 "p904_i0",          ## 08 Number of days/week of vigorous physical activity 10+ minutes
                 "p914_i0",          ## 09 Duration of vigorous activity
                 "p1309_i0",         ## 10 Fresh fruit intake
                 "p1319_i0",         ## 11 Dried fruit intake
                 "p1289_i0",         ## 12 Cooked vegetable intake
                 "p1299_i0",         ## 13 Salad / raw vegetable intake
                 "p1329_i0",         ## 14 Oily fish intake
                 "p1339_i0",         ## 15 Non-oily fish intake
                 "p1349_i0",         ## 16 Processed meat intake
                 "p1369_i0",         ## 17 Beef intake
                 "p1379_i0",         ## 18 Lamb/mutton intake
                 "p1389_i0",         ## 19 Pork intake
                 "p1438_i0",         ## 20 Bread intake
                 "p1448_i0",         ## 21 Bread type
                 "p1458_i0",         ## 22 Cereal intake
                 "p1468_i0",         ## 23 Cereal type
                 "p20117_i0",        ## 24 Alcohol 
                 "p1180_i0",         ## 25 Morning/evening person (chronotype)
                 "p1160_i0",         ## 26 Sleep duration
                 "p1200_i0",         ## 27 Sleeplessness/insomnia
                 "p1210_i0",         ## 28 Snoring
                 "p1220_i0",         ## 29 Daytime dozing/sleeping (narcolepsy)
                 "p20453"            ## 30 Ever taken cannabis
)
baseline_field <- c("p22189",                ## 31 Townsend deprivation index
                    "p738_i0",               ## 32 Income
                    paste0("p6138_i", 0:3),  ## 33-36 Education
                    paste0("p6142_i", 0:3)   ## 37-40 Employment status
)
habit_baseline_df <- fread2(TRAIT, select = c("eid", habit_field, baseline_field))

############## Construct SES ############## 
# Income, Education and Employment groups #
###########################################
## Categorize income groups
incm_cnd <- recode(habit_baseline_df[["p738_i0"]],
                  "Less than 18,000" = "1",
                  "18,000 to 30,999" = "2",
                  "31,000 to 51,999" = "3",
                  "52,000 to 100,000" = "4",
                  "Greater than 100,000" = "5") %>% as.integer()
## Categorize education groups (p6138_i0 ~ p6138_i3)
### set 6138_i0 as the highest education status
### 1-2: College or above
### 3-6: High school or equivalent
### -7: Less than high school (33853828)
edu_list <- str_split_i(habit_baseline_df[["p6138_i0"]], "\\|", 1)
edu_cnd <- recode(edu_list,
                  "College or University degree" = "3",
                  "A levels/AS levels or equivalent" = "3",
                  "O levels/GCSEs or equivalent" = "2",
                  "CSEs or equivalent" = "2",
                  "NVQ or HND or HNC or equivalent" = "2",
                  "Other professional qualifications eg: nursing, teaching" = "2",
                  "None of the above" = "1") %>% as.integer()
## Categorize employment groups (p6142_i0 ~ p6142_i3)
### Employed Status including those in paid employment or self-employed, retired, 
### doing unpaid or voluntary work, or being full or part time students)
emp_list <- str_split_i(habit_baseline_df[["p6142_i0"]], "\\|", 1)
emp_cnd <- recode(emp_list,
                   "In paid employment or self-employed" = "2",
                   "Retired" = "2",
                   "Doing unpaid or voluntary work" = "2",
                   "Full or part-time student" = "2",
                   "Looking after home and/or family" = "1",
                   "Unable to work because of sickness or disability" = "1",
                   "None of the above" = "1") %>% as.integer()
ses_df <- data.frame(Income = incm_cnd, 
                     Education = edu_cnd, 
                     Employment = emp_cnd,
                     row.names = habit_baseline_df$eid)
## Estimate SES
SES <- rep(NA, nrow(ses_df))
nna_idx <- which(rowSums(is.na(ses_df)) == 0)
set.seed(20250207)
SES_LCA3 <- poLCA(cbind(Income, Education, Employment) ~ 1, 
                  data = ses_df[nna_idx,],
                  nclass = 3, 
                  maxiter = 10000, 
                  graphs = F, 
                  tol = 1e-6, 
                  na.rm = T, 
                  probs.start = NULL, 
                  nrep = 1, 
                  verbose = T, 
                  calc.se = T)
SES[nna_idx] <- SES_LCA3$predclass
high_idx <- which.max(SES_LCA3$probs$Income[, 5] + SES_LCA3$probs$Income[, 4])
low_idx <- which.max(SES_LCA3$probs$Income[, 1])
medium_idx <- setdiff(1:3, c(high_idx, low_idx))
SES <- factor(SES, 
              levels = c(low_idx, medium_idx, high_idx), 
              labels = c("Low", "Medium", "High"))
#   high    Low Medium 
# 83383 114037 219188 


################## Construct lifestyle ################## 
# Nosmoke, Activity, Diet, Noalcohol, Sleep, Nocannabis #
######################################################### 
## Define no current smoke
nosmoke_cnd <- ifelse(habit_baseline_df[["p20116_i0"]] == "Never", 1,
                      ifelse(habit_baseline_df[["p20116_i0"]] %in% c("Previous", "Current"), 0, NA))
### Former smokers who had quit smoking more than 30 years
nosmoke_year <- as.numeric(habit_baseline_df[["p21003_i0"]]) - 
  as.numeric(habit_baseline_df[["p2897_i0"]])
nosmoke_cnd[which(nosmoke_year > 30)] <- 1

## Define activity
### frequency: moderate at least 5 days a week and vigorous once a week
num_moderate <- as.integer(habit_baseline_df[["p884_i0"]])
num_vigorous <- as.integer(habit_baseline_df[["p904_i0"]])
num_activity <- num_moderate > 5 & num_vigorous > 1
### 150 minutes moderate activity per week
time_moderate <- ifelse(as.numeric(habit_baseline_df[["p894_i0"]]) > 150, T, F)
### 75 minutes vigorous activity per week 
time_vigorous <- ifelse(as.numeric(habit_baseline_df[["p914_i0"]]) > 75, T, F)
### activity score
act_mat <- cbind(num_activity, time_moderate, time_vigorous)
activity_cnd <- rep(0, nrow(act_mat))
activity_cnd[rowSums(act_mat, na.rm = T) > 0] <- 1
activity_cnd[rowSums(is.na(act_mat)) == 3] <- NA

## Define healthy diet
### Fruits: per serving = fresh fruit- 1 piece OR dried fruit- 5 pieces
fruit_mat <- cbind(as.numeric(habit_baseline_df[["p1309_i0"]]), 
                   as.numeric(habit_baseline_df[["p1319_i0"]])/5) 
fruit_cnd <- rowSums(fruit_mat) >= 4 # 4 servings/day
fruit_cnd[rowSums(fruit_mat, na.rm = T) >= 4] <- T
### Vegetable: per serving = cooked/raw vegetables- 3 heaped tablespoons
veg_mat <- cbind(as.numeric(habit_baseline_df[["p1289_i0"]])/3, 
                 as.numeric(habit_baseline_df[["p1299_i0"]])/3)
veg_cnd <- rowSums(veg_mat) >= 4 # 4 servings/day
veg_cnd[rowSums(veg_mat,na.rm = T) >= 4] <- T
### Fish: 2 habit_baseline_df/week
fish_mat <- lapply(c("p1329_i0", "p1339_i0"), function(x){
  recode(habit_baseline_df[[x]],
         "Never" = "0", 
         "Less than once a week" = "1", 
         "Once a week" = "2",
         "2-4 times a week" = "3", 
         "5-6 times a week" = "4", 
         "Once or more daily" = "5") %>% as.integer()
}) %>% Reduce("cbind", .)
fish_cnd <- fish_mat[,1] >= 3 | fish_mat[,2] >= 3 # Oily fish intake 2 times/week OR Non-oily fish intake 2 times/week
fish_cnd[fish_mat[,1] ==2 & fish_mat[,2] == 2] <- T # include 2 + 2 (Once a week)
### processed red meat: 1 times/week
pmeat_cnd <- ifelse(habit_baseline_df[["p1349_i0"]] %in% 
                      c("Never", "Less than once a week", "Once a week"), T,
                    ifelse(habit_baseline_df[["p1349_i0"]] %in% 
                             c("2-4 times a week", "5-6 times a week", "Once or more daily"), F, NA))
### non-processed red meat: 1.5 times/week
npmeat_mat <- lapply(c("p1369_i0", "p1379_i0", "p1389_i0"), function(x){
  recode(habit_baseline_df[[x]],
         "Never" = "0", 
         "Less than once a week" = "1", 
         "Once a week" = "2",
         "2-4 times a week" = "3", 
         "5-6 times a week" = "4", 
         "Once or more daily" = "5") %>% as.integer()
}) %>% Reduce("cbind", .)
npmeat_cnd <- rowSums(meat_mat, na.rm = T) <= 3  # 1 + 2 OR 1 * 3
npmeat_cnd[npmeat_mat[,1] >=3 | npmeat_mat[,1] >=3 | npmeat_mat[,1] >=3] <- F # exclude 3 + 0 + 0 (2 times/week)
### Whole grains: 3 servings/day
grains_mat <- matrix(0, nrow(habit_baseline_df), 2)
bread_idx <- which(habit_baseline_df[["p1448_i0"]] == "Wholemeal or wholegrain") # wholemeal/wholegrain bread- 1 slice/day
grains_mat[bread_idx, 1] <- habit_baseline_df[["p1438_i0"]][bread_idx] %>% as.numeric()
cereal_idx <- which(habit_baseline_df[["p1468_i0"]] %in% c("Bran cereal (e.g. All Bran, Branflakes)", 
                                                           "Oat cereal (e.g. Ready Brek, porridge)", 
                                                           "Muesli")) # bran/oat/muesli cereal- 1 bowl/day 
grains_mat[cereal_idx, 2] <- habit_baseline_df[["p1458_i0"]][cereal_idx] %>% as.numeric()
grains_cnd <- rowSums(grains_mat, na.rm = T) >= 3
### Healthy diet: At least 4 of the 6 food groups
diet_mat <- data.frame(fruit_cnd, veg_cnd, fish_cnd, 
                       pmeat_cnd, npmeat_cnd, grains_cnd)
diet_cnd <- rowSums(diet_mat) >= 4
diet_cnd[rowSums(diet_mat, na.rm = T) >= 4] <- T # Can satisfy even when missing = F
diet_cnd[rowSums(is.na(diet_mat)) + rowSums(diet_mat,na.rm = T) < 4] <- F # Never satisfy even when missing = T

## Define never alcohol use
noalcohol_cnd <- ifelse(habit_baseline_df[["p20117_i0"]] == "Never", 1,
                      ifelse(habit_baseline_df[["p20117_i0"]] %in% c("Previous", "Current"), 0, NA))

## Define sleep behavior
# "1180_i0",         ## Morning/evening person (chronotype)
# "1160_i0",         ## Sleep duration
# "1200_i0",         ## Sleeplessness/insomnia
# "1210_i0",         ## Snoring
# "1220_i0",         ## Daytime dozing/sleeping (narcolepsy)
## chronotype (morning preference)
# 1	Definitely a 'morning' person
# 2	More a 'morning' than 'evening' person
chrono_cnd <- ifelse(habit_baseline_df[["p1180_i0"]] %in% 
                       c("Definitely a 'morning' person", "More a 'morning' than 'evening' person"), 1,
                     ifelse(habit_baseline_df[["p1180_i0"]] %in% 
                              c("Definitely an 'evening' person", "More an 'evening' than a 'morning' person"), 0, NA))
## Sleep duration 7-8h
duration_cnd <- ifelse(is.na(habit_baseline_df[["p1160_i0"]] %>% as.numeric()), NA,
                       ifelse(habit_baseline_df[["p1160_i0"]] %>% as.numeric() <= 8 &
                                habit_baseline_df[["p1160_i0"]] %>% as.numeric() >= 7, 1, 0) )
## 1: never or rarely insomnia
insomnia_cnd <- ifelse(habit_baseline_df[["p1200_i0"]] %in% c("Never/rarely"), 1,
                     ifelse(habit_baseline_df[["p1200_i0"]] %in% c("Sometimes", "Usually"), 0, NA))

## 2: No complaints of snoring
snoring_cnd <- ifelse(habit_baseline_df[["p1210_i0"]] == "No", 1,
                       ifelse(habit_baseline_df[["p1210_i0"]] == "Yes", 0, NA))
# no frequently narcolepsy
# 0	Never/rarely
# 1	Sometimes
narcolepsy_cnd <- ifelse(habit_baseline_df[["p1220_i0"]] %in% c("Never/rarely", "Sometimes"), 1,
                       ifelse(habit_baseline_df[["p1220_i0"]] %in% c("Often", "All of the time"), 0, NA))
# sleep mat
sleep_mat <- data.frame(chrono_cnd, duration_cnd, insomnia_cnd,
                        snoring_cnd, narcolepsy_cnd)
sleep_cnd <- rowSums(sleep_mat) >= 4
# Can satisfy even when missing = F
sleep_cnd[rowSums(sleep_mat,na.rm = T) >= 4] <- T
# Never satisfy even when missing = T
sleep_cnd[rowSums(is.na(sleep_mat)) + rowSums(sleep_mat,na.rm = T) < 4] <- F

## Define never use cannabis
cannabis_cnd <- ifelse(habit_baseline_df[["p20453"]] %in% c("No"), 1,
                       ifelse(grepl("Yes", habit_baseline_df[["p20453"]]), 0, NA))

## Build life factors data
lf_df <- data.frame(Nosmoke = nosmoke_cnd, 
                    Activity = activity_cnd, 
                    Diet = diet_cnd, 
                    Noalcohol = noalcohol_cnd,
                    Sleep = sleep_cnd,
                    Nocannabis = cannabis_cnd)

lifescore1 <- rep(NA, nrow(lf_df))
lifescore1[rowSums(lf_df, na.rm = T) >= 4] <- 3 
lifescore1[rowSums(lf_df[,c(1:5)], na.rm = T) == 3] <- 3
lifescore1[rowSums(lf_df, na.rm = T) >= 2 &
             (rowSums(lf_df, na.rm = T) + rowSums(is.na(lf_df))) < 4 ] <- 2 
lifescore1[(rowSums(lf_df, na.rm = T) + rowSums(is.na(lf_df))) < 2 | 
             rowSums(lf_df, na.rm = T) == 0] <- 1 

lifescore2 <- rep(2, nrow(lf_df))
lifescore2[rowSums(lf_df, na.rm = T) >= 4] <- 3 
lifescore2[rowSums(lf_df[,c(1:5)], na.rm = T) >= 3] <- 3 
lifescore2[rowSums(lf_df == 0, na.rm = T) >= 4] <- 1

lifescore3 <- rep(2, nrow(lf_df))
lifescore3[rowSums(lf_df, na.rm = T) >= 4] <- 3 
lifescore3[rowSums(lf_df[,c(1:5)], na.rm = T) >= 3] <- 3 
lifescore3[rowSums(lf_df == 0, na.rm = T) >= 5] <- 1

# Summarize strata variables
Age_60 <- ifelse(habit_baseline_df$p21022 <= 60, "Age <= 60", "Age > 60") %>%
  factor(., levels = c("Age <= 60", "Age > 60"))
Age_65 <- ifelse(habit_baseline_df$p21022 <= 65, "Age <= 65", "Age > 65") %>%
  factor(., levels = c("Age <= 65", "Age > 65"))
age_range <- range(habit_baseline_df$p21022)
label_Age_50_60 <- c("Age <= 50", "50 < Age <= 60", "Age > 60")
Age_50_60 <- cut(habit_baseline_df$p21022,
                 breaks = c(age_range[1], 50, 60, age_range[2]),
                 labels = label_Age_50_60,
                 ordered_result = T)
label_Age_45_55_65 <- c("Age <= 45", "45 < Age <= 55", "55 < Age <= 65", "Age > 65")
Age_45_55_65 <- cut(habit_baseline_df$p21022,
                    breaks = c(age_range[1], 45, 55, 65, age_range[2]),
                    labels = label_Age_45_55_65,
                    ordered_result = T)
strata_df <- data.frame(eid = habit_baseline_df[["eid"]],
                        Townsend = habit_baseline_df[["p22189"]],
                        Income = incm_cnd, 
                        Education = edu_cnd, 
                        Employment = emp_cnd,
                        SES = SES,
                        Nosmoke = nosmoke_cnd, 
                        Activity = activity_cnd, 
                        Diet = diet_cnd, 
                        Noalcohol = noalcohol_cnd,
                        Sleep = sleep_cnd,
                        Nocannabis = cannabis_cnd,
                        Age_60 = Age_60, 
                        Age_65 = Age_65, 
                        Age_50_60 = Age_50_60, 
                        Age_45_55_65 = Age_45_55_65, 
                        Lifescore1 = lifescore1,
                        Lifescore2 = lifescore2,
                        Lifescore3 = lifescore2)
strata_data <- alply(c(1: length(sample_id)), 1, function (ss){

  if (grepl("test", names(sample_id)[ss])){

    strata_var_s <- strata_df[match(sample_id[[ss]], strata_df$eid), ]
    saveRDS(strata_var_s[, -1], file = paste0(PATH, "out_pheno/",
                                              names(sample_id)[ss], "/strata_df.rds"))
  }
  return(ss)
})
