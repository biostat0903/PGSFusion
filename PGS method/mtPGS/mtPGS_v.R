#usr/bin/bash!

library(bigreadr)
library(GECKO)
library(optparse)
# Parameter setting
args_list <- list(
  make_option("--summ1", type = "character", default = NULL,
              help = "INPUT: Summary1 path", metavar = "character"),
  make_option("--summ2", type = "character", default = NULL,
              help = "INPUT: Summary2 path", metavar = "character"),
  make_option("--nsin", type = "numeric", default = NULL,
              help = "annotation", metavar = "character"),
  make_option("--output", type = "character", default = NULL,
              help = "INPUT: a output path", metavar = "character"),
  make_option("--anc", type = "character", default = NULL,
              help = "INPUT: ancestry", metavar = "character")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

#opt <- list(summ1='/root/case/3/bmi_giant.assoc.txt',
            #summ2='/root/case/3/ldl_eur.txt',
            #nsin=0,
            #output='/root/pgsfusion_test/job-mtPGS_file')

##################### calculate var and cov ##########################
stat1 <- fread2(opt$summ1)
stat2 <- fread2(opt$summ2)

n1 <- stat1[2,4] + stat1[2,5]
n2 <- stat2[2,4] + stat2[2,5]

z1 <- stat1[,9]/stat1[,10]
z2 <- stat2[,9]/stat2[,10]

stat1 <- data.frame(stat1[,c(1,3,2,6,7)],n1,z1,stat1[,11])
stat2 <- data.frame(stat2[,c(1,3,2,6,7)],n2,z2,stat2[,11])

colnames(stat1) <- c('chr','bp','SNP','A1','A2','N','Z','P')
colnames(stat2) <- c('chr','bp','SNP','A1','A2','N','Z','P')

n1in <- round(mean(stat1$N))
n2in <- round(mean(stat2$N))
nsin <- opt$nsin
Weightin <- T
Fix_Vein <- ifelse(nsin==0 , T , F)
Test <- T
###################### read ldscore file ##########################
ldscore <- fread2(paste0('/disk/reference_pgsfusion/',opt$anc,'_w_ld_chr/ldscore.txt'))

snp <- ldscore[which(ldscore[,2] %in% stat1[,3] & ldscore[,2] %in% stat2[,3]),2]

ldscore <- ldscore[which(ldscore[,2] %in% snp),]
stat2 <- stat2[which(stat2[,3] %in% snp),]
stat1 <- stat1[which(stat1[,3] %in% snp),]
Result <- GECKO_R(stat1,stat2,n1in,n2in,nsin,ldscore,Weightin,Fix_Vein,Test)

################## vg and ve #######################
v_e <- matrix(data=c(Result[1,1],Result[1,2],Result[1,2],Result[1,3]),
              nrow=2,ncol=2,byrow=T)
v_g <- matrix(data=c(Result[1,5],Result[1,6],Result[1,6],Result[1,7]),
              nrow=2,ncol=2,byrow=T)

write.table(v_e,file=paste0(opt$output,'/v_e.txt'),
            col.names=F,row.names=F,quote=F)
write.table(v_g,file=paste0(opt$output,'/v_g.txt'),
            col.names=F,row.names=F,quote=F)