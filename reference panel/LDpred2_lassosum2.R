######################
# Make LD shrinkage panel for LDpred2-grid, LDpred2-auto and lassosum2

# Load packages
library(bigsnpr)

# Set parameters
OMNI_PATH <- "/public/home/Datasets/interpolated_OMNI"
NCORES <- 8

DATA_PATH <- c("/public/home/Datasets/1000GP/AFR/hm3_imp/", 
               "/public/home/Datasets/1000GP/EAS/hm3_imp/", 
               "/public/home/biostat03/project/pgsfusionProject/ref_geno/")

for (dd in 1: 3){
  setwd(DATA_PATH[dd])
  
  # Load genotype data
  snp_readBed("./merge.bed")
  # obj.bigSNP <- snp_attach("./hm3_imp/mergeout.rds")
  obj.bigSNP <- snp_attach("./merge.rds")
  G <- snp_fastImputeSimple(obj.bigSNP$genotypes)
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  POS2 <- snp_asGeneticPos(CHR, POS, dir = OMNI_PATH, ncores = NCORES)
  
  # Estimate LD structure
  for (chr in 1:22) {
    
    print(chr)
    ind.chr <- which(CHR == chr)
    corr <- snp_cor(G, ind.col = ind.chr, size = 3 / 1000,
                    infos.pos = POS2[ind.chr], ncores = NCORES)
    corr <- as_SFBM(corr0,
                    paste0("./tmp-data1/LD_sparse_chr", chr),
                    compact = TRUE)
    saveRDS(corr, file = paste0("./LD_with_block_chr", chr, ".rds"))
  
    if (chr == 1) {
      
      ld <- Matrix::colSums(corr^2)
     
    } else {
      
      ld <- c(ld, Matrix::colSums(corr^2))
    }
  }
  
  # Construct SNP information data
  maf <- snp_MAF(G)
  map_hm3_plus <- data.frame(chr = CHR,
                             pos = POS,  
                             gd_pos = POS2,  
                             rsid = obj.bigSNP$map$marker.ID,
                             af_ref = maf,
                             a0 = obj.bigSNP$map$allele2,
                             a1 = obj.bigSNP$map$allele1,
                             ldsc = ld)
  saveRDS(map_hm3_plus, file = "map_hm3_plus.rds")
}
