library(rhdf5)

arg <- commandArgs(TRUE)

DATA_DIR <- "/gpfs/chencao/zhenghuili/data/geno/pgsfusion/ldblk_1kg/ldblk_1kg_eur/"
REGION_DIR <- "/gpfs/chencao/zhenghuili/data/geno/pgsfusion/ldblk_1kg/ldblk_1kg_eur/region/"

region <- read.table(paste0(REGION_DIR, "chr", arg[1], ".bed"))
regions <- nrow(region)

OUT <- paste0("ldblk_1kg_chr", arg[1], ".hdf5")
h5createFile(OUT)
off <- 0
for (n in seq_len(regions)) {
    tryCatch({
        ldblk <- read.table(paste0(DATA_DIR, arg[1], "_", n, ".matrix"), header = FALSE, row.names = NULL)
        ldblk <- as.matrix(ldblk)
        colnames(ldblk) < NULL
        snps <- read.table(paste0(DATA_DIR, arg[1], "_", n, ".snplist"), header = FALSE, row.names = NULL)$V1
        snps <- matrix(snps, nrow = 1)
        snps <- array(snps)
        h5createGroup(OUT, paste0("blk_", n - off))
        h5write(ldblk, OUT, paste0("blk_", n - off, "/ldblk"))
        h5write(snps, OUT, paste0("blk_", n - off, "/snplist"))
        # print(off)
    }, warning = function(e) {
        off <<- off + 1
    })
}
h5closeAll()
