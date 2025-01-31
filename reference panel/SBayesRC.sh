
GCTB=/root/biosoft/gctb_2.05beta_Linux/gctb


ma_file="MA_file"                # GWAS summary data in COJO format
genotype="YOUR_GENOTYPE_Prefix"  # genotype prefix as LD reference (PLINK format), with {CHR} to spefify multiple
outDir="YOUR_OUTPUT"             # Output folder that would be created automatically
threads=4                        # Number of CPU cores for eigen decomposition
#---usually don't need change bellow
genoCHR="1"                       # If more than 1 genotype file, input range (e.g. "1-22") here.
refblock=""                      # Text file to define LD blocks, by default to use our GRCH37 coordination 
tool="gctb"                      # Command line to run gctb for generating the full LD matrix


##############################################
# Code
# Step1: generate the LD block information and script
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist
Rscript -e "SBayesRC::LDstep1(mafile='$ma_file', genoPrefix='$genotype', \
            outDir='$outDir', genoCHR='$genoCHR', blockRef='$refblock', log2file=TRUE)"