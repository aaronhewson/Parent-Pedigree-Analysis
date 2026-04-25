# Parent pedigree analysis of unidentified cultivars in the Jim Dunckley & Plant and Food Research Orchards, using Mendelian inconsistent error counts.
# Part of Aaron Hewson's MSc research.
# Methodology and code adapted from Muranty et al. 2020, https://doi.org/10.1186/s12870-019-2171-6. 

# Script requires data to be in PLINK binary (bed) format, and to use PLINK to estimate IBD.
# PLINK version used: PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)

# Prior to using this script, ensure SNP data is in PLINK .ped/.map format. Remove all duplicates and all triploids.

# Format data with PLINK --------------------------------------------------
# In IBD calculation, a minimum value X is chosen to reduce output size. Parent-offspring pair should have IBD = ~0.5.
# Cutoff of 0.4 is used to provide some leniency. This should capture all true parent-offspring pairs.

# Set working directory - containing data files and PLINK application
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Parent Pedigree Analysis/Inputs/Dip")

# Changing format from .ped to .bed
system("plink --file Dip_No_Dupe --make-bed --out Dip_Data")

# Run IBD calculation to generate .genome file, with a minimum IBD threshold
system("plink --bfile Dip_Data --genome full --min 0.4")

# Write IBD results to a .txt file, for manual investigation as required.
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "Dip_IBD.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Set file input/output paths ---------------------------------------------

# Directory and rootfilename of SNP data
# i.e. where PLINK files are "xxx.bed", "xxx.bim", "xxx.fam", the rootfilename is "xxx"
SNPdata.filepath <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/Parent Pedigree Analysis/Inputs/Dip"
SNPdata.rootfilename <- "Dip_Data"

# PLINK .genome filename (file should be same directory as SNP data) 
plink.genome.filename <- "plink.genome"

# Output filenames for duo and trio results (.txt files). 
# Can specify which directory these are saved in.
filename.trio.output <- "C:/Users/curly/Desktop/Apple Genotyping/Results/Parent Pedigree Analysis/Diploid Results/trios.txt"
filename.duo.output <- "C:/Users/curly/Desktop/Apple Genotyping/Results/Parent Pedigree Analysis/Diploid Results/duos.txt"

# Directory containing supplementary script "rutilsHMuranty_pedigree.R" 
function.filepath  <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/Parent Pedigree Analysis/Muranty_Scripts"


# Load packages and supplementary script ----------------------------------

# Note - requires snpStats
# if not installed, run:
#    install.packages("BiocManager")
#    BiocManager::install("snpStats")

#Load packages
library(snpStats)

#Load supplementary script
source(paste(function.filepath, "rutilsHMuranty_pedigree.R", sep = "/"))


# Load SNP data and PLINK IBD results -------------------------------------

data <- read.plink(paste(SNPdata.filepath, SNPdata.rootfilename, sep = "/"))

genome.df <- read.table(plink.genome.filename,
                        header = T, stringsAsFactors = FALSE)


# Extract SNP data --------------------------------------------------------

genot.mat <- matrix(as.integer(data$genotypes@.Data),
                    nrow = nrow(data$genotypes@.Data),
                    dimnames = dimnames(data$genotypes@.Data))
genot.mat[genot.mat == 0] <- NA
genot.mat <- genot.mat - 1
genotAB.mat <- gt012toAB(genot.mat)
rm(genot.mat)
nb.markers <- ncol(genotAB.mat)

genome.df <- genome.df[order(genome.df$PI_HAT, decreasing = TRUE),]


# Test all duos present in PLINK IBD results ------------------------------

# As "--min 0.4" was used in PLINK, this will be  duos with PI_HAT larger than 0.4

duos.to.test.df <- genome.df[, c("IID1","IID2")]
duos.to.test.df <- duos.to.test.df[order(duos.to.test.df$IID1),]

duos.to.test.df <- cbind(duos.to.test.df,
                         menderr = rep(NA, nrow(duos.to.test.df)))

all.genot1 <- unique(duos.to.test.df$IID1)
nb.genot1 <- length(all.genot1)

for(ind1 in seq(nb.genot1)) {
  name.ind1 <- all.genot1[ind1]
  AA.ind1 <- genotAB.mat[name.ind1, ] == "AA"
  BB.ind1 <- genotAB.mat[name.ind1, ] == "BB"
  all.ind2 <- duos.to.test.df$IID2[duos.to.test.df$IID1 == name.ind1]
  if(length(all.ind2) > 1) {
    duos.to.test.df$menderr[duos.to.test.df$IID1 == name.ind1] <-
      apply(genotAB.mat[all.ind2, AA.ind1] == "BB", 1, sum) +
      apply(genotAB.mat[all.ind2, BB.ind1] == "AA", 1, sum)
  } else {
    duos.to.test.df$menderr[duos.to.test.df$IID1 == name.ind1] <-
      sum(genotAB.mat[all.ind2, AA.ind1] == "BB") +
      sum(genotAB.mat[all.ind2, BB.ind1] == "AA")
  }
}
rm(ind1, name.ind1, AA.ind1, BB.ind1, all.ind2)

duos.to.test.df <- duos.to.test.df[order(duos.to.test.df$menderr),]

# Plot ranked numbers of Mendelian to decide the threshold for declaring parent-offspring duos
plot(duos.to.test.df$menderr, ylab = "nb mendelian errors")

zoom.region <- 25000:35000 # Adjust the region of the plot to zoom in to better decide the threshold
plot(duos.to.test.df$menderr[zoom.region], ylab = "nb mendelian errors")

# Define threshold for parent-offspring duos and apply it
threshold.PO <- nb.markers * 0.03 # adjust according to data
duos.to.test.df <- cbind(duos.to.test.df,
                         CP.infer = rep(0, nrow(duos.to.test.df)))
duos.to.test.df$CP.infer[duos.to.test.df$menderr < threshold.PO] <- 1


# Build a list of parent-child duos ---------------------------------------

nb.C.P.by.indiv <-
  table(c(duos.to.test.df$IID1[duos.to.test.df$CP.infer == 1],
          duos.to.test.df$IID2[duos.to.test.df$CP.infer == 1]))

CP.list <- vector("list", length = length(nb.C.P.by.indiv))
names(CP.list) <- names(nb.C.P.by.indiv)

for(indiv in names(nb.C.P.by.indiv))
  CP.list[[indiv]] <-
  sort(c(duos.to.test.df[duos.to.test.df$CP.infer == 1 &
                           duos.to.test.df$IID1 == indiv, "IID2"],
         duos.to.test.df[duos.to.test.df$CP.infer == 1 &
                           duos.to.test.df$IID2 == indiv, "IID1"]))


# Test potential trios for Mendelian errors -------------------------------

trios.test.df <- data.frame(
  array(NA, dim = c(0, 4),
        dimnames = list(NULL, c("indiv", "Parent1", "Parent2", "nb.inc.all.mk"))),
  stringsAsFactors = FALSE)

for(indiv in #1:10) {
    seq_along(CP.list)) {
  nb.C.P.this.ind <- length(CP.list[[indiv]])
  if(nb.C.P.this.ind > 1) {
    indiv.name <- names(CP.list)[indiv]
    indiv.gt <- genotAB.mat[indiv.name,]
    for(par1 in seq(nb.C.P.this.ind - 1)) {
      par1.name <- CP.list[[indiv]][par1]
      par1.gt <- genotAB.mat[par1.name,]
      for(par2 in seq(from = par1 + 1, to = nb.C.P.this.ind)) {
        par2.name <- CP.list[[indiv]][par2]
        par2.gt <- genotAB.mat[par2.name,]
        mend.inc <- mendelian.checking.by.ind(
          indiv.gt, par1.gt, par2.gt)
        trio.test.here <- data.frame(
          indiv = indiv.name,
          Parent1 = par1.name,
          Parent2 = par2.name,
          nb.inc.all.mk = sum(mend.inc),
          stringsAsFactors = FALSE)
        trios.test.df <- rbind(trios.test.df, trio.test.here)
      }
    }
  }
}
rm(indiv, nb.C.P.this.ind, indiv.name, indiv.gt, par1, par1.name, par1.gt,
   par2, par2.name, par2.gt, mend.inc, trio.test.here)


# Export duo and trio results ---------------------------------------------

write.table(trios.test.df,
            file = filename.trio.output,
            sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE )

write.table(duos.to.test.df,
            file = filename.duo.output,
            sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE )
