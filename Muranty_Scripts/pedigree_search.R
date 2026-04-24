# Hélène Muranty, IRHS, 2019-09-12
# pedigree identification based on Mendelian inconsistent error counts

## note: require snpStats
# to install it
#    install.packages("BiocManager")
#    BiocManager::install("snpStats")

# REQUIRED DATA & FORMAT
# the script is built for SNP marker data in plink binary format
# (see http://zzz.bwh.harvard.edu/plink/data.shtml#bed)
# first run plink on these data to estimate pairwise IBD
# with the "--minX" option to reduce output file size and number of duos to test (e.g. with X = 0.4)
# (see http://zzz.bwh.harvard.edu/plink/ibdibs.shtml#genome)

# settings
# put here directory and rootfilename of your SNP data
# i.e. if your three files are "xxx.bed", "xxx.bim", "xxx.fam", the rootfilename is "xxx"
SNPdata.filepath <- "..."
SNPdata.rootfilename <- "..."
# put here the name of the output file obtained by running pairwise IBD estimation with plink
plink.genome.filename <- "..."
# this file is supposed to be in the same directory as the SNP data
# put here the filenames for trio and duo results (text files)
# the result files will be in the working directory
filename.trio.output <- "..."
filename.duo.output <- "..."
# put here the directory where you saved the file rutilsHMuranty_pedigree.R
function.filepath  <- "..."

# load libraries and functions
library(snpStats) # read SNP data in plink binary format
# source the auxiliary functions
source(paste(function.filepath, "rutilsHMuranty_pedigree.R", sep = "/"))

# 1. load the snp data and the results of pairwise IBD estimations
data <- read.plink(paste(SNPdata.filepath, SNPdata.rootfilename, sep = "/"))

genome.df <- read.table(plink.genome.filename,
                        header = T, stringsAsFactors = FALSE)

# 2. extract genotypic data
genot.mat <- matrix(as.integer(data$genotypes@.Data),
                    nrow = nrow(data$genotypes@.Data),
                    dimnames = dimnames(data$genotypes@.Data))
genot.mat[genot.mat == 0] <- NA
genot.mat <- genot.mat - 1
genotAB.mat <- gt012toAB(genot.mat)
rm(genot.mat)
nb.markers <- ncol(genotAB.mat)

genome.df <- genome.df[order(genome.df$PI_HAT, decreasing = TRUE),]

# 3. test all duos present in the results of plink pairwise IBD estimation
# thus with a PI_HAT larger than X if the --minX option was used
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

# plot ranked numbers of Mendelian to decide the threshold for declaring parent-offspring duos
plot(duos.to.test.df$menderr, ylab = "nb mendelian errors")

zoom.region <- 1000:1500 # adjust the region of the plot to zoom in to better decide the threshold
plot(duos.to.test.df$menderr[zoom.region], ylab = "nb mendelian errors")

# define threshold for parent-offspring duos and apply it to your data
threshold.PO <- nb.markers * 0.005 # adjust according to your data
duos.to.test.df <- cbind(duos.to.test.df,
                              CP.infer = rep(0, nrow(duos.to.test.df)))
duos.to.test.df$CP.infer[duos.to.test.df$menderr < threshold.PO] <- 1

# 4. build a list of parent-child infered relationship
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

# 5. test trios for inconsistencies : a child and its two parents
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

# 6. export the results
write.table(trios.test.df,
            file = filename.trio.output,
            sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE )

write.table(duos.to.test.df,
            file = filename.duo.output,
            sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE )
