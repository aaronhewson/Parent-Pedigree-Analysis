# Hélène Muranty, IRHS, 2017-10-27
# identification of grand-parent pairs for identified parent-offspring duos
# not included in identified trios, based on Mendelian inconsistent error counts
# initial step: create objects indicating status of identified duos
# this is OK with 25K SNPs and ~400 duos

## note: require snpStats
# to install it
#    install.packages("BiocManager")
#    BiocManager::install("snpStats")

# REQUIRED DATA & FORMAT
# the script is built for SNP marker data in plink binary format
# (see http://zzz.bwh.harvard.edu/plink/data.shtml#bed)

# settings
# put here directory and rootfilename of your SNP data
# i.e. if your three files are "xxx.bed", "xxx.bim", "xxx.fam", the rootfilename is "xxx"
SNPdata.filepath <- "..."
SNPdata.rootfilename <- "..."
# put here the name of the file containing the list of duos to test
# i.e. those not included in trios
# this file should be in the same directory as the SNP data
duos.to.test.filename <- "..."
# put here the name of the file with the status objects, should be of ".RData" type
# this file will be exported in the directory where the SNP data are
duos.status.filename <- "..."
# put here the directory where you saved the file rutilsHMuranty_pedigree.R
function.file.path  <- "..."

# load libraries and functions
library(snpStats) # read SNP data in plink binary format
# source the auxiliary functions
source(paste(function.file.path, "rutilsHMuranty_pedigree.R", sep = "/"))

# 1. load the snp data
# snp data
data <- read.plink(paste(SNPdata.filepath, SNPdata.rootfilename, sep = "/"))

genot.mat <- matrix(as.integer(data$genotypes@.Data),
                    nrow = nrow(data$genotypes@.Data),
                    dimnames = dimnames(data$genotypes@.Data))
genot.mat[genot.mat == 0] <- NA
genot.mat <- genot.mat - 1
genotAB.mat <- gt012toAB(genot.mat)
rm(genot.mat)

duos.to.test.df <- read.table(
    file = paste(SNPdata.filepath,duos.to.test.filename, sep = "/"),
    header = TRUE,
    stringsAsFactors = FALSE)
nb.duos.to.test <- nrow(duos.to.test.df)

#build empty objects
duos.status.GP.notBB <- array(NA, dim = c(2 * nb.duos.to.test, ncol(genotAB.mat)),
    dimnames = list(c(paste(duos.to.test.df[, "IID1"],
                            duos.to.test.df[, "IID2"], sep = "-"),
                      paste(duos.to.test.df[, "IID2"],
                            duos.to.test.df[, "IID1"], sep = "-")),
                    colnames(genotAB.mat)))
duos.status.GP.notAA <- array(NA, dim = c(2 * nb.duos.to.test, ncol(genotAB.mat)),
      dimnames = list(c(paste(duos.to.test.df[, "IID1"],
                              duos.to.test.df[, "IID2"], sep = "-"),
                        paste(duos.to.test.df[, "IID2"],
                              duos.to.test.df[, "IID1"], sep = "-")),
                      colnames(genotAB.mat)))
nb.inform.mk <- rep(NA, 2 * nb.duos.to.test)
names(nb.inform.mk) <- c(paste(duos.to.test.df[, "IID1"],
                               duos.to.test.df[, "IID2"], sep = "-"),
                         paste(duos.to.test.df[, "IID2"],
                               duos.to.test.df[, "IID1"], sep = "-"))
start.time.out <- Sys.time()
for(duo.index in #1:100) {
                 seq(nb.duos.to.test)) {
#  start.time <- Sys.time()
  iid1.name <- duos.to.test.df[duo.index, "IID1"]
  iid2.name <- duos.to.test.df[duo.index, "IID2"]
  iid1.data <- genotAB.mat[iid1.name,]
  iid2.data <- genotAB.mat[iid2.name,]
# individual named in column IID1 is the parent
# individual named in column IID2 is the offspring,
  duos.status.GP.notBB[duo.index,] <-
    (iid2.data == "AA" & iid1.data != "BB") |
    (iid2.data == "AB" & iid1.data == "BB")
  duos.status.GP.notAA[duo.index,] <-
    (iid2.data == "BB" & iid1.data != "AA") |
    (iid2.data == "AB" & iid1.data == "AA")
  nb.inform.mk[duo.index] <- sum(duos.status.GP.notBB[duo.index,]) +
                               sum(duos.status.GP.notAA[duo.index,])
# individual named in column IID1 is the offspring,
# individual named in column IID2 is the parent
  duos.status.GP.notBB[duo.index + nb.duos.to.test,] <-
    (iid1.data == "AA" & iid2.data != "BB") |
    (iid1.data == "AB" & iid2.data == "BB")
  duos.status.GP.notAA[duo.index + nb.duos.to.test,] <-
    (iid1.data == "BB" & iid2.data != "AA") |
    (iid1.data == "AB" & iid2.data == "AA")
  nb.inform.mk[duo.index + nb.duos.to.test] <-
    sum(duos.status.GP.notBB[duo.index + nb.duos.to.test,]) +
      sum(duos.status.GP.notAA[duo.index + nb.duos.to.test,])
#  finish.time <- Sys.time()
#  cat("duo", iid1.name, "-", iid2.name,
#      "examined in", finish.time-start.time, fill = T)
}
finish.time.out <- Sys.time()
cat("all duos examined in", finish.time.out - start.time.out, fill = T)

object.size(duos.status.GP.notBB)
object.size(duos.status.GP.notAA)
object.size(nb.inform.mk)

save(duos.status.GP.notBB, duos.status.GP.notAA, nb.inform.mk,
     file = paste(SNPdata.filepath, duos.status.filename, sep = "/"))
