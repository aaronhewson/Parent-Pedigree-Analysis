# Hélène Muranty, IRHS, 2019-09-19
# identification of grand-parent pairs for identified parent-offspring duos
# not included in identified trios, based on Mendelian inconsistent error counts
# this script enables to divide the work in several parts
# by working on groups of individuals as candidate GP1
# a kind of parallelization for dummies
# some messages indicating the progression are put on the user interface
# you could save these to a file when running in non-interactive interface

## note: require snpStats
# to install it
#    install.packages("BiocManager")
#    BiocManager::install("snpStats")

# REQUIRED DATA & FORMAT
# the script is built for SNP marker data in plink binary format
# (see http://zzz.bwh.harvard.edu/plink/data.shtml#bed)
# first run create_status_duos.R to obtain one of the input file

# settings
# put here directory and rootfilename of your SNP data
# i.e. if your three files are "xxx.bed", "xxx.bim", "xxx.fam", the rootfilename is "xxx"
SNPdata.filepath <- "..."
SNPdata.rootfilename <- "..."
# put here the name of the file containing the status results for the duos
# you want to test, obtained with create_status_duos.R
# this file should be in the same directory as the SNP data
duos.status.filename <- "..."
# put here the directory where you want to export the results
results.file.path <- "..."
# put here two parts for the result file
result.filename.part1 <- "..."
result.filename.part2 <- "..."
# put here the directory where you saved the file rutilsHMuranty_pedigree.R
function.file.path  <- "..."

# define an acceptable error threshold for exporting a result
error.threshold <- 0.025

# define a part of the individuals to test as GP1
# for testing a total of ~1400 individuals, I used values
# start.i  1 21 41  71 101 151 251 401 601  901 1201
# start.j 20 40 70 100 150 250 400 600 900 1200 nb.genot -1
# nbgenot defined below
start.i <- 1
end.i <- 20

# load libraries and functions
library(snpStats) # read SNP data in plink binary format
# source the auxiliary functions
source(paste(function.file.path, "rutilsHMuranty_pedigree.R", sep = "/"))

# 1. load the snp data
# snp data
data <- read.plink(paste(SNPdata.filepath, SNPdata.rootfilename, sep = "/"))


duos.to.test.df <- read.table(
    file = paste(SNPdata.filepath, "duos_not_in_trios_SNP_FBo_GF_plus_253K_FBKey_red.txt", sep = "/"),
    header = TRUE,
    stringsAsFactors = FALSE)
load(file = paste(datafile.path, duos.status.filename, sep = "/"))

# 2. extract genotypic data
genot.mat <- matrix(as.integer(data$genotypes@.Data),
                    nrow = nrow(data$genotypes@.Data),
                    dimnames = dimnames(data$genotypes@.Data))
genot.mat[genot.mat == 0] <- NA
genot.mat <- genot.mat - 1
genotAB.mat <- gt012toAB(genot.mat)
rm(genot.mat)

all.indiv <- rownames(genotAB.mat)
nb.genot <- length(all.indiv)

nb.duos.to.test <- nrow(duos.to.test.df)


# define the result file name
result.filename <- paste(result.filename.part1, start.i, "_", end.i,
                          result.filename.part2, ".txt", sep = "")
result.file <- paste(results.file.path, results.filename, sep = "/")
# create an empty result file
GP.duo.test.df <- data.frame(
  array(NA, dim = c(0, 6),
        dimnames = list(NULL, c("offspring", "parent", "GP1", "GP2", "nb.inform.mk", "nb.inc.all.mk"))),
  stringsAsFactors = FALSE)
write.table(GP.duo.test.df, file = result.file, quote = FALSE)

# start loop on GP1
start.time.out <- Sys.time()
for(ind.i in seq(from = start.i, to = end.i)) {
  start.time.in <- Sys.time()
  ind.i.id <- all.indiv[ind.i]
  ind.i.genot.mat <- genotAB.mat[ind.i.id, ]
# all other individuals are potential GP2
  all.ind.j <- seq(from = ind.i + 1, to = #end.i + 1)
                                          nb.genot)
  all.ind.j.id <- all.indiv[all.ind.j]
  all.ind.j.genot.mat <- genotAB.mat[all.ind.j.id,]
  if(length(all.ind.j) > 1) {
# syntax when several GP2 are to be tested
    pairs.status.BB.here <- apply(all.ind.j.genot.mat, 1, function(x, y) {
      x == "BB" & y == "BB" }, y = ind.i.genot.mat)
    pairs.status.AA.here <- apply(all.ind.j.genot.mat, 1, function(x, y) {
      x == "AA" & y == "AA" }, y = ind.i.genot.mat)
    pairs.status.inform.here <- apply(all.ind.j.genot.mat, 1, function(x, y) {
      x != "--" & y != "--" }, y = ind.i.genot.mat)
    colnames(pairs.status.BB.here) <- paste(ind.i.id, all.ind.j.id, sep = "-X-")
    colnames(pairs.status.AA.here) <- paste(ind.i.id, all.ind.j.id, sep = "-X-")
    colnames(pairs.status.inform.here) <- paste(ind.i.id, all.ind.j.id, sep = "-X-")
    for(duo.index in #1:100) {
                     seq(nb.duos.to.test)) {
# individual named in column IID2 is the offspring,
# individual named in column IID1 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID2"]
# continue only if the focal GP1 is not the focal offspring
      if(ind.i.id != offspring.ID) {
# remove the focal offspring from the list of potential GP2
        offspring.in.all.ind.j <- which(all.ind.j.id == offspring.ID)
        if(length(offspring.in.all.ind.j) > 0) {
          if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) > 1) {
# syntax when the focal offspring has to be removed from the list of GP2 to test
# but there are still several GP2 to test
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                               - offspring.in.all.ind.j], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                   - offspring.in.all.ind.j], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID1"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[-offspring.in.all.ind.j][where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) == 1) {
# syntax when only one GP2 is tested after removing the focal offspring
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID1"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id[- offspring.in.all.ind.j],
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        } else { # the focal offspring is not in the list of GP2 to test
          if((length(all.ind.j.id)) > 1) { # there are several GP2 to test
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                  ], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                  ], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                               ], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                   ], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID1"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if(length(all.ind.j.id) == 1) { # there is only one GP2 to test
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                ], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                ], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                                ], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                ], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID1"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id,
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        }
      }

# individual named in column IID1 is the offspring,
# individual named in column IID2 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID1"]
# continue only if the focal GP1 is not the focal offspring
      if(ind.i.id != offspring.ID) {
# remove the focal offspring from the list of potential GP2
        offspring.in.all.ind.j <- which(all.ind.j.id == offspring.ID)
        if(length(offspring.in.all.ind.j) > 0) {
# syntax when the focal offspring has to be removed from the list of GP2 to test
# but there are still several GP2 to test
          if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) > 1) {
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                               - offspring.in.all.ind.j], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                   - offspring.in.all.ind.j], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID2"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[-offspring.in.all.ind.j][where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) == 1) {
# syntax when only one GP2 is tested after removing the focal offspring
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID2"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id[- offspring.in.all.ind.j],
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        } else { # the focal offspring is not in the list of GP2 to test
          if((length(all.ind.j.id)) > 1) { # there are several GP2 to test
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                  ], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                  ], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                               ], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                   ], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID2"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if(length(all.ind.j.id) == 1) { # there is only one GP2 to test
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                ], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                ], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                ], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                ], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID2"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id,
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        }
      }
    }
  } else { # there is only one GP2 to test
    pairs.status.BB.here <- all.ind.j.genot.mat == "BB" & ind.i.genot.mat == "BB"
    pairs.status.AA.here <- all.ind.j.genot.mat == "AA" & ind.i.genot.mat == "AA"
    pairs.status.inform.here <- all.ind.j.genot.mat != "--" & ind.i.genot.mat != "--"
    for(duo.index in #1:100) {
                     seq(nb.duos.to.test)) {
# individual named in column IID2 is the offspring,
# individual named in column IID1 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID2"]
# continue only if neither the focal GP1 not the focal GP2 are the focal offspring
      if((ind.i.id != offspring.ID) & (all.ind.j.id != offspring.ID)) {
        nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index,]], na.rm = T)
        nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index,]], na.rm = T)
        nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index,]], na.rm = T) +
                       sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index,]], na.rm = T)
        error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
        if(error.rate < error.threshold) {
          GP.duo.test.here <- data.frame(
            offspring = offspring.ID,
            parent = duos.to.test.df[duo.index, "IID1"],
            GP1 = ind.i.id,
            GP2 = all.ind.j.id,
            nb.inform.mk = nb.inform,
            nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
            stringsAsFactors = FALSE)
          write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                      row.names = FALSE, quote = FALSE, append = T)
        }
      }
# individual named in column IID1 is the offspring,
# individual named in column IID2 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID1"]
# continue only if neither the focal GP1 not the focal GP2 are the focal offspring
      if((ind.i.id != offspring.ID) & (all.ind.j.id != offspring.ID)) {
        nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[
                                duo.index + nb.duos.to.test,]], na.rm = T)
        nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[
                                duo.index + nb.duos.to.test,]], na.rm = T)
        nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[
                                      duo.index + nb.duos.to.test,]], na.rm = T) +
                       sum(pairs.status.inform.here[duos.status.GP.notAA[
                                      duo.index + nb.duos.to.test,]], na.rm = T)
        error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
        if(error.rate < error.threshold) {
          GP.duo.test.here <- data.frame(
            offspring = offspring.ID,
            parent = duos.to.test.df[duo.index, "IID2"],
            GP1 = ind.i.id,
            GP2 = all.ind.j.id,
            nb.inform.mk = nb.inform,
            nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
            stringsAsFactors = FALSE)
          write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                      row.names = FALSE, quote = FALSE, append = T)
        }
      }
    }
  }
  finish.time.in <- Sys.time()
  cat("all pairs involving", ind.i.id,
      "examined in", finish.time.in - start.time.in, fill = T)
}
finish.time.out <- Sys.time()

cat("all GP pairs examined in", finish.time.out - start.time.out, fill = T)
