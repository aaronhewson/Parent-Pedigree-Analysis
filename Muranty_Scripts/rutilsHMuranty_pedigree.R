# Hélène Muranty, IRHS, 2019-09-12
# auxiliary functions for pedigree analysis

gt012toAB <- function(gt) {
# objective: reformat genotypic data from 0,1,2,NA format to AA, AB, BB, -- format 
# testing 
# gt012toAB(cbind(c(0,0,0,1,1,NA), c(1,2,1,2,NA,2), c(2,1,0,0,1,2)))
# should give
# cbind(c("AA","AA","AA","AB","AB","--"), c("AB","BB","AB","BB","--","BB"), c("BB","AB","AA","AA","AB","BB"))

  if(!(is.matrix(gt))) {
    stop("gt012toAB function intented for matrix\n")
  }
  gtAB <- array("--", dim = dim(gt), dimnames = dimnames(gt))
  gtAB[!is.na(gt) & gt == 0] <- "AA"
  gtAB[!is.na(gt) & gt == 1] <- "AB"
  gtAB[!is.na(gt) & gt == 2] <- "BB"
  return(gtAB)
}

mendelian.checking.by.ind <- function(indiv.data, parent1.data, parent2.data) {
# objective: identify Mendelian inconsistent errors in parent-offspring trios
 if((length(indiv.data) != length(parent1.data)) | 
    (length(indiv.data) != length(parent2.data)) |
    (length(parent1.data) != length(parent2.data)) ) {
    stop("data for candidate offspring, parent1 and parent2 should have the same length")
  }
  if(any(is.na(match(c(indiv.data, parent1.data, parent2.data), c("AA", "AB", "BB", "--"))))) {
    stop("this function is intented for data in AA, AB, BB, -- format")
  }
  indiv.mend.error <- rep(NA, length(indiv.data))
  names(indiv.mend.error) <- names(indiv.data)
  
  if(any(parent1.data == "AA" & parent2.data == "AA")) {
    indiv.mend.error[parent1.data == "AA" & parent2.data == "AA"] <-
      (indiv.data[parent1.data == "AA" & parent2.data == "AA"] == "AB") |
      (indiv.data[parent1.data == "AA" & parent2.data == "AA"] == "BB")
  }
  if(any(parent1.data == "AA" & parent2.data == "AB")) {
    indiv.mend.error[parent1.data == "AA" & parent2.data == "AB"] <-
      indiv.data[parent1.data == "AA" & parent2.data == "AB"] == "BB"
  }
  if(any(parent1.data == "AA" & parent2.data == "BB")) {
    indiv.mend.error[parent1.data == "AA" & parent2.data == "BB"] <-
      (indiv.data[parent1.data == "AA" & parent2.data == "BB"] == "AA") |
      (indiv.data[parent1.data == "AA" & parent2.data == "BB"] == "BB")
  }
  if(any(parent1.data == "AA" & parent2.data == "--")) {
    indiv.mend.error[parent1.data == "AA" & parent2.data == "--"] <-
      indiv.data[parent1.data == "AA" & parent2.data == "--"] == "BB"
  }
  if(any(parent1.data == "AB" & parent2.data == "AA")) {
    indiv.mend.error[parent1.data == "AB" & parent2.data == "AA"] <-
      indiv.data[parent1.data == "AB" & parent2.data == "AA"] == "BB"
  }
  if(any(parent1.data == "AB" & 
           (parent2.data == "AB" | parent2.data == "--"))) {
    indiv.mend.error[parent1.data == "AB" & 
                    (parent2.data == "AB" | parent2.data == "--")] <- FALSE
  }
  if(any(parent1.data == "AB" & parent2.data == "BB")) {
    indiv.mend.error[parent1.data == "AB" & parent2.data == "BB"] <-
      indiv.data[parent1.data == "AB" & parent2.data == "BB"] == "AA"
  }
  if(any(parent1.data == "BB" & parent2.data == "AA")) {
    indiv.mend.error[parent1.data == "BB" & parent2.data == "AA"] <-
      (indiv.data[parent1.data == "BB" & parent2.data == "AA"] == "AA") |
      (indiv.data[parent1.data == "BB" & parent2.data == "AA"] == "BB")
  }
  if(any(parent1.data == "BB" & parent2.data == "AB")) {
    indiv.mend.error[parent1.data == "BB" & parent2.data == "AB"] <-
      indiv.data[parent1.data == "BB" & parent2.data == "AB"] == "AA"
  }
  if(any(parent1.data == "BB" & parent2.data == "BB")){
    indiv.mend.error[parent1.data == "BB" & parent2.data == "BB"] <- 
      (indiv.data[parent1.data == "BB" & parent2.data == "BB"] == "AA") |
      (indiv.data[parent1.data == "BB" & parent2.data == "BB"] == "AB")
  }
  if(any(parent1.data == "BB" & parent2.data == "--")){
    indiv.mend.error[parent1.data == "BB" & parent2.data == "--"] <- 
      (indiv.data[parent1.data == "BB" & parent2.data == "--"] == "AA") 
  }
  if(any(parent1.data == "--" & parent2.data == "AA")) {
    indiv.mend.error[parent1.data == "--" & parent2.data == "AA"] <-
      (indiv.data[parent1.data == "--" & parent2.data == "AA"] == "BB")
  }
  if(any(parent1.data == "--" & 
           (parent2.data == "AB" | parent2.data == "--"))) {
    indiv.mend.error[parent1.data == "--" & 
                    (parent2.data == "AB" | parent2.data == "--")] <- FALSE
  }
  if(any(parent1.data == "--" & parent2.data == "BB")) {
    indiv.mend.error[parent1.data == "--" & parent2.data == "BB"] <-
      (indiv.data[parent1.data == "--" & parent2.data == "BB"] == "AA")
  }
  
  return(indiv.mend.error)
}

