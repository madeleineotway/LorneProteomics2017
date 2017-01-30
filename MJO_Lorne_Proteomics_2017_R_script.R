###############################
# MJO
# Lorne Proteomics 2017
# motway@cmri.org.au
# R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
# Bioconductor version 3.4 (BiocInstaller 1.24.0)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
###############################

#install.packages("VennDiagram")
library(VennDiagram)

#install.packages("Peptides")
library(Peptides)

#install.packages("fBasics")
library(fBasics)

# source("http://bioconductor.org/biocLite.R")
# biocLite("plyr")
library(plyr)

###########################

# Import data frames for peptides and FDR
pep.all <-
  read.table(
    "/Users/madeleineotway/Desktop/MJO_Lorne_Proteomics_2017_peptide_areas.txt",
    sep = "\t", header = T)
fdr.all <-
  read.table(
    "/Users/madeleineotway/Desktop/MJO_Lorne_Proteomics_2017_FDR.txt",
    sep = "\t", header = T)

# Create and rename columns in new data frame
pep.lmh <- data.frame("protein" = pep.all[,1], 
                      "peptide" = pep.all[,2],
                      "precursor" = pep.all[,3],
                      "low.1" = pep.all[,6],
                      "low.2" = pep.all[,7],
                      "low.3" = pep.all[,8],
                      "low.4" = pep.all[,9],
                      "low.5" = pep.all[,10],
                      "low.6" = pep.all[,11],
                      "med.1" = pep.all[,12],
                      "med.2" = pep.all[,13],
                      "med.3" = pep.all[,14],
                      "med.4" = pep.all[,15],
                      "med.5" = pep.all[,16],
                      "med.6" = pep.all[,17],
                      "high.1" = pep.all[,18],
                      "high.2" = pep.all[,19],
                      "high.3" = pep.all[,20],
                      "high.4" = pep.all[,21],
                      "high.5" = pep.all[,22],
                      "high.6" = pep.all[,23])

# Create and rename columns in new data frame
fdr.lmh <- data.frame("protein" = fdr.all[,1], 
                      "peptide" = fdr.all[,2],
                      "precursor" = fdr.all[,4],
                      "low.1" = fdr.all[,8],
                      "low.2" = fdr.all[,9],
                      "low.3" = fdr.all[,10],
                      "low.4" = fdr.all[,11],
                      "low.5" = fdr.all[,12],
                      "low.6" = fdr.all[,13],
                      "med.1" = fdr.all[,14],
                      "med.2" = fdr.all[,15],
                      "med.3" = fdr.all[,16],
                      "med.4" = fdr.all[,17],
                      "med.5" = fdr.all[,18],
                      "med.6" = fdr.all[,19],
                      "high.1" = fdr.all[,20],
                      "high.2" = fdr.all[,21],
                      "high.3" = fdr.all[,22],
                      "high.4" = fdr.all[,23],
                      "high.5" = fdr.all[,24],
                      "high.6" = fdr.all[,25],
                      "Decoy" = fdr.all[,7])
                        

# Create new data frame with all non decoy proteins
fdr.decoy <- fdr.lmh[!(fdr.lmh$Decoy=="True"),]

# Match proteins, peps and precursors between FDR and pep lists
fdr.decoy$fdr.pep <-
  ifelse(
    fdr.decoy$protein %in% pep.lmh$protein &
      fdr.decoy$peptide %in% pep.lmh$peptide &
      fdr.decoy$precursor %in% pep.lmh$precursor, 1, 0)

# Create new data frame with only the rows that match between FDR and pep lists
fdr.match <- fdr.decoy[!(fdr.decoy$fdr.pep == "0"),]

# Sort data frames by protein, then pep, then precursor
sort.fdr <- arrange(fdr.match, protein, peptide, precursor)
sort.pep <- arrange(pep.lmh, protein, peptide, precursor)

# Create new data frame for peptide FDR
pep.fdr.01 <- sort.pep

# Function for FDR calculation at 1%
fdr.01 <- function(x, y, table.in) {
  # create pep.fdr.01 from sort.pep
  # table.in = pep.fdr.01
  # x= first sample vector
  # y= last sample vector
  
  table.in <- pep.fdr.01
  table.in [,x] <- ifelse(sort.fdr[,x] <= 0.01, sort.pep[,x], NA)
  pep.fdr.01 <- table.in
  
  for (i in x:y) {
    table.in <- pep.fdr.01
    table.in [,i] <- ifelse(sort.fdr[,i] <= 0.01, sort.pep[,i], NA)
    pep.fdr.01 <- table.in 
  }
  return(pep.fdr.01)
}

pep.fdr.01 <- fdr.01(4, 21, pep.fdr.01)


# New data frame per group, with group specific missing values removed
low.mv <-
  pep.fdr.01[(rowSums(is.na(pep.fdr.01[4:9])) <= 0) ,]
low.mv[10:21] <- NULL
med.mv <-
  pep.fdr.01[(rowSums(is.na(pep.fdr.01[10:15])) <= 0) ,]
med.mv[16:21] <- NULL
med.mv[4:9] <- NULL
high.mv <-
  pep.fdr.01[(rowSums(is.na(pep.fdr.01[16:21])) <= 0) ,]
high.mv[4:15] <- NULL

# Calculate and remove all peptides with a CV >20%
low.mv$mean <- rowMeans(low.mv[4:9])
low.mv$sd <- rowStdevs(low.mv[4:9])
low.mv$cv <- (low.mv$sd / low.mv$mean * 100)
low.mv$cv.20 <- ifelse ((low.mv$cv < 20), 0, 1)

med.mv$mean <- rowMeans(med.mv[4:9])
med.mv$sd <- rowStdevs(med.mv[4:9])
med.mv$cv <- (med.mv$sd / med.mv$mean * 100)
med.mv$cv.20 <- ifelse ((med.mv$cv < 20), 0, 1)

high.mv$mean <- rowMeans(high.mv[4:9])
high.mv$sd <- rowStdevs(high.mv[4:9])
high.mv$cv <- (high.mv$sd / high.mv$mean * 100)
high.mv$cv.20 <- ifelse ((high.mv$cv < 20), 0, 1)

low.cv.20 <- data.frame(low.mv[low.mv$cv.20 != 1,])
med.cv.20 <- data.frame(med.mv[med.mv$cv.20 != 1,])
high.cv.20 <- data.frame(high.mv[high.mv$cv.20 != 1,])

# Venn diagram of proteins
pro.venn.list <-
  list("Low (2 μg)" = low.cv.20$protein,
       "Medium (4 μg)" = med.cv.20$protein,
       "High (6 μg)" = high.cv.20$protein)

venn.diagram(
  pro.venn.list,
  filename = "/Users/madeleineotway/Desktop/Pro_Venn_cv_pc.png",
  imagetype = "png",
  print.mode = "percent",
  sigdigs = 1,
  fill = c("darkseagreen3", "lightgoldenrod2", "pink2"),
  alpha = 0.5,
  cex = 1.5,
  cat.col = "black",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.dist = c(0.04, 0.04, 0.03),
  cat.pos = 0 )

# Venn diagram of peptides
pep.venn.list <-
  list("Low (2 μg)" = low.cv.20$peptide,
       "Medium (4 μg)" = med.cv.20$peptide,
       "High (6 μg)" = high.cv.20$peptide)

venn.diagram(
  pep.venn.list,
  filename = "/Users/madeleineotway/Desktop/Pep_Venn_cv_pc.png",
  imagetype = "png",
  print.mode = "percent",
  sigdigs = 1,
  fill = c("darkseagreen3", "lightgoldenrod2", "pink2"),
  alpha = 0.5,
  cex = 1.5,
  cat.col = "black",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.dist = c(0.04, 0.04, 0.03),
  cat.pos = 0 )

# New data frame with the frequency of peptides per protein
low.pro.count <- count(low.cv.20$protein)
low.pro.freq <- data.frame(pep.per.pro = c(1, 2, 3, 4, 5, 6), 
                           freq = c((length(which(low.pro.count$freq == 1))),
                                    (length(which(low.pro.count$freq == 2))),
                                    (length(which(low.pro.count$freq == 3))),
                                    (length(which(low.pro.count$freq == 4))),
                                    (length(which(low.pro.count$freq == 5))),
                                    (length(which(low.pro.count$freq == 6))) ))

med.pro.count <- count(med.cv.20$protein)
med.pro.freq <- data.frame(pep.per.pro = c(1, 2, 3, 4, 5, 6), 
                           freq = c((length(which(med.pro.count$freq == 1))),
                                    (length(which(med.pro.count$freq == 2))),
                                    (length(which(med.pro.count$freq == 3))),
                                    (length(which(med.pro.count$freq == 4))),
                                    (length(which(med.pro.count$freq == 5))),
                                    (length(which(med.pro.count$freq == 6))) ))

high.pro.count <- count(high.cv.20$protein)
high.pro.freq <- data.frame(pep.per.pro = c(1, 2, 3, 4, 5, 6), 
                           freq = c((length(which(high.pro.count$freq == 1))),
                                    (length(which(high.pro.count$freq == 2))),
                                    (length(which(high.pro.count$freq == 3))),
                                    (length(which(high.pro.count$freq == 4))),
                                    (length(which(high.pro.count$freq == 5))),
                                    (length(which(high.pro.count$freq == 6))) ))

# Plot the disriputiopn of the proteins
plot(x = c(low.pro.freq$pep.per.pro,
           med.pro.freq$pep.per.pro,
           high.pro.freq$pep.per.pro),
     y = c(low.pro.freq$freq,
           med.pro.freq$freq,
           high.pro.freq$freq),
     col = c(rep("seagreen", 6), rep("goldenrod", 6), rep("hotpink", 6)),
     pch = 19, cex = 1.5,
     xlab = "# Peptides/Protein",
     ylab = "Frequency")
par(bty = "l")
