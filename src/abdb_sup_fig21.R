#install.packages("ggplot2")

library("ggplot2")
library("tibble")
library("viridis")
AA = "LEU"
chain = "H"
CDR = "CDR-H2"

## some pretty colors
library(MASS)  # in case it is not already loaded 
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

setwd("C:/Users/pprobert/Desktop/TempWorld/");

scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedMotifsPointsOnlyForR.txt", header=TRUE, colClasses=c("character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain", "AA", "Position", "phi","psi"));

listAAs = levels(factor(scatter.data$"AA"))	#"ALA" "ARG" "ASN" "ASP" "CYS" "GLN" "GLU" "GLY" "HIS" "ILE" "LEU"  "LYS" "MET" "PHE" "PRO" "SER" "THR" "TRP" "TYR" "VAL"

scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
if(length(scatter.psi) > 0){
  png(paste("PPIbindersNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
  print(paste("PPIbindersNoGLY.png", sep = ""))
  z <- kde2d(scatter.phi, scatter.psi, n=50)
  tryCatch({
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("PPI binders (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    legend("topright", paste("n=", length(scatter.psi)), bty="n")
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedABDBMotifsPointsOnlyForR.txt", header=TRUE, colClasses=c("character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain", "AA", "Position", "phi","psi"));

scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
if(length(scatter.psi) > 0){
  png(paste("ABDBbindersNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
  print(paste("ABDBbindersNoGLY.png", sep = ""))
  z <- kde2d(scatter.phi, scatter.psi, n=50)
  tryCatch({
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB binders (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    legend("topright", paste("n=", length(scatter.psi)), bty="n")
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedABDBALLMotifsOnlyForR.txt", header=TRUE, colClasses=c("character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain", "AA", "Position", "phi","psi"));

scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
if(length(scatter.psi) > 0){
  png(paste("ABDBMotifsNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
  print(paste("ABDBMotifsNoGLY.png", sep = ""))
  z <- kde2d(scatter.phi, scatter.psi, n=50)
  tryCatch({
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB Motifs (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    legend("topright", paste("n=", length(scatter.psi)), bty="n")
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedAllFolderPPIcutForR.txt", header=TRUE, colClasses=c("character","factor","numeric","character","numeric","numeric"), col.names=c("ABDBname", "AA", "Position", "chain", "phi","psi"));

scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
if(length(scatter.psi) > 0){
  png(paste("PPImotifsNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
  print(paste("PPImotifsNoGLY.png", sep = ""))
  z <- kde2d(scatter.phi, scatter.psi, n=50)
  tryCatch({
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("PPI motifs (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    legend("topright", paste("n=", length(scatter.psi)), bty="n")
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/analyzedTabs.txt", header=FALSE, colClasses=c("character","character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain","CDR", "AA", "Position", "phi","psi"));


CDRvec = c("CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3")
scatter.psi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
if(length(scatter.psi) > 0){
  png(paste("CDRsNoGly.png", sep = ""), width = 8000, height = 7000, res = 1000)
  print(paste("CDRsNoGly.png", sep = ""))
  z <- kde2d(scatter.phi, scatter.psi, n=50)
  tryCatch({
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB CDR residues (Non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    legend("topright", paste("n=", length(scatter.psi)), bty="n")
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

FRvec = c( "HFR1" ,  "HFR2" ,  "HFR3" ,  "HFR4" ,  "LFR1" ,  "LFR2" ,  "LFR3" ,  "LFR4" )
scatter.psi <- scatter.data[which(scatter.data[,"CDR"] %in% FRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"CDR"] %in% FRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
if(length(scatter.psi) > 0){
  png(paste("FRsNoGly.png", sep = ""), width = 8000, height = 7000, res = 1000)
  print(paste("FRsNoGly.png", sep = ""))
  z <- kde2d(scatter.phi, scatter.psi, n=50)
  tryCatch({
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB FR residues (Non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    legend("topright", paste("n=", length(scatter.psi)), bty="n")
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



