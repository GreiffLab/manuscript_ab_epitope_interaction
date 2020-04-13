# network aa physchem
setwd('/Users/rahmadakbar/greifflab/aims/aimugen/src')
# detach('package:themeakbar', unload = TRUE)
library('themeakbar')
library('janitor')
library('ggplot2')
library('pheatmap')
library('beeswarm')
library(data.table)
library('dplyr')
library('rlist')
library('scales')
library('reshape2')
library('VennDiagram')
library('tidyr')
library('RCy3')
library(RColorBrewer)
library('igraph')
library('tidyverse')
library('viridis')
library('poweRlaw')
# library('fossil')
library('gtools')
library('vegan')
library('iNEXT')
library('wesanderson')
library('Peptides')
library('ggfortify')

names(wes_palettes)
theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
my_darjee = colorRampPalette(wes_palette('Darjeeling1'))(8)
print(my_darjee)
abcolor = '#6464FF'
abcolor2 = '#1010FF'
abcolor3 = '#3399CC'
agcolor = '#C0C0C0'
blues = colorRampPalette(brewer.pal(8,'Blues'))(14)[2:10]
oranges = brewer.pal(6,'Oranges')
greys = brewer.pal(6,'Greys')
blues2 = brewer.pal(6,'Blues')
blues3 = c('#9010FF', '#1010FF', '#C8C8FF')
blues4 = c(blues3, blues2[2:6])
cdrfr_colors = c(oranges,blues)
rewcolors = c(greys[3:5], c('#9010FF', '#1010FF', '#C8C8FF'))
rewcolors2  = c(greys[3:4], blues3[1:2])
blues5 = c('#1010FF','#1010FF','#1010FF',greys[3],'#1010FF')
outfigdir = 'abdb_figures_2019'
ragcolor_rgb = c(col2rgb('lightsalmon'))/255
ragcolor = rgb(ragcolor_rgb[1], ragcolor_rgb[2], ragcolor_rgb[3])
rabcolor_rgb = c(col2rgb('lightpink'))/255
rabcolor = rgb(rabcolor_rgb[1], rabcolor_rgb[2], rabcolor_rgb[3])

ppi_mcolor_rgb = c(col2rgb('orange'))/255
ppi_mcolor = rgb(ppi_mcolor_rgb[1], ppi_mcolor_rgb[2], ppi_mcolor_rgb[3])
ppi_pcolor_rgb = c(col2rgb('darkorange'))/255
ppi_pcolor = rgb(ppi_pcolor_rgb[1], ppi_pcolor_rgb[2], ppi_pcolor_rgb[3])

print(ppi_mcolor)
print(ppi_pcolor)

my_ccols = c('salmon',
             'deepskyblue1', 'deepskyblue2', 'deepskyblue3', 'deepskyblue4',
             'orange1', 'orange2', 'orange3', 'orange4')


my_ccols2 = c(abcolor3, 'lightpink',
              abcolor3, 'lightpink',
              abcolor3, 'lightpink',
              abcolor3, 'lightpink')


my_ccols3 = c(ggplot2::alpha('lightpink', 1), ggplot2::alpha(ppi_mcolor,0.5),
              ggplot2::alpha('lightpink', 1), ggplot2::alpha(ppi_mcolor,0.5), 
              ggplot2::alpha('lightpink',1), ggplot2::alpha(ppi_mcolor,0.5), 
              ggplot2::alpha('red', 0.5), ggplot2::alpha(ppi_mcolor, 1))


outpdf = function(infile, tag, width=8, height=8){
  # assumes infile has extensins (.csv etc)
  # uses the first bit of the file as name and add the tag
  # uses ggsave, opens pdf plot after saving
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_%s.pdf', outfigdir,inname,tag)
  ggsave(outname, width=width, height = height)
  system(sprintf('open %s', outname))
  print(sprintf('opening %s', outname))
  
}

outpng = function(infile, tag, width=8, height=8){
  # assumes infile has extensins (.csv etc)
  # uses the first bit of the file as name and add the tag
  # uses ggsave, opens pdf plot after saving
  inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
  outname = sprintf('%s/%s_%s.png', outfigdir,inname,tag)
  ggsave(outname, width=width, height = height)
  system(sprintf('open %s', outname))
  print(sprintf('opening %s', outname))
}



#install.packages("ggplot2")


# ### from P. Robert
# 
# library("ggplot2")
# library("tibble")
# library("viridis")
# AA = "LEU"
# chain = "H"
# CDR = "CDR-H2"
# 
# ## some pretty colors
# library(MASS)  # in case it is not already loaded 
# library(RColorBrewer)
# k <- 11
# my.cols <- rev(brewer.pal(k, "RdYlBu"))
# 
# setwd("C:/Users/pprobert/Desktop/TempWorld/");
# 
# scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedMotifsPointsOnlyForR.txt", header=TRUE, colClasses=c("character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain", "AA", "Position", "phi","psi"));
# 
# listAAs = levels(factor(scatter.data$"AA"))	#"ALA" "ARG" "ASN" "ASP" "CYS" "GLN" "GLU" "GLY" "HIS" "ILE" "LEU"  "LYS" "MET" "PHE" "PRO" "SER" "THR" "TRP" "TYR" "VAL"
# 
# scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
# scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
# if(length(scatter.psi) > 0){
#   png(paste("PPIbindersNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
#   print(paste("PPIbindersNoGLY.png", sep = ""))
#   z <- kde2d(scatter.phi, scatter.psi, n=50)
#   tryCatch({
#     plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("PPI binders (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
#     contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
#     legend("topright", paste("n=", length(scatter.psi)), bty="n")
#     dev.off()
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# 
# scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedABDBMotifsPointsOnlyForR.txt", header=TRUE, colClasses=c("character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain", "AA", "Position", "phi","psi"));
# 
# scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
# scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
# if(length(scatter.psi) > 0){
#   png(paste("ABDBbindersNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
#   print(paste("ABDBbindersNoGLY.png", sep = ""))
#   z <- kde2d(scatter.phi, scatter.psi, n=50)
#   tryCatch({
#     plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB binders (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
#     contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
#     legend("topright", paste("n=", length(scatter.psi)), bty="n")
#     dev.off()
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# 
# scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedABDBALLMotifsOnlyForR.txt", header=TRUE, colClasses=c("character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain", "AA", "Position", "phi","psi"));
# 
# scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
# scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
# if(length(scatter.psi) > 0){
#   png(paste("ABDBMotifsNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
#   print(paste("ABDBMotifsNoGLY.png", sep = ""))
#   z <- kde2d(scatter.phi, scatter.psi, n=50)
#   tryCatch({
#     plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB Motifs (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
#     contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
#     legend("topright", paste("n=", length(scatter.psi)), bty="n")
#     dev.off()
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# 
# scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/AnalyzedAllFolderPPIcutForR.txt", header=TRUE, colClasses=c("character","factor","numeric","character","numeric","numeric"), col.names=c("ABDBname", "AA", "Position", "chain", "phi","psi"));
# 
# scatter.psi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
# scatter.phi <- scatter.data[which(scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
# if(length(scatter.psi) > 0){
#   png(paste("PPImotifsNoGLY.png", sep = ""), width = 4000, height = 3500, res = 500)
#   print(paste("PPImotifsNoGLY.png", sep = ""))
#   z <- kde2d(scatter.phi, scatter.psi, n=50)
#   tryCatch({
#     plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("PPI motifs (non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
#     contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
#     legend("topright", paste("n=", length(scatter.psi)), bty="n")
#     dev.off()
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# 
# 
# scatter.data <- read.table("C:/Users/pprobert/Desktop/Main/D-ABDB_Paper/Ramachandran/PhiPsiPlots/analyzedTabs.txt", header=FALSE, colClasses=c("character","character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain","CDR", "AA", "Position", "phi","psi"));
# 
# 
# CDRvec = c("CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3")
# scatter.psi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
# scatter.phi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
# if(length(scatter.psi) > 0){
#   png(paste("CDRsNoGly.png", sep = ""), width = 8000, height = 7000, res = 1000)
#   print(paste("CDRsNoGly.png", sep = ""))
#   z <- kde2d(scatter.phi, scatter.psi, n=50)
#   tryCatch({
#     plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB CDR residues (Non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
#     contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
#     legend("topright", paste("n=", length(scatter.psi)), bty="n")
#     dev.off()
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# FRvec = c( "HFR1" ,  "HFR2" ,  "HFR3" ,  "HFR4" ,  "LFR1" ,  "LFR2" ,  "LFR3" ,  "LFR4" )
# scatter.psi <- scatter.data[which(scatter.data[,"CDR"] %in% FRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
# scatter.phi <- scatter.data[which(scatter.data[,"CDR"] %in% FRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
# if(length(scatter.psi) > 0){
#   png(paste("FRsNoGly.png", sep = ""), width = 8000, height = 7000, res = 1000)
#   print(paste("FRsNoGly.png", sep = ""))
#   z <- kde2d(scatter.phi, scatter.psi, n=50)
#   tryCatch({
#     plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("ABDB FR residues (Non Glycine)"), pch=19, cex=.4, cex.lab=2.0, cex.axis=2.0, cex.main = 2.0)
#     contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
#     legend("topright", paste("n=", length(scatter.psi)), bty="n")
#     dev.off()
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# ### end from P. Robert

rama_clusters  = function(){
  infile = 'abdb_outfiles_2019/RamachandranClusters.csv'
  df = read_csv(infile)
  print(df)
  cnames = colnames(df)
  outdf = data_frame()
  n = dim(df)[1]
  for (cname in cnames[2:length(cnames)]){
    percent_content =  df[,cname]
    colnames(percent_content) = 'pcnt'
    print(percent_content)
    tempdf = data.frame(percent = percent_content, 'source' = rep(cname, n), 'rama_region' = df$rama_region )
    outdf = rbind(outdf, tempdf)
    print(tempdf)
    # print(outdf)
  }
  outdf$pcnt2 = as.numeric(gsub(",", ".", gsub("\\.", "", outdf$pcnt)))
  print(outdf)
  # stop()
  ggplot(data=outdf) + 
    geom_bar(mapping = aes(x=rama_region, y = pcnt2, fill= source), stat = 'identity', position = 'dodge') + 
    scale_fill_manual(values = c(abcolor3, ppi_mcolor, alpha(abcolor3,0.5), abcolor, abcolor2, ppi_pcolor)) + 
    labs(x = 'Ramachandran regions', y = 'Percentage') + 
    theme(axis.title = element_text(size=25),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          axis.text = element_text(size = 20))
  
  outpdf(infile, 'bar', width = 9)
}

# run stuff
rama_clusters()



