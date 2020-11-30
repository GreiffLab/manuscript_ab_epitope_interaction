install.packages("ggplot2")

setwd("C:/Users/pprobert/Desktop/Main/MoveFromKarlskrona/Ramachandran/Ramachandran")
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

scatter.data <- read.table("analyzedTabs.txt", header=FALSE, colClasses=c("character","character","character","factor","numeric","numeric","numeric"), col.names=c("ABDBname","Chain","CDR", "AA", "Position", "phi","psi"));

listAAs = levels(factor(scatter.data$"AA"))	#"ALA" "ARG" "ASN" "ASP" "CYS" "GLN" "GLU" "GLY" "HIS" "ILE" "LEU"  "LYS" "MET" "PHE" "PRO" "SER" "THR" "TRP" "TYR" "VAL"
listchains = levels(factor(scatter.data$"Chain"))
listCDR = levels(factor(scatter.data$"CDR"))

listCDR = c("CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3", "HFR2" ,  "HFR3" ,  "HFR4" ,  "LFR1" ,  "LFR2" ,  "LFR3" ,  "LFR4" )
listAAs = c("ALA", "ARG", "ASN" ,"ASP", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",  "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "VAL")
listAAs = c("CYS")
# Problems with TYR and CYS
# For each chain and heavy/light			
for (AA in listAAs){
	for (chain in listchains){
		for (CDR in listCDR){
			scatter.psi <- scatter.data[which(scatter.data[,"Chain"]==chain & scatter.data[,"CDR"]==CDR & scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
			scatter.phi <- scatter.data[which(scatter.data[,"Chain"]==chain & scatter.data[,"CDR"]==CDR & scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
			if(length(scatter.psi) > 0){
				png(paste(AA ,  chain , CDR , ".png", sep = ""), width = 8000, height = 7000, res = 1000)
				print(paste(AA ,  chain , CDR , ".png", sep = ""))
				z <- kde2d(scatter.phi, scatter.psi, n=50)
				tryCatch({
					plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste(AA, " in ", CDR), pch=19, cex=.4)
					contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
					legend("topright", paste("n=", length(scatter.psi)), bty="n")
					dev.off()
					}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			}
		}
	}
}

# For each AA, CDR vs FR

CDRvec = c("CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3")
scatter.psi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]!='GLY' & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
if(length(scatter.psi) > 0){
  png(paste("CDRsNoGly.png", sep = ""), width = 8000, height = 7000, res = 1000)
  print(paste("CDRsNoGly.png", sep = ""))
  z <- kde2d(scatter.phi, scatter.psi, n=50)
  tryCatch({
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("all CDRs, no GLY"), pch=19, cex=.4)
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
    plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste("all FRs No GLY"), pch=19, cex=.4)
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    legend("topright", paste("n=", length(scatter.psi)), bty="n")
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


 
for (AA in listAAs){
	CDRvec = c("CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3")
	scatter.psi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
	scatter.phi <- scatter.data[which(scatter.data[,"CDR"] %in% CDRvec & scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
	if(length(scatter.psi) > 0){
		png(paste(AA ,  "CDRs.png", sep = ""), width = 8000, height = 7000, res = 1000)
		print(paste(AA ,  "CDRs.png", sep = ""))
		z <- kde2d(scatter.phi, scatter.psi, n=50)
		tryCatch({
			plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste(AA, " in all CDRs"), pch=19, cex=.4)
			contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
			legend("topright", paste("n=", length(scatter.psi)), bty="n")
			dev.off()
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	FRvec = c( "HFR1" ,  "HFR2" ,  "HFR3" ,  "HFR4" ,  "LFR1" ,  "LFR2" ,  "LFR3" ,  "LFR4" )
	scatter.psi <- scatter.data[which(scatter.data[,"CDR"] %in% FRvec & scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
	scatter.phi <- scatter.data[which(scatter.data[,"CDR"] %in% FRvec & scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
	if(length(scatter.psi) > 0){
		png(paste(AA ,  "FRs.png", sep = ""), width = 8000, height = 7000, res = 1000)
		print(paste(AA ,  "FRs.png", sep = ""))
		z <- kde2d(scatter.phi, scatter.psi, n=50)
		tryCatch({
			plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste(AA, " in all FRs"), pch=19, cex=.4)
			contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
			legend("topright", paste("n=", length(scatter.psi)), bty="n")
			dev.off()
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	scatter.psi <- scatter.data[which(scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
	scatter.phi <- scatter.data[which(scatter.data[,"AA"]==AA & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
	if(length(scatter.psi) > 0){
		png(paste(AA ,  "All.png", sep = ""), width = 8000, height = 7000, res = 1000)
		print(paste(AA ,  "All.png", sep = ""))
		z <- kde2d(scatter.phi, scatter.psi, n=50)
		tryCatch({
			plot(scatter.phi, scatter.psi, xlab=expression(phi), xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), ylab=expression(psi),  main= paste(AA, " in all dataset CDR+FR"), pch=19, cex=.4)
			contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
			legend("topright", paste("n=", length(scatter.psi)), bty="n")
			dev.off()
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
}
## https://stats.stackexchange.com/questions/31726/scatterplot-with-contour-heat-overlay
## compute 2D kernel density, see MASS book, pp. 130-131
			par(pty="s")

			#legend("topleft", paste(AA, " in ", CDR), bty="n")
			#abline(h=mean(scatter.phi), v=mean(scatter.psi), lwd=2)




scatter.psi <- scatter.data[which(scatter.data[,"AA"]=="LEU" & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"psi"]
scatter.phi <- scatter.data[which(scatter.data[,"AA"]=="LEU" & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"phi"]
scatter.chain <- scatter.data[which(scatter.data[,"AA"]=="LEU" & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"Chain"]
scatter.CDR <- scatter.data[which(scatter.data[,"AA"]=="LEU" & !is.nan(scatter.data[,"psi"])& !is.nan(scatter.data[,"phi"])),"CDR"]

# basic plot
par(pty="s")
plot(x=scatter.phi, y=scatter.psi, xlim=c(-3.15,3.15), ylim=c(-3.15,3.15), main="General", pch=20, xlab=expression(phi), ylab=expression(psi), cex=0.1, asp=1.0)



#"https://stackoverflow.com/questions/32206623/what-does-level-mean-in-ggplotstat-density2d"
# https://www.sharpsightlabs.com/blog/density-plot-in-r/
# Prepare for 2D plotting: create a class with x,y points

df <- tibble(x_variable = scatter.phi, y_variable = scatter.psi)

ggplot(df, aes(x = x_variable, y = y_variable)) + stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') + geom_bin2d() + theme_bw()



ggplot(df, aes(x = x_variable, y = y_variable)) + stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') +  scale_fill_viridis()  + geom_point(size=0.5)


ggplot(df, aes(x = x_variable, y = y_variable)) +
  stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') +  scale_fill_viridis() + geom_point()

ggplot(df, aes(x = x_variable, y = y_variable)) +
  stat_density2d(aes(fill = ..level..), contour = F, geom = 'tile') +  scale_fill_viridis() + geom_point(size=0.5)


// doesnt work, polygon is out
ggplot(df, aes(x = x_variable, y = y_variable)) +
  stat_density2d(aes(fill = ..level.., alpha=..level..), size=2, bins=5, geom = 'polygon') +  scale_fill_viridis() + geom_point()

colour = factor(scatter.chain)

// thsi is bullshit
geom_polygon(data=df, aes(x=x_variable, y=y_variable, group=scatter.CDR), color="white", fill="grey80") 
install.packages("viridis")

df <- tibble(x_variable = scatter.psi, y_variable = scatter.phi)
base_plot <- ggplot(df, aes(x = x_variable, y = y_variable)) + 
  geom_point()
base_plot + stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)

filled.contour(x=mid.points, y=mid.points, z=grid, levels=c(0,0.002,0.02,1), col=c("#FFFFFF","#FFE8C5","#FFCC7F"), main="Glycine (Symmetric)", asp=1.0, xlab=expression(phi), ylab=expression(psi))