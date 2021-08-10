#!/usr/bin/env Rscript
#############################################################
#Author: Matias Mendeville
#############################################################

plotQDNAseq <- function(set, path.to.profile){
	writeTo <- path.to.profile

	for (i in 1:length(sampleNames(set))) {
		png(paste(writeTo,sampleNames(set)[i], ".png", sep=""), width=297, height=210, units='mm', res=150)
		plot(set[,i], dotres=1, ylim=c(-3,5), ylab=expression(normalized~log[2]~read~count), main=paste(sampleNames(set)[i], sep=""),
		losscol = "#d1d1ff", gaincol = "#ffd1d1", delcol="#d1d1ff", ampcol="#ffd1d1" )
		dev.off()
	}
}
