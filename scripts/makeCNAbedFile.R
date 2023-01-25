#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis
# Function to create a BED file of all CNAs per sample
## Date: 8 August 2017, December 2017
# Author: Matias Mendeville
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
suppressMessages(library(QDNAseq))
suppressMessages(library(denstrip))
suppressMessages(library(CGHcall))

source("scripts/functions.R")
source("scripts/addCytobands.R")
source("scripts/CGHcallPlus.R")

recalled <- snakemake@input[["recalled"]]
beddir <- snakemake@params[["beddir"]]
cytoband_data<-snakemake@params[["cytobands"]]
max.focal.size.mb<-snakemake@params[["max_focal_size_mb"]]
failed <- snakemake@params[["failed"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=TRUE, split=FALSE)

##############################################################################################################

QCN.reCalled <- readRDS(recalled)

# create BED files

makeCNAbedFile <- function(CalledQDNAseqReadCounts, max.focal.size.mb=3, beddir,cytoband_data){
  for (i in 1:length(sampleNames(CalledQDNAseqReadCounts))){

		#
		sample <- CalledQDNAseqReadCounts[,i]
		sname <- sampleNames(sample)

		# print
		print( paste('Creating BED files of all CNAs in sample: ', sname, sep='') )

		# get all call and segment values of CNAs in the sample
		allCNA.ix <- which(calls(sample)!=0)
		allCNAs <- calls(sample)[ allCNA.ix ]
		seg.allCNA <- segmented(sample)[ allCNA.ix ]

		# coordinates of all CNAs
		reg.allCNAs <- fData(sample)[ allCNA.ix, 1:3]

		# Put together in BED format
		allCNAsPerBin <- cbind(reg.allCNAs, allCNAs, seg.allCNA)
		colnames(allCNAsPerBin)[4] <- 'call'
		colnames(allCNAsPerBin)[5] <- 'segment'

    bedfile<-paste(beddir, sname, '_allCNAsPerBin.bed', sep="")
		# save as BED file; with or without column names..?
		write.table(allCNAsPerBin, bedfile, sep='\t', row.names=F, col.names=T, quote=F)

		# create BED file for all called regions:
		# collapse all adjacent bins based on unique segment values
		coll.calls.bed <- c()
		#all_segments <- uniquecreate.file()(allCNAsPerBin$segment)
    if(length(unique(allCNAsPerBin$segment)) !=0){
			for (j in 1:length(unique(allCNAsPerBin$segment))){
			# first segment

			SEG.val <- as.numeric( unique(allCNAsPerBin$segment)[j] )

			# if the unique segment value occurs more often, then for each one create an entry
			intervals <- seqToIntervals(which(allCNAsPerBin$segment== SEG.val ))

				for (z in 1:nrow(intervals)){

					start.ix 	<- intervals[z,1]
					end.ix 		<- intervals[z,2]

					SEG.start 	<- as.numeric( allCNAsPerBin[start.ix,2] )
					SEG.end 	<- as.numeric( allCNAsPerBin[end.ix,3] )
					SEG.size 	<- as.numeric( allCNAsPerBin[end.ix,3] - allCNAsPerBin[start.ix,2] )

					SEG.chromo 	<- as.numeric( allCNAsPerBin[end.ix,1] )
					SEG.call 	<- as.numeric( allCNAsPerBin[end.ix,4] )

					# check if chromosome is the same; if not the bin spans more chromos and should be separated
					if(allCNAsPerBin[start.ix, 1] != allCNAsPerBin[end.ix, 1]){
					       print(paste(j, 'WARNING: Segment ',  SEG.val, ' spans more than 1 chromosome', sep=''))}

					SEG.coll 	<- cbind(SEG.chromo, SEG.start, SEG.end, SEG.size, SEG.val, SEG.call)

					coll.calls.bed <- rbind(coll.calls.bed, SEG.coll)
				}
			}

			# reorder dataframe
			CNAbed <- coll.calls.bed
			CNAbed <- CNAbed [ order(CNAbed[,1], CNAbed[,2]), ]
			# if there is only 1 CNA in sample, transform the dataframe
			if	(nrow(coll.calls.bed )==1)  {
			CNAbed <- as.data.frame(t(CNAbed))
			}	else {
			CNAbed <- as.data.frame(CNAbed)
			}

			CNAbed[,4] <- round( CNAbed[,4]/1000000, digits = 2)
			CNAbed[,5] <- round( CNAbed[,5], digits=2)


			# add cytoband info
			CNAbed <- addCytobands(CNAbed, cytoband_data)
			CNAbed <- CNAbed[, c(1:4,7,5,6)]
    }else{CNAbed<- data.frame(matrix(ncol = 7, nrow = 0))}
      colnames(CNAbed) <- c('chromo', 'start', 'end', 'cytobands','sizeInMb', 'segment', 'call')
      fileName <- paste(beddir, sname, '_CNAs.bed', sep='')
      write.table(CNAbed, fileName, sep='\t', row.names=F, col.names=T, quote=F)

			# make a separate BED file with only focal CNAs (< 3 Mb)
			focal.ix <- which(CNAbed[,4] <= max.focal.size.mb)
			focalCNAs <- CNAbed[focal.ix, ]

			focalfileName <- paste(beddir, sname, '_focalCNAs.bed', sep='')
			write.table(focalCNAs, focalfileName, sep='\t', row.names=F, col.names=F, quote=F)
	}

}

makeCNAbedFile(QCN.reCalled, max.focal.size.mb, beddir, cytoband_data)

#create output for failed samples - for snakemake compatibility.
failed_samples<-read.table(failed, stringsAsFactors=FALSE, header=TRUE)
if(length(failed_samples[,1]>0)){for(file in failed_samples[,1]){
    file.create(paste(beddir, file, '_allCNAsPerBin.bed', sep=""))
    file.create(paste(beddir, file, '_CNAs.bed', sep=''))
    file.create(paste(beddir, file, '_focalCNAs.bed', sep=''))
}}
