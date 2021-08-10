#!/usr/bin/env Rscript
##############################################################################

#Author: Tjitske logs
#date: Feb 2018

#This script is a small wrapper around ACE to perform the postanalysisloop in the snakemake pipeline
##############################################################################
suppressMessages(library(QDNAseq))
source('scripts/ACE.R')

inputfile <-snakemake@input[["segmented"]]
outputdir<-snakemake@params[["outputdir"]]
failed <- snakemake@params[["failed"]]
fitpickertable<-snakemake@input[["fitpicker"]]
ploidies<-as.integer(snakemake@wildcards[["ploidy"]])

imagetype <- snakemake@config[["ACE"]][["imagetype"]]
trncname<- snakemake@config[["ACE"]][["trncname"]]

copyNumbersSegmented <- readRDS(inputfile)

postanalysisloop(copyNumbersSegmented , modelsfile=fitpickertable, imagetype=imagetype, outputdir=outputdir, trncname=trncname)

failed_samples<-read.table(failed, stringsAsFactors=FALSE, header=TRUE)
if(length(failed_samples[,1]>0)){for(file in failed_samples[,1]){
    file.create(paste(outputdir,"segmentfiles/",file,"_segments.tsv",sep=""))
}}
