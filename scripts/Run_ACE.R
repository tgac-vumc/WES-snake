#!/usr/bin/env Rscript
##############################################################################

#Author: Tjitske logs
#date: Dec 2017

#This script is a small wrapper around ACE to work in the snakemake pipeline
##############################################################################

suppressMessages(library(QDNAseq))
source('scripts/ACE.R')

ploidies<-as.integer(snakemake@wildcards[["ploidy"]])
inputfile <-snakemake@input[["segmented"]]
outputdir<-snakemake@params[["outputdir"]]
failed <- snakemake@params[["failed"]]
log<-snakemake@log[[1]]

imagetype <- snakemake@config[["ACE"]][["imagetype"]]
method<-snakemake@config[["ACE"]][["method"]]
penalty<-as.numeric(snakemake@config[["ACE"]][["penalty"]])
cap<-as.integer(snakemake@config[["ACE"]][["cap"]])
trncname<- snakemake@config[["ACE"]][["trncname"]]
printsummaries<- snakemake@config[["ACE"]][["printsummaries"]]

copyNumbersSegmented <- readRDS(inputfile)

parameters <- data.frame(options = c("ploidies","imagetype","method","penalty","cap","trncname","printsummaries"),
                         values = c(paste0(ploidies,collapse=", "),imagetype,method,penalty,cap,trncname,printsummaries))

write.table(parameters, file=log, quote = FALSE, sep = "\t", na = "", row.names = FALSE)

ploidyplotloop(copyNumbersSegmented ,outputdir , ploidies,imagetype,method,penalty,cap,trncname,printsummaries)

#create output for failed samples - for snakemake compatibility.
failed_samples<-read.table(failed, stringsAsFactors=FALSE, header=TRUE)
if(length(failed_samples[,1]>0)){for(file in failed_samples[,1]){
    file.create(paste(outputdir, ploidies,"N/", file,"/summary_",file,".",imagetype,sep=""))
}}
