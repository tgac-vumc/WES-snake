#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis segmentBins
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
source("scripts/plotQDNAseq.R")

dewaved<- snakemake@input[["dewaved"]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
segmented <- snakemake@output[["segmented"]]
profiles <- snakemake@params[["profiles"]]
failed<- snakemake@params[["failed"]]
min_used_reads<-snakemake@params[["minimal_used_reads"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=TRUE , split=FALSE)

# Adjust segmentation settings based on binsize
if (bin==15) {SDundo=0.75; alph=1e-15}
if (bin==30) {SDundo=0.75; alph=1e-15}
if (bin==100) {SDundo=0.10; alph=1e-20} # default = 1e-20 # 0.01 used for PELLL_FS8_a0.01
if (bin==1000) {SDundo=0.10; alph=1e-20}
# TODO: mogelijk hogere alpha nodig, om meer segmenten te krijgen: 0.01 (1e-2). Deze setting used in tmp_100kbp

# load data
QCN.fcnsd <- readRDS(dewaved)

QCN.fcnsds <- segmentBins(QCN.fcnsd[,QCN.fcnsd$used.reads > min_used_reads ], undo.splits='sdundo', undo.SD=SDundo, alpha=alph, transformFun="sqrt")
QCN.fcnsdsn <- normalizeSegmentedBins(QCN.fcnsds)

saveRDS(QCN.fcnsdsn, segmented)
plotQDNAseq(QCN.fcnsdsn, profiles)

littledata<-colnames(QCN.fcnsd[,QCN.fcnsd$used.reads <= min_used_reads ])
if(length(littledata>0)){for(file in littledata){file.create(paste(profiles, file,".png",sep=""))
}}

write.table(littledata, file=failed )
