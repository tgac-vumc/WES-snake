#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis to make CGH regions
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
suppressMessages(library(CGHcall))
suppressMessages(library(CGHregions))
suppressMessages(library(CGHtest))

source('scripts/CGHcallPlus.R')

recalled <- snakemake@input[["recalled"]]
RegionsCGH<-snakemake@output[["RegionsCGH"]]
profiles <- snakemake@output[["profiles"]]

averr <- snakemake@params[["averr"]]  #0.0075	 # default = 0.01

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T, split=FALSE)
##############################################################################################################
# CGH regions & frequency plot regions
##############################################################################################################

QCN.reCalled <- readRDS(recalled)

# Make CGH
CGH_reCalledRCs <- makeCgh(QCN.reCalled)

# perform CGHregions
reCalledRegionsCGH <- CGHregions(CGH_reCalledRCs, averror=averr)

# save data:
saveRDS(reCalledRegionsCGH, RegionsCGH)

# frequency plot reCalled (CGH) regions
#filename <- paste('profiles/freqPlot/freqPlotREGIONS_', bin , 'kbp.png', sep='')
png(profiles, res=300, width=14, height=7, unit='in')
frequencyPlot(reCalledRegionsCGH)
dev.off()
