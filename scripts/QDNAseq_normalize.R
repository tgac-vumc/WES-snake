#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis to normalize bins and dewave
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))

source("scripts/functions.R")
source("scripts/plotQDNAseq.R")

binReadCounts <- snakemake@input[["binReadCounts"]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
corrected <- snakemake@output[["corrected"]]
profiles <- snakemake@params[["profiles"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T , split=FALSE)
##############################################################################################################
# Correct, Normalize & Dewave raw data
##############################################################################################################
QRC <- readRDS(binReadCounts)

QRC.f <- applyFilters(QRC, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE)
QRC.f <- estimateCorrection(QRC.f)
QCN.fc <- correctBins(QRC.f)
QCN.fcn <- normalizeBins(QCN.fc)
QCN.fcns <- smoothOutlierBins(QCN.fcn)

saveRDS(QCN.fcns, corrected)
plotQDNAseq(QCN.fcns, profiles)
