#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis from fastq to called readcounts
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################

suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
suppressMessages(library(R.cache))
setCacheRootPath(path="../.Rcache")

bam <- snakemake@input[["bams"]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
genome <- snakemake@params[["genome"]]
binReadCounts <- snakemake@output[["binReadCounts"]]

##############################################################################################################
# Get bin annotations and bin read counts
##############################################################################################################

bins <- getBinAnnotations(bin, genome=genome)
QRC <- binReadCounts(bins, bamfiles=bam, cache=TRUE)

sub("(_[ACGT]+)?(_S\\d+)?(_L\\d{3})?_R\\d{1}_\\d{3}(\\.f(ast)?q\\.gz)?$", "", sampleNames(QRC)) -> samples

if (sum(duplicated(samples)) > 0) {
        QRC <- poolRuns(QRC, samples=samples, force=TRUE)
}
saveRDS(QRC, binReadCounts)
