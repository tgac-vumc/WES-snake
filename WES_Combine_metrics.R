#!/usr/bin/env R
#########################################################

#Author: Tjitske de Vries
#date:   08-05-2018
#Name:   Combine_metrics.R

#this script is to merge the important output of multiple metrics files.

# The combined metrics file is not a complete overview of all metrics obtained.

 #########################################################

HSmetrics<-snakemake@input[["HSmetrics"]]
Alignment_met<-snakemake@input[["AlignmentMetrics"]]
insert_size_metrics<-snakemake@input[["InsertsizeMetrics"]]
#sv_metrics<-snakemake@input[["sv_metrics"]]
sample=snakemake@wildcards[["sample"]]
output=snakemake@output[["combined"]]



#header_HSmetrics<-read.table(HSmetrics , skip=6, nrows=1, stringsAsFactors=F)
#metrics<-read.table(HSmetrics, skip=7, nrows=1, col.names=header_HSmetrics[1,1:53])
metrics<-read.delim(HSmetrics, skip=6, nrows=1, header=T)

header_alignmetrics<-read.table(Alignment_met , skip=6, nrows=1, stringsAsFactors=F)
#alignmetrics<-read.table(Alignment_met, skip=9, nrows=1, col.names=header_alignmetrics[1,1:22])
alignmetrics<-read.delim(Alignment_met, skip=8, nrows=1, col.names=header_alignmetrics)

#header_insertmetrics<-read.table(insert_size_metrics , skip=6, nrows=1, stringsAsFactors=F)
#insert_metrics<-read.table(insert_size_metrics, skip=7, nrows=1, col.names=header_insertmetrics[1,1:18])
insert_metrics<-read.delim(insert_size_metrics, skip=6, nrows=1, header=T)

#header_svmetrics<-read.table(sv_metrics , skip=1, nrows=1, stringsAsFactors=F)
#svmetrics<-read.table(sv_metrics, skip=3, nrows=1, col.names=header_svmetrics[1,1:14])
#svmetrics<-read.delim(sv_metrics, skip=2)


all<-cbind(metrics,alignmetrics, insert_metrics)
all$SAMPLE<-sample
# this is part of the SV metrics which is only part of gridss
#all$PCT_STRUCTURAL_VAR_READS=all$STRUCTURAL_VARIANT_READS/all$TOTAL_READS
#all$PCT_INDEL_READS=all$INDEL_READS/all$TOTAL_READS
#all$PCT_SPLIT_READS=all$SPLIT_READS/all$TOTAL_READS
#all$PCT_SOFT_CLIPPED_READS=all$SOFT_CLIPPED_READS/all$TOTAL_READS

summary=all[,c("SAMPLE","TOTAL_READS","PF_UNIQUE_READS","PCT_PF_UQ_READS","PF_UQ_READS_ALIGNED", "PCT_PF_UQ_READS_ALIGNED", "PCT_OFF_BAIT",
             "MEAN_BAIT_COVERAGE", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE","PCT_USABLE_BASES_ON_BAIT", "PCT_USABLE_BASES_ON_TARGET", 
            "FOLD_ENRICHMENT", "ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_100X" ,"PCT_CHIMERAS","PCT_ADAPTER","PF_INDEL_RATE", 
            "READS_ALIGNED_IN_PAIRS", "PCT_READS_ALIGNED_IN_PAIRS", "MEAN_READ_LENGTH","MEDIAN_INSERT_SIZE","MEDIAN_ABSOLUTE_DEVIATION", "MEAN_INSERT_SIZE", 
            "STANDARD_DEVIATION", "GC_DROPOUT","AT_DROPOUT","PCT_EXC_DUPE", "PCT_EXC_MAPQ", "PCT_EXC_BASEQ", "PCT_EXC_OVERLAP", "PCT_EXC_OFF_TARGET")]

write.table(summary, file=output, col.names=TRUE, row.names=FALSE, quote=FALSE)



# this is part of the SV metrics which is only part of gridss

#"STRUCTURAL_VARIANT_READS", "PCT_STRUCTURAL_VAR_READS","STRUCTURAL_VARIANT_READ_PAIRS","INDEL_READS","PCT_INDEL_READS","SPLIT_READS","PCT_SPLIT_READS","SOFT_CLIPPED_READS","PCT_SOFT_CLIPPED_READS","UNMAPPED_READS","DISCORDANT_READ_PAIRS","UNMAPPED_MATE_READS","STRUCTURAL_VARIANT_READ_ALIGNMENTS","INDEL_READ_ALIGNMENTS","SPLIT_READ_ALIGNMENTS","SOFT_CLIPPED_READ_ALIGNMENTS","DISCORDANT_READ_PAIR_ALIGNMENTS","UNMAPPED_MATE_READ_ALIGNMENTS",
