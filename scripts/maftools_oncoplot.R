if(exists('snakemake')) {
        summary_math_name =  snakemake@params[["summary_name"]]
        path_math =          snakemake@params[["path_maf"]]
        summaryplot_name =   snakemake@output[["summaryPlot_name"]]
        oncoplot_name =      snakemake@output[["oncoplot_name"]]
} else{
	stop("currently only support Snakemake")
}


library(maftools)


file_list <- list.files(path = path_math, pattern = ".maf", full.names = TRUE)
#file_list <- list.files(path = ".", pattern = ".maf")

merged <- merge_mafs(file_list)

write.mafSummary(maf = merged, basename = summary_math_name)
#write.mafSummary(maf = merged, basename = "merged_all")

variants <- merged@variants.per.sample
mutations <- merged@data

#plotmafSummary(maf = tnbc, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

pdf(summaryplot_name)
plotmafSummary(maf = merged, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


pdf(oncoplot_name)
oncoplot(merged, top = 20, showTumorSampleBarcodes = F)
dev.off()


#TEST
save(merged, file = 'merged.txt')


