##########################################################################################
#!/usr/bin/env Rscript

# Add CYTOBAND INFORMATION of regions for HG19 format #
# input BED file: with required column names: 'chromo', 'start', 'end'
# output BED file with cytoband in extra column

# Date: 15 August 2017
# Matias Mendeville

#December2017 adapted to fit in Snakemake pipeline

##########################################################################################
addCytobands <- function(BED, cytoband_data){

	# Cytoband download: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
	cyto <- read.table(cytoband_data, sep='\t', header=F, 	stringsAsFactors=F)
	colnames(cyto) <- c('chr', 'start', 'end', 'cytoband', 'bla')
	# convert 'chr1' to 1; and X into 23 and Y into 24
	new.chr <- c()
	for ( i in 1:nrow(cyto)){
	chr <- strsplit(cyto$chr[i], 'chr')[[1]][2]
	new.chr <- c(new.chr, chr)
	}
	y.ix <- which(new.chr=='Y')
	x.ix <- which(new.chr=='X')
	new.chr[y.ix] <- 24
	new.chr[x.ix] <- 23
	cyto$chr <- as.numeric(new.chr)

	## rename BED file columns:
	colnames(BED) <- c('chromo', 'start', 'end')

	##
	cytobands <- c()

	for ( i in 1:nrow(BED)){
		# get cytoband of region-start
		v1 <- BED$chromo[i] == cyto$chr
		v2 <- BED$start[i] > cyto$start
		v3 <- BED$start[i] < cyto$end
		cyto.start <- cyto$cytoband[which(v1&v2&v3)]

		# get cytoband of region-end
		v4 <- BED$chromo[i] == cyto$chr
		v5 <- BED$end[i] > cyto$start
		v6 <- BED$end[i] < cyto$end
		cyto.end <- cyto$cytoband[which(v4&v5&v6)]

		# add to vector with all cyto starts and ends
		region.cyto <- paste(cyto.start, cyto.end, sep='-')
		cytobands <- c(cytobands, region.cyto)
	}

	BED <- cbind(BED, cytobands)

} # EOF
