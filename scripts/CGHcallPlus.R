#this script contain the functions: make_cghRawPlus, frequencyPlot, segmentDataWeighted,  CGHregionsPlus,
#regioningPlus, repdata , WECCA.heatmapPlus, mark.genes , mark.bed,  add.cytobands , add.genes, plot.profiles

suppressMessages(library(CGHcall))
suppressMessages(library(CGHregions))
suppressMessages(library(WECCA))
suppressMessages(library(matrixStats))
suppressMessages(library(QDNAseq))

# originally: QDNAseqReadCounts instead of QDNAseqSignals
setMethod('plot', signature(x='cghRaw', y='missing'),
  getMethod('plot', signature=c(x='QDNAseqSignals', y='missing')))
setMethod('plot', signature(x='cghSeg', y='missing'),
  getMethod('plot', signature=c(x='QDNAseqSignals', y='missing')))
setMethod('plot', signature(x='cghCall', y='missing'),
  getMethod('plot', signature=c(x='QDNAseqSignals', y='missing')))

setMethod("frequencyPlot", signature(x="cghCall", y="missing"), frequencyPlotCalls)

.CGHcallPlus <- new.env()
evalq({

setMethod("frequencyPlot", signature(x="cghRegions", y="missing"),
function (x, y, main='Frequency Plot', gaincol='red', losscol='blue', misscol=NA, build='GRCh37',... ) #TLos changed colors loss bleu, gain red
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  chrom.lengths <- CGHbase:::.getChromosomeLengths(build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  calls <- regions(x)
  loss.freq <- rowMeans(calls < 0)
  gain.freq <- rowMeans(calls > 0)
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='frequency', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  if (length(chrom.ends) > 1)
    for (j in names(chrom.ends)[-length(chrom.ends)])
      abline(v=chrom.ends[j], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  mtext(paste(nrow(x), 'regions'), side=3, line=0, adj=0)
})

make_cghRawPlus <-
function (input)
{
    if (class(input) == "character")
        input <- read.table(input, header = T, sep = "\t", fill = T,
            quote = "")
    if (class(input[, 2]) == "factor")
        input[, 2] <- as.character(input[, 2])
    if (class(input[, 2]) == "character") {
        input[, 2] <- sub("^chr", "", input[, 2])
        input[input[, 2] == "X", 2] <- "23"
        input[input[, 2] == "Y", 2] <- "24"
        input[input[, 2] == "MT", 2] <- "25"
        input[, 2] <- as.integer(input[, 2])
    }
    if (any(duplicated(input[, 1]))) {
        replicate.probes <- unique(input[, 1][duplicated(input[, 1])])
        uniques <- input[!input[, 1] %in% replicate.probes, ]
        replicates <- input[input[, 1] %in% replicate.probes, ]
        replicates.avg <- aggregate(replicates[, -1],
            list(probe=replicates[, 2]), median, na.rm=TRUE)
        input <- rbind(uniques, replicates.avg)
    }
    input <- input[order(input[, 2], input[, 3]), ]
    copynumber <- as.matrix(input[, 5:ncol(input)])
    rownames(copynumber) <- input[, 1]
    if (ncol(copynumber) == 1)
        colnames(copynumber) <- colnames(input)[5]
    annotation <- data.frame(Chromosome = input[, 2], Start = input[,
        3], End = input[, 4], row.names = input[, 1])
    metadata <- data.frame(labelDescription = c("Chromosomal position",
        "Basepair position start", "Basepair position end"),
        row.names = c("Chromosome", "Start", "End"))
    dimLabels <- c("featureNames", "featureColumns")
    annotation <- new("AnnotatedDataFrame", data = annotation,
        dimLabels = dimLabels, varMetadata = metadata)
    result <- new("cghRaw", copynumber = copynumber, featureData = annotation)
}
environment(make_cghRawPlus) <- environment(CGHbase:::make_cghRaw)
make_cghRaw <- make_cghRawPlus

segmentDataWeighted <-
function (input, method = "DNAcopy", ...)
{
    if (method == "DNAcopy") {
        CNA.object <- DNAcopy::CNA(copynumber(input), chromosomes(input),
            bpstart(input), data.type = "logratio")
        cat("Start data segmentation .. \n")
        segmented <- segment(CNA.object, ...)
        numclone <- segmented$output$num.mark
        smrat <- segmented$output$seg
        numsmrat <- cbind(smrat, numclone)
        repdata <- function(row) {
            rep(row[1], row[2])
        }
        makelist <- apply(numsmrat, 1, repdata)
        joined <- unlist(makelist)
        rm(makelist)
        joined <- matrix(joined, ncol = ncol(input), byrow = FALSE)
        joined <- CGHcall:::.assignNames(joined, input)
        result <- CGHcall:::.segFromRaw(input, joined)
    }
    result
}

CGHregionsPlus <- function(input, ...) {
  regions <- CGHregions:::CGHregions(input, ...)
  # End positions of regions should be the end position of the last data point of that region,
  # but instead CGHregions returns the start position of the last data point.
  # Check if that is indeed the case:
  if (class(input) == 'cghCall') {
    if (sum(regions@featureData@data$End %in% input@featureData@data$Start) > sum(regions@featureData@data$End %in% input@featureData@data$End))
      for (row in rownames(regions@featureData@data))
        regions@featureData@data[row, 'End'] <- input@featureData@data[input@featureData@data$Chromosome == regions@featureData@data[row, 'Chromosome'] & input@featureData@data$Start == regions@featureData@data[row, 'End'], 'End'][1]
  }
  regions
}
environment(CGHregionsPlus) <- environment(CGHregions:::CGHregions)
CGHregions <- CGHregionsPlus

regioningPlus <- function (cghdata.called, threshold = 0.00001, cghdata.regions = NULL)
{
    find.reg.modus <- function(x) {
        if (nrow(x) == 1)
            return(x)
        splitter <- list()
        splitter[[1]] <- c(1)
        index.temp <- 1
        j <- 1
        for (i in 1:(dim(x)[1] - 1)) {
            if (all(x[i, ] == x[i + 1, ])) {
                index.temp <- c(index.temp, i + 1)
                splitter[[j]] <- index.temp
            }
            else {
                index.temp <- i + 1
                j <- j + 1
                splitter[[j]] <- index.temp
            }
        }
        region.details <- NULL
        for (i in 1:length(splitter)) {
            region.details <- rbind(region.details, c(min(splitter[[i]]),
                max(splitter[[i]])))
        }
        modus <- which.max(region.details[, 2] - region.details[,
            1] + 1)
        return(x[region.details[modus[1], 1], ])
    }
    cat("CGHregions of hard call data...")
    if (is.null(cghdata.regions))
        cghdata.regions <- CGHregionsPlus(cghdata.called, averror = threshold)
    cat("...done", "\n")
    print(paste("threshold used:", threshold, sep = " "))
    calls.annotation <- pData(featureData(cghdata.called))
    regions.annotation <- pData(featureData(cghdata.regions))
    cat("Map regions to clones...")
    reg.to.clones <- list()
    counter <- 0
    for (chr in unique(regions.annotation[, 1])) {
        reg.ann.temp <- regions.annotation[regions.annotation[,
            1] == chr, 1:4]
        for (r in 1:dim(reg.ann.temp)[1]) {
            counter <- counter + 1
            A1 <- which(calls.annotation[, 1] == chr)
            A2 <- which(calls.annotation[, 2] >= reg.ann.temp[r,
                2])
            A3 <- which(calls.annotation[, 3] <= reg.ann.temp[r,
                3])
            reg.to.clones[[counter]] <- intersect(intersect(A1,
                A2), A3)
        }
    }
    cat("...done", "\n")
    cghdata.probs <- numeric()
    for (i in 1:dim(calls(cghdata.called))[2]) {
        cghdata.probs <- cbind(cghdata.probs, cbind(probloss(cghdata.called)[,
            i], probnorm(cghdata.called)[, i], probgain(cghdata.called)[,
            i], probamp(cghdata.called)[, i]))
    }
    cat("Calculate mode soft call signature for each region...")
    cghdata.regprobs <- numeric()
    for (i in 1:length(reg.to.clones)) {
        cghdata.regprobs <- rbind(cghdata.regprobs, find.reg.modus(cghdata.probs[reg.to.clones[[i]],
            , drop = FALSE]))
    }
    cat("...done", "\n")
    softcalls.samplenames <- character()
    for (i in 1:dim(calls(cghdata.called))[2]) {
        if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] ==
            3) {
            softcalls.samplenames <- c(softcalls.samplenames,
                paste(c("probloss_", "probnorm_", "probgain_"),
                  colnames(regions(cghdata.regions))[i], sep = ""))
        }
        if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] ==
            4) {
            softcalls.samplenames <- c(softcalls.samplenames,
                paste(c("probloss_", "probnorm_", "probgain_",
                  "probamp_"), colnames(regions(cghdata.regions))[i],
                  sep = ""))
        }
    }
    colnames(cghdata.regprobs) <- softcalls.samplenames
    rownames(cghdata.regprobs) <- rownames(regions(cghdata.regions))
    regdata <- list()
    regdata$ann <- regions.annotation
    regdata$hardcalls <- regions(cghdata.regions)
    regdata$softcalls <- cghdata.regprobs
    return(regdata)
}
environment(regioningPlus) <- environment(WECCA:::regioning)
regioning <- regioningPlus

WECCA.heatmapPlus <- function (cghdata.regioned, dendrogram, build='GRCh37',
  ...) {
  nclasses <- sort(unique(as.numeric(cghdata.regioned$hardcalls)))
  cols <- c('lightgreen', 'darkgreen', 'lightgray', 'darkslategray')
  chr.color <- rep(1, nrow(cghdata.regioned$ann))
  centromeres <- CGHbase:::.getCentromere(build)
  for (chr in unique(cghdata.regioned$ann$Chromosome))
    chr.color[cghdata.regioned$ann$Chromosome == chr &
      (cghdata.regioned$ann$Start + cghdata.regioned$ann$End) / 2 >
      centromeres[chr]] <- 2
  even <- cghdata.regioned$ann$Chromosome %% 2 == 0
  chr.color[even] <- chr.color[even] + 2
  chr.color <- cols[chr.color]

  Y <- rep(FALSE, dim(cghdata.regioned$hardcalls)[1])
  for (i in 2:(dim(cghdata.regioned$ann)[1])) {
      if ((cghdata.regioned$ann[i - 1, 1] != cghdata.regioned$ann[i,
          1])) {
          Y[i] <- TRUE
      }
  }
  Y[1] <- TRUE
  begin.chr <- rep("", dim(cghdata.regioned$ann)[1])
  begin.chr[Y] <- cghdata.regioned$ann[Y, 1]
  color.coding <- c("-2"="darkblue", "-1"="blue", "0"="black", "1"="red",
    "2"="darkred")[as.character(nclasses)]  #TLos changed colors loss bleu, gain red.
  heatmap(cghdata.regioned$hardcalls, Colv = as.dendrogram(dendrogram),
      Rowv=NA, col=color.coding, labRow=begin.chr, RowSideColors=chr.color,
      scale="none", ...)
}
environment(WECCA.heatmapPlus) <- environment(WECCA:::WECCA.heatmap)
WECCA.heatmap <- WECCA.heatmapPlus

mark.genes <- function(symbols, chrs=1:24, build='GRCh37', side=3, line=-1, ...) {
  genes <- AnnotationDbi::select(Homo.sapiens::Homo.sapiens, keys=keys(Homo.sapiens::Homo.sapiens,
    keytype='SYMBOL'), cols=c('CHRLOC', 'CHRLOCEND'), keytype='SYMBOL')
  genes$CHRLOCCHR[genes$CHRLOCCHR == 'X'] <- '23'
  genes$CHRLOCCHR[genes$CHRLOCCHR == 'Y'] <- '24'
  genes$CHRLOCCHR[genes$CHRLOCCHR == 'MT'] <- '25'
  genes <- genes[genes$CHRLOCCHR %in% as.character(1:25),]
  genes$CHRLOCCHR <- as.integer(genes$CHRLOCCHR)
  if (length(side) == 1)
    side <- rep(side, length(symbols))
  if (length(line) == 1)
    line <- rep(line, length(symbols))

  chrom           <- genes$CHRLOCCHR
  pos             <- abs(genes$CHRLOC)
  pos2            <- abs(genes$CHRLOCEND)
  uni.chrom <- chrs
  chrom.lengths <- CGHbase:::.getChromosomeLengths(build)[as.character(uni.chrom)]
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
  }
  genes$pos <- pos
  genes$pos2 <- pos2
  for (i in 1:length(symbols)) {
    matches <- genes[genes$SYMBOL == symbols[i],]
    rect(matches$pos, par('usr')[3], matches$pos2, par('usr')[4], col='#88888888', border='#88888888')
    axis(side=side[i], at=matches$pos+(matches$pos2-matches$pos)/2, labels=rep(symbols[i], nrow(matches)), tick=FALSE, line=line[i], cex.axis=.75, ...)
  }
}

mark.bed <- function(bed, chrs=1:24, build='GRCh37', col='#88888888') {
  if (class(bed) == 'character')
    bed <- read.table(bed, sep='\t', as.is=TRUE)
  colnames(bed) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand')
  bed$chromosome <- sub('^chr', '', bed$chromosome)
  bed$chromosome[bed$chromosome == 'X'] <- '23'
  bed$chromosome[bed$chromosome == 'Y'] <- '24'
  bed$chromosome[bed$chromosome %in% c('M', 'MT')] <- '25'
  bed$chromosome <- as.integer(bed$chromosome)
  # bed$start <- bed$start + 1

  chrom           <- bed$chromosome
  pos             <- bed$start
  pos2            <- bed$end
  uni.chrom <- chrs
  chrom.lengths <- CGHbase:::.getChromosomeLengths(build)[as.character(uni.chrom)]
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
  }
  bed$pos <- pos
  bed$pos2 <- pos2
  rect(bed$pos, -5.3, bed$pos2, 5.3, col=col, border=col)
}

add.cytobands <- function(dat, genome.build = 'GRCh37') {
  bands <- read.table(paste('http://www.cangem.org/download.php?platform=CG-PLM-6&flag=', genome.build, sep=''), sep='\t', header=TRUE, as.is=TRUE)
  colnames(bands) <- tolower(colnames(bands))
  colnames(bands)[colnames(bands)=='chr'] <- 'chromosome'
  rownames(bands) <- bands[,1]
  bands$chromosome[bands$chromosome=='X'] <- '23'
  bands$chromosome[bands$chromosome=='Y'] <- '24'
  bands$chromosome[bands$chromosome=='MT'] <- '25'
  bands$chromosome <- as.integer(bands$chromosome)
  dat$cytoband <- NA
  tmp <- colnames(dat)
  colnames(dat) <- tolower(colnames(dat))
  for (band in rownames(bands)) {
    index <- !is.na(dat$chromosome) & dat$chromosome == bands[band, 'chromosome'] &
                  !is.na(dat$start) & dat$start >= bands[band, 'start'] &
                  !is.na(dat$start) & dat$start <= bands[band, 'end']
    if (length(index)>0)
      dat[index, 'startband'] <- bands[band, 'band']
    index <- !is.na(dat$chromosome) & dat$chromosome == bands[band, 'chromosome'] &
                    !is.na(dat$end) & dat$end >= bands[band, 'start'] &
                    !is.na(dat$end) & dat$end <= bands[band, 'end']
    if (length(index)>0)
      dat[index, 'endband'] <- bands[band, 'band']
  }
  dat$startband[is.na(dat$startband)] <- 'unknown'
  dat$endband[is.na(dat$endband)] <- 'unknown'
  dat$cytoband <- paste(dat$startband, '-', dat$endband, sep='')
  dat$cytoband[dat$startband==dat$endband] <- dat$startband[dat$startband==dat$endband]
  dat$startband <- NULL
  dat$endband <- NULL
  colnames(dat) <- tmp
  dat
}

add.genes <- function(dat, genome.build = 'GRCh37') {
  if (!exists('genes'))
    genes <<- read.table(paste('http://www.cangem.org/download.php?platform=CG-PLM-26&flag=', genome.build, sep=''), sep='\t', header=TRUE, row.names=1, as.is=TRUE)
  # genes <- genes[-grep('pseudo', genes$type),]
  genes <- genes[genes$type %in% c('protein_coding', 'miRNA'),]
  genes$chromosome[genes$chromosome=='X'] <- '23'
  genes$chromosome[genes$chromosome=='Y'] <- '24'
  genes$chromosome[genes$chromosome=='MT'] <- '25'
  genes$chromosome <- as.integer(genes$chromosome)
  colnames(dat) <- tolower(colnames(dat))
  for (region in rownames(dat)) {
    index <- genes$chromosome == dat[region, 'chromosome'] &
                    genes$end  > dat[region, 'start'] &
                  genes$start  < dat[region, 'end']
    dat[region, 'genes'] <- sum(index)
    dat[region, 'symbols'] <- paste(genes[index, 'symbol'], collapse=';')
  }
  dat
}

plot.profiles <- function(cgh, directory, byChr=FALSE) {
  tmp <- sampleNames(cgh)
  if (!file.exists(directory))
    dir.create(directory)
  if ('filter' %in% colnames(fData(cgh))) {
    chrs <- unique(chromosomes(cgh)[fData(cgh)$filter])
  } else {
    chrs <- unique(chromosomes(cgh))
  }
  for (i in 1:length(sampleNames(cgh))) {
    if (byChr) {
      for (chr in chrs) {
        png(file.path(directory, paste(tmp[i], '-chr', chr, '.png', sep='')), width=297, height=210, units='mm', res=150)
        plot(cgh[chromosomes(cgh) == chr,i], ylab=expression(normalized~log[2]~read~count), dotres=1)
        dev.off()
      }
    } else {
      png(file.path(directory, paste(tmp[i], '.png', sep='')), width=297, height=210, units='mm', res=150)
      plot(cgh[,i], ylab=expression(normalized~log[2]~read~count), dotres=1)
      dev.off()
    }
  }
}

}, envir=.CGHcallPlus)
attach(.CGHcallPlus)

# EOF
