#!/usr/bin/env Rscript
library(rtracklayer)
library(csaw)
library(Rsamtools)
library(edgeR)
library(mosaics)

#  for testing:
chip_bam <- "/home/petr/PycharmProjects/cenh3_chip_seq_pipeline/output/mapped_reads/chip.all.sorted.bam"
cb <- BamFile(chip_bam, index="/home/petr/PycharmProjects/cenh3_chip_seq_pipeline/output/mapped_reads/chip.all.sorted.bam.csi")

input_bam <- "/home/petr/PycharmProjects/cenh3_chip_seq_pipeline/output/mapped_reads/input.all.sorted.bam"
ib <- BamFile(input_bam, index="/home/petr/PycharmProjects/cenh3_chip_seq_pipeline/output/mapped_reads/input.all.sorted.bam.csi")
genome_chrom_sizes <- "/home/petr/PycharmProjects/cenh3_chip_seq_pipeline/data/CEN6_ver_220406.fasta.chromsizes"
window_size <- 1000
step_size <- 100

gs <- (read.table(genome_chrom_sizes, header = FALSE, stringsAsFactors = FALSE, row.names = 1))
GS <- gs[,1]; names(GS) <- rownames(gs)

chip <- import(chip_bam)
input <- import(input_bam)

get_density <- function(x, tw=1000000){
  cvg <- coverage(x)
  bins <- tileGenome(seqlengths(x), tilewidth = tw)
  d <- binnedAverage(unlist(bins), cvg, "coverage")
  d
}

GS_GR <- GRanges(seqnames = names(GS), ranges = IRanges(start = 1, end = GS))

sliding_windows <- slidingWindows(GS_GR, window_size, step_size)

counts <- windowCounts(list(cb, ib), spacing = 500, width = 2000, ext=141, filter = 2)

counts_norm <- normOffsets(counts)

chip_counts <- assay(counts_norm)[,1]
input_counts <- assay(counts_norm)[,2]

chip_counts <- assay(counts)[,1]
input_counts <- assay(counts)[,2]




chip_pseudo1 <- chip_counts[seq(1, length(chip_counts), by=2)]
chip_pseudo2 <- chip_counts[seq(2, length(chip_counts), by=2)]
input_pseudo1 <- input_counts[seq(1, length(input_counts), by=2)]
input_pseudo2 <- input_counts[seq(2, length(input_counts), by=2)]
pseudo_counts <- cbind(chip_pseudo1, chip_pseudo2, input_pseudo1, input_pseudo2)
colnames(pseudo_counts) <- c("chip_pseudo1", "chip_pseudo2", "input_pseudo1", "input_pseudo2")

# Create a DGEList object
yp <- DGEList(counts=pseudo_counts)

yp <- calcNormFactors(yp)
design <- model.matrix(~factor(c(1, 1, 0, 0))) # Design matrix for pseudoreplicates
yp <- estimateDisp(yp, design)

fit <- glmQLFit(yp, design)
result <- glmQLFTest(fit, coef=2)

enriched_regions <- topTags(result, n=Inf, p.value=0.05)


exp_design <- model.matrix(~ factor(c("ChiP", "Input")))
y <- asDGEList(counts_norm)

y <- estimateDisp(y, exp_design)

fit <- glmQLFit(y, exp_design)

total_chip <- sum(chip_counts)
total_input <- sum(input_counts)

fold_change <- chip_counts / (input_counts + 1)

par(mfrow=c(4,1))
lims=c(100000, 130000)
plot(log10(fold_change), pch=18, cex=0.5, col="#00000050", xlim=lims)
plot(log10(input_counts+1), pch=18, cex=0.5, col="#00000050", xlim=lims)
plot(log10(chip_counts+1), pch=18, cex=0.5, col="#00000050", xlim=lims)
plot(-log10(pvals), pch=18, cex=0.5, col="#00000050", xlim=lims, ylim=c(0,50))


par(mfrow=c(1,1))
plot(log10(input_counts +1 ), log10(chip_counts + 1), pch = 18, cex=0.5, col="#00000020")
abline(0,1, col = 'red')

ma <- chip_counts/(chip_counts + input_counts + 1)

length(fold_change)
hist(fold_change, breaks=1000, main="Histogram of fold changes", xlab="Fold change")

# Vectorized function to calculate p-values using Fisher's exact test
fisher_test_vectorized <- function(chip, input, total_chip, total_input) {
  mapply(function(chip_count, input_count) {
    fisher.test(matrix(c(chip_count, total_chip - chip_count, input_count, total_input - input_count), nrow=2))$p.value
  }, chip, input)
}
p_values <- fisher_test_vectorized(chip_counts, input_counts, total_chip, total_input)


setwd('tmp')

constructBins(chip_bam, fileFormat = 'bam')
constructBins(input_bam, fileFormat = 'bam')


contingency_tables <- mapply(chip_counts, input_counts,FUN=function(chip, input) matrix(c(chip, total_chip-chip, input, total_input-input), nrow = 2), SIMPLIFY = FALSE)
library(parallel)
pvals <- unlist(mclapply(contingency_tables, function(x) fisher.test(x, alternative = 'greater')$p.value, mc.cores = 4))


ch_file <- "chip.all.sorted.bam_fragL200_bin200.txt"
in_file <- "input.all.sorted.bam_fragL200_bin200.txt"
binTFBS <- readBins(type=c('chip', 'input'), fileName = c(ch_file, in_file))

 plot( binTFBS, plotType="input" )


fitTFBS <- mosaicsFit( binTFBS, analysisType="IO", bgEst="rMOM" )

plot(fitTFBS)

 peakTFBS <- mosaicsPeak( fitTFBS, signalModel="1S", FDR=0.05,  maxgap=1000, minsize=5000, thres=5 )
