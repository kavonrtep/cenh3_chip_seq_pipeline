#!/usr/bin/env Rscript
library(optparse)


# parse command line arguments
option_list <- list(
  make_option(c("-c", "--chip"), action = "store", type = "character",
              help = "chip bam file", default = NA),
  make_option(c("-i", "--input"), action = "store", type = "character",
              help = "input bam file", default = NA),
  make_option(c("-p", "--prefix"), action = "store", type = "character",
              help = "prefix for output files", default = NA),
  make_option(c("-t", "--threads"), action = "store", type = "integer",
              help = "number of threads", default = 10),
  make_option(c("-n", "--normalized_only"), action = "store_true", default = FALSE,
              help = "only calculate normalized data")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, args = commandArgs(TRUE))

# check if chip, input and prefix are provided
if (is.na(opt$chip) || is.na(opt$input) || is.na(opt$prefix)){
  stop("chip, input and prefix are required")
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Rbeast)
  library(Rsamtools)
  library(csaw)
  library(parallel)
  library(edgeR)
})

chip_bam_file <- opt$chip
input_bam_file <- opt$input
prefix <- opt$prefix

# check if directories with prefix exist, if not create them
if (!dir.exists(dirname(prefix))){
  dir.create(dirname(prefix), recursive = TRUE)
}


if (FALSE){
  # for testing
  chip_bam_file <- "/mnt/ceph/454_data/Pisum_pangenome/assemblies/K-2524_2024-06-11/analysis/ChIP-seq_CENH3/P43-6/mapped_reads/chip.all.sorted.bam"
  input_bam_file <- "/mnt/ceph/454_data/Pisum_pangenome/assemblies/K-2524_2024-06-11/analysis/ChIP-seq_CENH3/P43-6/mapped_reads/input.all.sorted.bam"
}


chip_bam_index <- paste0(chip_bam_file, ".bai")
if (!file.exists(chip_bam_index)){
  chip_bam_index <- paste0(chip_bam_file, ".csi")
  if (!file.exists(chip_bam_index)){
    stop("Index file for chip bam not found")
  }
}

input_bam_index <- paste0(input_bam_file, ".bai")
if (!file.exists(input_bam_index)){
  input_bam_index <- paste0(input_bam_file, ".csi")
  if (!file.exists(input_bam_index)){
    stop("Index file for input bam not found")
  }
}
window_size <- 2000
window_size10k <- 10000
ncpu <- opt$threads


chip_BAM <- BamFile(chip_bam_file, index=chip_bam_index)
input_BAM <- BamFile(input_bam_file, index=input_bam_index)

message("Reading BAM files and counting reads for 2k windows")
counts <- windowCounts(list(chip_BAM, input_BAM),
                       spacing = window_size, width = window_size, filter = 1)


message("Reading BAM files and counting reads for 10k windows")
counts10k <- windowCounts(list(chip_BAM, input_BAM),
                          spacing = window_size10k, width = window_size10k, filter = 1)


# normalize the counts using trimmed mean of M-values
message("Normalizing counts")
count_list <- DGEList(assay(counts))
count_list <- calcNormFactors(count_list, method="TMM")
normalized_counts <- cpm(count_list)

counts10k_list <- DGEList(assay(counts10k))
counts10k_list <- calcNormFactors(counts10k_list, method="TMM")
normalized_counts10k <- cpm(counts10k_list)



gr <- rowRanges(counts)
gr$chip <- assay(counts)[,1]
gr$input <- assay(counts)[,2]
gr$log2fc <- log2(gr$chip + 1) - log2(gr$input + 1)
gr$log2fc_weighted <- gr$log2fc * log(gr$chip + gr$input + 1)
gr$chip_norm <- normalized_counts[,1]
gr$input_norm <- normalized_counts[,2]
min_value <- min(normalized_counts[normalized_counts != 0])
gr$log2fc_norm <- log2(gr$chip_norm + min_value) - log2(gr$input_norm + min_value)
gr$log2fc_weighted_norm <- gr$log2fc_norm * log(gr$chip_norm + gr$input_norm + min_value)
gr_part <- split(gr, seqnames(gr))

gr10k <- rowRanges(counts10k)
gr10k$chip <- assay(counts10k)[,1]
gr10k$input <- assay(counts10k)[,2]
gr10k$log2fc <- log2(gr10k$chip + 1) - log2(gr10k$input + 1)
gr10k$log2fc_weighted <- gr10k$log2fc * log(gr10k$chip + gr10k$input + 1)
gr10k$chip_norm <- normalized_counts10k[,1]
gr10k$input_norm <- normalized_counts10k[,2]
min_value <- min(normalized_counts10k[normalized_counts10k != 0])
gr10k$log2fc_norm <- log2(gr10k$chip_norm + min_value) - log2(gr10k$input_norm + min_value)
gr10k$log2fc_weighted_norm <- gr10k$log2fc_norm * log(gr10k$chip_norm + gr10k$input_norm + min_value)
gr10k_part <- split(gr10k, seqnames(gr10k))





# calculate beast on each chromosome or scaffold, min size is 10 intervals
gr10k_part_filtered <- gr10k_part[sapply(gr10k_part, function(x)length(x$log2fc_weighted) > 10)]
if (!opt$normalized_only){
  message("Calculating BEAST on 10k windows")
  beast10k_out <- mclapply(gr10k_part_filtered,
                           FUN = function (x)beast(x$log2fc_weighted,
                                                   season = "none",
                                                   tcp.minmax = c(0,200),
                                                   quiet = TRUE),
                           mc.cores=10)

  # add trend$Y to the original data gr10k_part_filtered
  for (i in seq_along(gr10k_part_filtered)){
    gr10k_part_filtered[[i]]$score <- beast10k_out[[i]]$trend$Y
  }
  # merge the data back
  gr10k_bw <- unlist(gr10k_part_filtered)
  # normalize the score to median
  gr10k_bw$score <- gr10k_bw$score - median(gr10k$log2fc_weighted)
  # save the data
  message("Saving 10k window results")
  out1 <- paste0(prefix, "_10k.bw")
  export(gr10k_bw, out1, format="bigWig")
}


# calculate beast on each chromosome or scaffold on normalized data, min size is 10 intervals
message("Calculating BEAST on 10k windows on normalized data")
beast10k_out_norm <- mclapply(gr10k_part_filtered,
                              FUN = function (x)beast(x$log2fc_weighted_norm,
                                                      season = "none",
                                                      tcp.minmax = c(0,200),
                                                      quiet = TRUE),
                              mc.cores=10)
# add trend$Y to the original data gr10k_part_filtered
for (i in seq_along(gr10k_part_filtered)){
  gr10k_part_filtered[[i]]$score <- beast10k_out_norm[[i]]$trend$Y
}

# merge the data back
gr10k_bw_norm <- unlist(gr10k_part_filtered)
# normalize the score to median
gr10k_bw_norm$score <- gr10k_bw_norm$score - median(gr10k$log2fc_weighted_norm)
# save the data
message("Saving 10k window results on normalized data")
out1_norm <- paste0(prefix, "_10k_norm.bw")
export(gr10k_bw_norm, out1_norm, format="bigWig")


# save coverage tracks for further analysis
chip_coverage <- gr10k_bw_norm
chip_coverage$score <- gr10k_bw_norm$chip
export(chip_coverage, paste0(prefix, "_chip_coverage_10k.bw"), format="bigWig")
input_coverage <- gr10k_bw_norm
input_coverage$score <- gr10k_bw_norm$input
export(input_coverage, paste0(prefix, "_input_coverage_10k.bw"), format="bigWig")
# save normalized data
chip_coverage_norm <- gr10k_bw_norm
chip_coverage_norm$score <- gr10k_bw_norm$chip_norm
export(chip_coverage_norm, paste0(prefix, "_chip_coverage_norm_10k.bw"), format="bigWig")
input_coverage_norm <- gr10k_bw_norm
input_coverage_norm$score <- gr10k_bw_norm$input_norm
export(input_coverage_norm, paste0(prefix, "_input_coverage_norm_10k.bw"), format="bigWig")





# calculate beast on each chromosome or scaffold, min size is 100 intervals
gr_part_filtered <- gr_part[sapply(gr_part, function(x)length(x$log2fc_weighted) > 100)]
if (!opt$normalized_only){
  message("Calculating BEAST on 2k windows")
  beast_out <- mclapply(gr_part_filtered,
                        FUN = function (x)beast(x$log2fc_weighted,
                                                season = "none",
                                                tcp.minmax = c(0,200),
                                                quiet = TRUE
                        ), mc.cores=10)
  # add trend$Y to the original data gr_part_filtered
  for (i in seq_along(gr_part_filtered)){
    gr_part_filtered[[i]]$score <- beast_out[[i]]$trend$Y
  }
  gr_bw <- unlist(gr_part_filtered)
  # normalize the score to median
  gr_bw$score <- gr_bw$score - median(gr$log2fc_weighted)
  out2 <- paste0(prefix, "_2k.bw")
  export(gr_bw, out2, format="bigWig")
}


# calculate beast on each chromosome or scaffold on normalized data, min size is 10 intervals
message("Calculating BEAST on 2k windows on normalized data")
beast_out_norm <- mclapply(gr_part_filtered,
                           FUN = function (x)beast(x$log2fc_weighted_norm,
                                                   season = "none",
                                                   tcp.minmax = c(0,200),
                                                   quiet = TRUE
                           ), mc.cores=10)
# add trend$Y to the original data gr_part_filtered
message("Saving 2k window results on normalized data")

for (i in seq_along(gr_part_filtered)){
  gr_part_filtered[[i]]$score <- beast_out_norm[[i]]$trend$Y
}
gr_bw_norm <- unlist(gr_part_filtered)
# normalize the score to median
gr_bw_norm$score <- gr_bw_norm$score - median(gr$log2fc_weighted_norm)
out2_norm <- paste0(prefix, "_2k_norm.bw")
export(gr_bw_norm, out2_norm, format="bigWig")

# save coverage tracks for further analysis
chip_coverage <- gr_bw_norm
chip_coverage$score <- gr_bw_norm$chip
export(chip_coverage, paste0(prefix, "_chip_coverage_2k.bw"), format="bigWig")
input_coverage <- gr_bw_norm
input_coverage$score <- gr_bw_norm$input
export(input_coverage, paste0(prefix, "_input_coverage_2k.bw"), format="bigWig")
# save normalized data
chip_coverage_norm <- gr_bw_norm
chip_coverage_norm$score <- gr_bw_norm$chip_norm
export(chip_coverage_norm, paste0(prefix, "_chip_coverage_norm_2k.bw"), format="bigWig")
input_coverage_norm <- gr_bw_norm
input_coverage_norm$score <- gr_bw_norm$input_norm
export(input_coverage_norm, paste0(prefix, "_input_coverage_norm_2k.bw"), format="bigWig")




save.image(paste0(prefix, ".RData"))
message("Done")

