#!/usr/bin/env Rscript
library(optparse)

if (FALSE){
  # for testing
  # make example data
  x <- c(rnorm(1000,2,0.5),
         rnorm(1000,0,0.5),
         rnorm(2000,2,0.5),
         rnorm(3000,0,0.5),
         rnorm(4000,2,0.5),
         rnorm(5000,0,0.5),
         rnorm(6000,2,0.5)
  )

}
calculate_intervals <- function(length_of_x, number_of_pieces, min_interval_size = NULL) {
  # Define the overlap fraction (66% overlap)
  overlap_fraction <- 0.66
  # Calculate the length of each interval (L)
  L_denominator <- overlap_fraction + (1 - overlap_fraction) * number_of_pieces
  interval_length <- length_of_x / L_denominator
  # Calculate the step size between intervals
  step_size <- interval_length * (1 - overlap_fraction)
  # Generate the starting indices
  starts <- 1 + step_size * (0:(number_of_pieces - 1))
  # Generate the ending indices
  ends <- starts + interval_length - 1
  # Ensure the ending indices do not exceed the length of the data vector
  ends <- pmin(ends, length_of_x)
  # Round the indices to the nearest integers
  starts <- floor(starts)
  ends <- ceiling(ends)
  # Adjust the intervals if the last interval is too small
  if (!is.null(min_interval_size)) {
    # Calculate the size of the last interval
    last_interval_size <- ends[length(ends)] - starts[length(starts)] + 1
    if (last_interval_size < min_interval_size) {
      # Merge the last interval with the previous one
      if (length(starts) > 1) {
        # Remove the last start and end
        starts <- starts[-length(starts)]
        ends <- ends[-length(ends)]
        # Extend the previous end to the last end
        ends[length(ends)] <- length_of_x
      } else {
        # If there is only one interval, set its end to length_of_x
        ends[1] <- length_of_x
      }
    }
  }
  return(list(starts = starts, ends = ends))
}

calculate_mean_between_change_points <- function(x, change_points){
  # calculate the mean between change points
  n <- length(x)
  # add firts and last change points - 0 and n
  change_points <- unique(sort(c(0, change_points, n)))
  mean_values <- numeric(length(change_points) - 1)
  width_of_intervals <- diff(change_points)
  for (i in seq_along(mean_values)){
    mean_values[i] <- mean(x[(change_points[i] + 1):change_points[i + 1]])
  }
  mean_values_x <- rep(mean_values, width_of_intervals)
  return(mean_values_x)
}

add_means_to_beast <- function(beast_out){

  ncp <- min(which.max(beast_out$trend$ncpPr) -1, length(beast_out$trend$cp))
  if (ncp == 0){
    message("No change points found")
    beast_out$trend$Y2 <- rep(mean(beast_out$data), length(beast_out$data))
    return(beast_out)
  }

  print(beast_out$trend$ncpPr)
  change_points <- beast_out$trend$cp[1:ncp]
  beast_out$trend$Y2 <- calculate_mean_between_change_points(beast_out$data, change_points)
  return(beast_out)
}


beast_wrapper <- function (x, ...){
  # if x is longer, split it into smaller parts and run beast on each part
  # there must be overlap between the parts to avoid edge effects - 30% of the window size
  args_list <- list(...)
  max_length_of_x <- 300000
  number_of_parts <- round(length(x) / max_length_of_x)
  if (number_of_parts <= 1){
    message('Running beast on the whole data')
    message('Length of x:', length(x))
    beast_args <- c(list(y = x), args_list)
    out_tmp <- do.call(beast, beast_args)
    # out_tmp <- beast(x,...)
    out_tmp <- add_means_to_beast(out_tmp)
    out_score <- out_tmp$trend$Y # 2
    out_sd <- out_tmp$trend$SD


    if (any(is.na(out_score))){
      message("NA in the result, trying to run beast again")
      message("#################################################")
      out_tmp <- do.call(beast, beast_args)
      # out_tmp <- beast(x,...)
      out_tmp <- add_means_to_beast(out_tmp)
      out_score <- out_tmp$trend$Y # 2
      out_sd <- out_tmp$trend$SD

      if (any(is.na(out_score))){
        message("NA again in the result")
        print(x)
        stop("NA in the BEAST results")
      }else{
        return(list(score = out_score, sd = out_sd))
      }
    }else{
      return(list(score = out_score, sd = out_sd))
    }
  }
  message('Splitting the data into ', number_of_parts, ' parts')
  message('Length of x:', length(x))
  message("--------------------------------")
  intervals <- calculate_intervals(length(x), number_of_parts, min_interval_size = 1000)
  starts <- intervals$starts
  ends <- intervals$ends


  beast_out <- list()
  for (i in seq_along(starts)){
    x_part <- x[starts[i]:ends[i]]
    beast_args <- c(list(y = x_part), args_list)
    beast_out[[i]] <- do.call(beast, beast_args)
    # beast_out[[i]] <- beast(x_part,...)
    if (any(is.na(beast_out[[i]]$trend$Y))){
      message("NA in the result, trying to run beast again")
      message("#################################################")
      beast_out[[i]] <- do.call(beast, beast_args)
      # beast_out[[i]] <- beast(x_part,...)
      if (any(is.na(beast_out[[i]]$trend$Y))){
        message("NA again in the result")
        print(x_part)
        stop("NA in the result")
      }
    }
    beast_out[[i]] <- add_means_to_beast(beast_out[[i]])
  }
  # concatenate the results,   in overlapping parts calculate the mean
  n_results <- rep(0,length(x))
  beast_merged_score <- rep(0,length(x))
  beast_merged_sd <- rep(0,length(x))
  for (i in seq_along(beast_out)){
    n_results[starts[i]:ends[i]] <- n_results[starts[i]:ends[i]] + 1
    beast_merged_score[starts[i]:ends[i]] <- beast_merged_score[starts[i]:ends[i]] + beast_out[[i]]$trend$Y #2
    beast_merged_sd[starts[i]:ends[i]] <- beast_merged_sd[starts[i]:ends[i]] + beast_out[[i]]$trend$SD
  }
  beast_merged_score <- beast_merged_score / n_results
  beast_merged_sd <- beast_merged_sd / n_results

  return(list(score = beast_merged_score, sd = beast_merged_sd))

}

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
              help = "only calculate normalized data"),
  make_option(c("-S", "--split"), action='store_true', default = FALSE,
              help="split the input into smaller parts and run beast on each part, use thi option if you analyze large genome" )
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
count_list <- calcNormFactors(count_list, method="upperquartile")
normalized_counts <- cpm(count_list)

counts10k_list <- DGEList(assay(counts10k))
counts10k_list <- calcNormFactors(counts10k_list, method="upperquartile")
normalized_counts10k <- cpm(counts10k_list)

# ratio of number of reads in chip and input

gr <- rowRanges(counts)
gr$chip <- assay(counts)[,1]
gr$input <- assay(counts)[,2]
norm_factor <- sum(gr$chip) / sum(gr$input)
gr$log2fc <- log2(gr$chip + 1) - log2(gr$input + 1)
gr$log2fc_weighted <- gr$log2fc * log(gr$chip + gr$input + 1)
gr$chip_norm <- normalized_counts[,1]
gr$input_norm <- normalized_counts[,2]

gr$chip_norm <- gr$chip/norm_factor
gr$input_norm <- gr$input

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

gr10k$chip_norm <- gr10k$chip/norm_factor
gr10k$input_norm <- gr10k$input



min_value <- min(normalized_counts10k[normalized_counts10k != 0])
gr10k$log2fc_norm <- log2(gr10k$chip_norm + min_value) - log2(gr10k$input_norm + min_value)
gr10k$log2fc_weighted_norm <- gr10k$log2fc_norm * log(gr10k$chip_norm + gr10k$input_norm + min_value)
gr10k_part <- split(gr10k, seqnames(gr10k))





# calculate beast on each chromosome or scaffold, min size is 10 intervals
gr10k_part_filtered <- gr10k_part[sapply(gr10k_part, function(x)length(x$log2fc_weighted) > 100)]
if (!opt$normalized_only){
  message("Calculating BEAST on 10k windows")
  if (opt$split){
    save.image(paste0(prefix, ".RData"))
    beast10k_out <- mclapply(gr10k_part_filtered,
                             FUN = function (x)beast_wrapper(x$log2fc_weighted,
                                                             season='none',
                                                             tcp.minmax=c(1,300),
                                                             torder.minmax=c(0,1),
                                                             tseg.min=5,
                                                             quiet = TRUE),
                             mc.cores=10)
    for (i in seq_along(gr10k_part_filtered)){
      gr10k_part_filtered[[i]]$score <- beast10k_out[[i]]$score
      gr10k_part_filtered[[i]]$sd <- beast10k_out[[i]]$sd
    }
  } else {
    beast10k_out <- mclapply(gr10k_part_filtered,
                             FUN = function (x)beast(x$log2fc_weighted,
                                                     season='none',
                                                     tcp.minmax=c(1,300),
                                                     torder.minmax=c(0,1),
                                                     tseg.min=5,
                                                     quiet = TRUE,
                             ),
                             mc.cores=10)

    # add trend$Y to the original data gr10k_part_filtered
    for (i in seq_along(gr10k_part_filtered)){
      gr10k_part_filtered[[i]]$score <- beast10k_out[[i]]$trend$Y
      gr10k_part_filtered[[i]]$sd <- beast10k_out[[i]]$trend$SD
    }
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
if (opt$split){

  save.image(file = paste0(prefix, "x.RData"))
  beast10k_out_norm <- mclapply(gr10k_part_filtered,
                                FUN = function (x)beast_wrapper(x$log2fc_weighted_norm,
                                                                 season = "none",
                                                                tcp.minmax=c(1,300),
                                                                torder.minmax=c(0,1),
                                                                tseg.min=5,
                                                                quiet = TRUE),
                                mc.cores=10)
  for (i in seq_along(gr10k_part_filtered)){
    gr10k_part_filtered[[i]]$score <- beast10k_out_norm[[i]]$score
    gr10k_part_filtered[[i]]$sd <- beast10k_out_norm[[i]]$sd
  }

} else {
  beast10k_out_norm <- mclapply(gr10k_part_filtered,
                                FUN = function (x)beast(x$log2fc_weighted_norm,
                                                        season = "none",
                                                        tcp.minmax=c(1,300),
                                                        torder.minmax=c(0,1),
                                                        tseg.min=5,
                                                        quiet = TRUE),
                                mc.cores=10)
# add trend$Y to the original data gr10k_part_filtered
  for (i in seq_along(gr10k_part_filtered)){
    gr10k_part_filtered[[i]]$score <- beast10k_out_norm[[i]]$trend$Y
    gr10k_part_filtered[[i]]$sd <- beast10k_out_norm[[i]]$trend$SD
  }
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
  if (opt$split){
    beast_out <- mclapply(gr_part_filtered,
                          FUN = function (x)beast_wrapper(x$log2fc_weighted,
                                                          season = "none",
                                                          tcp.minmax=c(1,300),
                                                          torder.minmax=c(0,1),
                                                          tseg.min=5,
                                                          quiet = TRUE
                          ), mc.cores=10)
    for (i in seq_along(gr_part_filtered)){
      gr_part_filtered[[i]]$score <- beast_out[[i]]$score
      gr_part_filtered[[i]]$sd <- beast_out[[i]]$sd
    }
  } else {
    beast_out <- mclapply(gr_part_filtered,
                          FUN = function (x)beast(x$log2fc_weighted,
                                                  season = "none",
                                                  tcp.minmax=c(1,300),
                                                  torder.minmax=c(0,1),
                                                  tseg.min=5,
                                                  quiet = TRUE
                          ), mc.cores=10)
    # add trend$Y to the original data gr_part_filtered
    for (i in seq_along(gr_part_filtered)){
      gr_part_filtered[[i]]$score <- beast_out[[i]]$trend$Y
      gr_part_filtered[[i]]$sd <- beast_out[[i]]$trend$SD
    }
  }
  gr_bw <- unlist(gr_part_filtered)
  # normalize the score to median
  gr_bw$score <- gr_bw$score - median(gr$log2fc_weighted)
  out2 <- paste0(prefix, "_2k.bw")
  export(gr_bw, out2, format="bigWig")
}


# calculate beast on each chromosome or scaffold on normalized data, min size is 10 intervals
message("Calculating BEAST on 2k windows on normalized data")
if (opt$split){
  beast_out_norm <- mclapply(gr_part_filtered,
                             FUN = function (x)beast_wrapper(x$log2fc_weighted_norm,
                                                             season = "none",
                                                             tcp.minmax=c(1,300),
                                                             torder.minmax=c(0,1),
                                                             tseg.min=5,
                                                             quiet = TRUE
                             ), mc.cores=10)
  for (i in seq_along(gr_part_filtered)){
    gr_part_filtered[[i]]$score <- beast_out_norm[[i]]$score
    gr_part_filtered[[i]]$sd <- beast_out_norm[[i]]$sd
  }
} else {
  beast_out_norm <- mclapply(gr_part_filtered,
                             FUN = function (x)beast(x$log2fc_weighted_norm,
                                                     season = "none",
                                                     tcp.minmax=c(1,300),
                                                     torder.minmax=c(0,1),
                                                     tseg.min=5,
                                                     quiet = TRUE
                             ), mc.cores=10)
  # add trend$Y to the original data gr_part_filtered
  message("Saving 2k window results on normalized data")
  for (i in seq_along(gr_part_filtered)){
    gr_part_filtered[[i]]$score <- beast_out_norm[[i]]$trend$Y
    gr_part_filtered[[i]]$sd <- beast_out_norm[[i]]$trend$SD
  }
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

