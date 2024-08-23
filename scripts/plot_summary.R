#!/usr/bin/env Rscript



# use rather optparse for argument parsing
library(optparse)
option_list <- list(
  make_option(c("-d", "--dir"), type="character",help="Directory with files to plot"),
  make_option(c("-c", "--chrom_sizes"), type="character", help="Path to chrom.sizes file"),
  make_option(c("-o", "--output"), type="character",help="Output file path - png file with plots")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
# all arguments are required - validate them
if (any(sapply(opt, is.null))){
  stop("All arguments are required")
}

dir_path <- opt$dir
chrom_sizes_path <- opt$chrom_sizes
output_path <- opt$output

library(rtracklayer)

# check files for plotting:



plot_type <- c(
  "chip_coverage_bs2000.bw" = 'p',
  "input_coverage_bs2000.bw" = 'p',
  "chip_vs_input_all.bs2000.bw" = 'r',
  "chip_vs_input_unique.bs2000.bw"= 'r',
  "epic2_all.bs2000.bedgraph"= 'r',
  "epic2_unique.bs2000.bedgraph"= 'r',
  "macs3_all_peaks.bedgraph"= 'r',
  "macs3_unique_peaks.bedgraph"= 'r',
  "peakBeast/peakBeast_10k_norm.bw"= 'r',
  "peakBeast/peakBeast_2k_norm.bw"= 'r'
)

plot_colors <- c(
  "chip_coverage_bs2000.bw" = '#0000FF20',
  "input_coverage_bs2000.bw" = '#0000FF20',
  "chip_vs_input_all.bs2000.bw" = '#0f5c18',
  "chip_vs_input_unique.bs2000.bw"= '#0f5c18',
  "epic2_all.bs2000.bedgraph"= '#5b05a6',
  "epic2_unique.bs2000.bedgraph"= '#5b05a6',
  "macs3_all_peaks.bedgraph"= '#704d02',
  "macs3_unique_peaks.bedgraph"= '#704d02',
  "peakBeast/peakBeast_10k_norm.bw"= '#87052640',
  "peakBeast/peakBeast_2k_norm.bw"= '#87052640'
)

# check that plot_colors and plot_type have the same keys in the same order
stopifnot(identical(names(plot_colors), names(plot_type)))


selected_fnames <- names(plot_type)

log_scale <- c(
  "chip_coverage_bs2000.bw",
  "input_coverage_bs2000.bw"
)



chr_sizes <- read.table(chrom_sizes_path, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE, col.names = c("seqname", "size"))
# sort by size - largest first
chr_sizes <- chr_sizes[order(chr_sizes$size, decreasing = TRUE),]
# add cumulative start position
cum_start <- c(0, cumsum(as.numeric(chr_sizes$size[-nrow(chr_sizes)])))
names(cum_start) <- chr_sizes$seqname
cum_end <- cumsum(as.numeric(chr_sizes$size))
names(cum_end) <- chr_sizes$seqname


# check if all files are present
selected_fnames_exist <- file.exists(file.path(dir_path, selected_fnames))

# use only existing files
selected_fnames_ok <- selected_fnames[selected_fnames_exist]

N_plots <- length(selected_fnames_ok)
png(output_path, width = 3000, height = 200 * N_plots, pointsize = 20)
par(mfrow = c(N_plots + 1, 1), mar = c(0, 4, 0, 1))
for (i in 1:N_plots){
  print(i)
  bw <- import(file.path(dir_path, selected_fnames_ok[i]))
  # recalculate coordinates  - using cumulatve starts
  S <- start(bw) + cum_start[as.vector(seqnames(bw))]
  E <- end(bw) + cum_start[as.vector(seqnames(bw))]
  H <- score(bw)
  color <- plot_colors[selected_fnames_ok[i]]
  if (selected_fnames_ok[i] %in% log_scale){
    H <- log2(H + 1)
  }

  plot(0, type = "n", xlim = c(0, max(cum_end)), ylim = c(min(H), max(H) * 1.2),
       ylab = "", axes = FALSE, xlab = "", xaxs = "i"
  )
  # line with cumulative end positions
  abline(v = cum_end, col = "#99999950", lty = 1, lwd = 3)
  #box()
  mtext(selected_fnames_ok[i], side = 3, line = -1.5)
  if (plot_type[selected_fnames_ok[i]] == 'p'){
    points(S, H, pch = 18, cex = 0.5, col = color, type='p')
  }
  if (plot_type[selected_fnames_ok[i]] == 'r'){
    rect(S, 0, E, H, col = color, border = color)
  }
  axis(2, las = 2)

}
# add sequence names to x axis - in the middle of the chromosome, no ticks
axis(1, at = (cum_end + cum_start)/2, labels = names(cum_end),
     las = 2, cex.axis = 2,
     lwd = 0, line = -1)
dev.off()



