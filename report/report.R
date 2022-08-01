library(DuoClustering2018)
library(optparse)
library(SingleCellExperiment)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-r", "--result"), type="character", default=NULL,
              help="result RData", metavar="character"),
  make_option(c("-w", "--workspace"), type="character", default="./",
              help="work directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$result)){
  print_help(opt_parser)
  stop("Result RData must be supplied (-i).", call.=FALSE)
}

setwd(opt$workspace)

load(opt$result)
ind=unique(df$method)

method_colors <- c(CIDR = "#332288", FlowSOM = "#6699CC", PCAHC = "#88CCEE", 
                   PCAKmeans = "#44AA99", pcaReduce = "#117733",
                   RtsneKmeans = "#999933", Seurat = "#DDCC77", SC3svm = "#661100", 
                   SC3 = "#CC6677", TSCAN = "grey34", ascend = "orange", SAFE = "black",
                   monocle = "red", RaceID2 = "blue",my_seurat="yellow")

include_colors=method_colors[ind]

perf <- plot_performance(df, method_colors = include_colors)
names(perf)
png("median_ari_vs_k.png")
perf$median_ari_vs_k
dev.off()

#perf$median_ari_heatmap_truek

stab <- plot_stability(df, method_colors = include_colors)
names(stab)
png("stability_vs_k.png")
stab$stability_vs_k
dev.off()

entr <- plot_entropy(df, method_colors = include_colors)
names(entr)
png("entropy_vs_k.png")
entr$entropy_vs_k
dev.off()

timing <- plot_timing(df, method_colors = include_colors, 
                      scaleMethod = "RtsneKmeans")
names(timing)
png("time_boxplot.png")
timing$time_boxplot
dev.off()
