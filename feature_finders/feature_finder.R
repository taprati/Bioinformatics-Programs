#!/usr/bin/env Rscript
# 
# Feature Finder
#

# Read in args
args = commandArgs(trailingOnly = TRUE)

# Assign args
input = args[1]
chrom = args[2]
type = args[3]
start = as.integer(args[4])
end = as.integer(args[5])

# Check for args
if (length(args) != 5) {
  cat("Usage: feature_finder.R input chromosome feature_type start_index end_index\n")
  stop("Invalid input. Not enough arguments!", call.=FALSE)
}

# Load the library.
suppressPackageStartupMessages(library(dplyr))

# Load in data
data <- read.delim(input, header=F, comment.char="#")

# Filter data by parameters
genes_of_interest <- data %>%
  filter(V1 == chrom, 
         V3 == type,
         (V4 >= start & V4 <= end) | (V5 >= start & V5 <= end) | (V4 <= start & V5 >= end)
         )

# Write genes of interest into output file
write.table(genes_of_interest, 
            file = "interesting_genes_R.out", 
            sep = "\t", 
            row.names = FALSE,
            col.names = FALSE, 
            quote = FALSE)

