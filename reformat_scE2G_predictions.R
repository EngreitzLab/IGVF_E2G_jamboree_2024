## Reformat pillar project predictors

library(optparse)

# Process input arguments --------------------------------------------------------------------------

# create arguments list
option_list = list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Path to transcripts input file", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "Path to output file", metavar = "character"),
  make_option(c("-c", "--cell_type"), type = "character", default = NULL,
              help = "Cell type", metavar = "character"),
  make_option(c("-m", "--method"), type = "character", default = NULL,
              help = "E2G method that produced the predictions", metavar = "character"),
  make_option(c("-v", "--version"), type = "character", default = NULL,
              help = "E2G method version", metavar = "character"),
  make_option(c("-t", "--threshold"), type = "character", default = NULL,
              help = "Used score threshold if applicable", metavar = "character") 
  
)

# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
  }
}

# check that all required parameters are provided
required_args <- c("input_file", "output_file", "method", "version")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

# Process file -------------------------------------------------------------------------------------

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load input file
pred <- fread(opt$input_file)

# get all score columns (all columns except EG-pair defining columns)
message("Reformatting predictions...")

# set cell type
if (!is.null(opt$cell_type)) {
  cell_type <- opt$cell_type
} else {
  cell_type <- unique(pred$CellType)
}

# create header lines
header <- c(
  paste("# Source:", opt$method),
  paste("# Version:", opt$version),
  "# GenomeBuild: GRCh38",
  "# URL: [add url]",
  "# Assays: 10x multiome",
  "# BiosampleAgnostic: False",
  paste("# BiosampleString:", cell_type)
)

# add threshold if applicable
if (!is.null(opt$threshold)) {
  header <- c(header, paste("# Threshold:", opt$threshold))
}

# add additional columns and extract output columns
pred <- pred %>% 
  mutate(ElementChr = chr, name = paste0(ElementChr, ":", start, "-", end), ElementStrand = ".") %>% 
  select(ElementChr, ElementStart = start, ElementEnd = end, ElementName = name,
         ElementClass = class, ElementStrand, GeneSymbol = TargetGene,
         GeneEnsemblID = TargetGeneEnsembl_ID, GeneTSS = TargetGeneTSS, Score = E2G.Score.qnorm)

# save to output file
message("Writing to output file...")
if (tools::file_ext(opt$output_file) == "gz") {
  
  # save to gzip compressed file
  tmp_file <- tools::file_path_sans_ext(opt$output_file)
  writeLines(header, con = tmp_file)
  fwrite(pred, file = tmp_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  system2("gzip", args = c("-f", tmp_file))
  
} else {
  
  # save to uncompressed file
  writeLines(header, con = opt$output_file)
  fwrite(pred, file = opt$output_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  
}

message("Done!")
