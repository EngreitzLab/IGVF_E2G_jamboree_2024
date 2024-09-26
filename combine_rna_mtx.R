## Combine multipe RNA matrix .mtx files into one. Used to combine data from multiple lanes

library(optparse)

## Parse command line arguments --------------------------------------------------------------------

# create arguments list
option_list = list(
  make_option(c("-i", "--input_files"), type = "character", default = NULL,
              help = "Comma separated list of input RNA matrix files", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "Output file for combined RNA matrix", metavar = "character")
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
required_args <- c("input_files", "output_file")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

## Define functions --------------------------------------------------------------------------------

# load one matrix file, assuming cell and gene name files are in the same direcotry
load_mtx <- function(mtx_file) {
  
  # get barcodes and genes files based on mtx_file name
  cell_type <- sub("cells_x_genes\\..+\\.(.+)\\.mtx", "\\1", basename(mtx_file))
  barcodes <- file.path(dirname(mtx_file), paste0("cells_x_genes.barcodes.", cell_type, ".txt"))
  genes <- file.path(dirname(mtx_file), paste0("cells_x_genes.genes.names.", cell_type, ".txt"))
  
  # load matrix and set column and row names
  mtx <- readMM(mtx_file)
  rownames(mtx) <- readLines(barcodes)
  colnames(mtx) <- readLines(genes)
  
  return(mtx)
  
}

## Combine RNA matrices ----------------------------------------------------------------------------

# required packages
suppressPackageStartupMessages({
  library(Matrix)
})

# load all input matrix files
message("Loading .mtx files...")
input_mtx <- unlist(strsplit(opt$input_files, split = ","))
input_mtx <- trimws(input_mtx)  # remove any leading and trailing white spaces
mtx <- lapply(input_mtx, FUN = load_mtx)

# check that column names are the same for all matrices
mtx_colnames <- lapply(mtx, FUN = colnames)
mtx_colnames_identical <- all(sapply(mtx_colnames, FUN = identical, mtx_colnames[[1]]))
if (mtx_colnames_identical == FALSE) {
  stop("Input matrix column names aren't identical.")
}

# combine into one matrix by appending rows
message("Combining .mtx files...")
mtx <- do.call(rbind, mtx)

# create output directory if needed
outdir <- dirname(opt$output_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# output files for barcode and gene names
cell_type <- sub("cells_x_genes\\..+\\.(.+)\\.mtx", "\\1", basename(opt$output_file))
barcodes_outfile <- file.path(outdir, paste0("cells_x_genes.barcodes.", cell_type, ".txt"))
genes_outfile <- file.path(outdir, paste0("cells_x_genes.genes.names.", cell_type, ".txt"))

# write matrix, column and row names to output files
message("Writing output files...")
invisible(writeMM(mtx, file = opt$output_file))
writeLines(rownames(mtx), con = barcodes_outfile)
writeLines(colnames(mtx), con = genes_outfile)

message("Done!")
