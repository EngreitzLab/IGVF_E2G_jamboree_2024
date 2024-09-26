## Extract data on all cell types from an IGVF multiome RNA matrix and write to separate matrix
## files, one per cell type

library(optparse)

## Parse command line arguments --------------------------------------------------------------------

# create arguments list
option_list = list(
  make_option(c("-m", "--matrix_file"), type = "character", default = NULL,
              help = "Path to input RNA matrix file to process.", metavar = "character"),
  make_option(c("-c", "--cell_annot_file"), type = "character", default = NULL,
              help = "File containing cell type annotations for each fragment.",
              metavar = "character"),
  make_option(c("-l", "--lane"), type = "character", default = NULL,
              help = "Sequencing lane id used to assign barcodes to lanes.", metavar = "character"),
  make_option(c("-s", "--sample_col"), type = "character", default = "sample",
              help = paste("Column name in cell type annotations table containing cell type or",
                           "sample information (default: 'sample')"),
              metavar = "character"), 
  make_option(c("-o", "--output_directory"), type = "character", default = ".",
              help = "Output directory in which matrix files will be written.",
              metavar = "character")
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
required_args <- c("matrix_file", "cell_annot_file")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

## Define functions --------------------------------------------------------------------------------

# function to extract barcodes for one cell barcode
extract_cell_type <- function(cell_type, mtx, cell_annot, sample_col, outfiles) {
  
  message("Extracting RNA data for cell type ", cell_type)
  
  # extract data on all cell barcodes for given cell type
  barcodes_cell_type <- pull(filter(cell_annot, !!sym(sample_col) == cell_type), 1)
  mtx_cell_type <- mtx[intersect(barcodes_cell_type, rownames(mtx)), ]
  
  # write matrix to output file
  writeMM(mtx_cell_type, file = outfiles[[cell_type]]$mtx)
  
  # write barcodes and genes to files (column and row names for matrix)
  writeLines(rownames(mtx_cell_type), con = outfiles[[cell_type]]$cells)
  writeLines(colnames(mtx_cell_type), con = outfiles[[cell_type]]$genes)
  
}

## Process matrix file -----------------------------------------------------------------------------

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
})

# create lane suffix
if (!is.null(opt$lane)) {
  lane_suffix <- paste0("-", opt$lane)
} else {
  lane_suffix <- NULL
}

# load RNA matrix and list of cell barcodes and genes
mtx <- readMM(opt$matrix_file)
cells <- readLines(file.path(dirname(opt$matrix_file), "cells_x_genes.barcodes.txt"))
genes <- readLines(file.path(dirname(opt$matrix_file), "cells_x_genes.genes.names.txt"))

# set column (genes) and row (cells) names
colnames(mtx) <- genes
rownames(mtx) <- paste0(cells, lane_suffix)

# load cell annotations
cell_annot <- fread(opt$cell_annot_file)

# create output directory if needed
outdir <- opt$output_directory
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# get all cell types in matrix
cell_types <- unique(cell_annot[[opt$sample_col]])

# create output files for all cell types
mtx_basename <- tools::file_path_sans_ext(basename(opt$matrix_file))
mtx_outfiles <- file.path(outdir, paste0(mtx_basename, ".", cell_types, ".mtx"))
cells_outfiles <- file.path(outdir, paste0("cells_x_genes.barcodes.", cell_types, ".txt"))
genes_outfiles <- file.path(outdir, paste0("cells_x_genes.genes.names.", cell_types, ".txt"))

# create list with all output files per cell type
outfiles <- tibble(cell_type = cell_types, mtx = mtx_outfiles, cells = cells_outfiles,
                   genes = genes_outfiles) %>% 
  split(x = .[, -1], f = .$cell_type)

# extract RNA matrices for each cell type
invisible(lapply(cell_types, FUN = extract_cell_type, mtx = mtx, cell_annot = cell_annot,
                 sample_col = opt$sample_col, outfiles = outfiles))
message("Done")
