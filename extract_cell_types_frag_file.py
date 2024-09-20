## Process fragment file and extract fragments for all cell types based on cell barcode annotations,
## save fragments for each cell type into separate files and create tabix indices

import argparse
import os
import gzip
from pathlib import Path
import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description = ("Split fragments file into multiple files based"
                                  "on annotated cell types for each cell barcode."))
parser.add_argument("-f", "--fragment_file", type = str, required = True,
                    help = "Input fragment file to process.")
parser.add_argument("-c", "--cell_annot_file", type = str, required = True,
                    help = "File containing cell type annotations for each fragment.")
parser.add_argument("-l", "--lane", type = str, required = False, default = None,
                    help = "Sequencing lane id used to assign barcodes to lanes.")
parser.add_argument("-o", "--output_directory", type = str, required = False, default = ".",
                    help = "Output directory in which fragment files will be written.")
parser.add_argument('--compress_and_index', action = 'store_true',
                    help = "Trigger sorting, compressing and tabix indexing of output files.")
args = parser.parse_args()

# create lane suffix
if args.lane is None:
  lane_suffix = ""
else:
  lane_suffix = "-" + str(args.lane)

# load cell annotation file
cell_annot = pd.read_csv(args.cell_annot_file)
cell_annot = cell_annot.rename(columns={'Unnamed: 0': 'barcode'})

# create dictionary with output file for every cell barcode based on annotated cell type
basename = Path(args.fragment_file).with_suffix('').stem
cell_annot['outfile'] = args.output_directory + "/" + basename + "." + cell_annot['sample'] + ".tsv"
output_files_dict = cell_annot.set_index('barcode')['outfile'].to_dict()

# create output directory if needed
os.makedirs(args.output_directory, exist_ok = True)

# open connections to all output files and create dictionary with file connections per output file
output_files = list(set(output_files_dict.values()))
file_connections = {}
for file_path in output_files:
    file = open(file_path, 'w')
    file_connections[file_path] = file

# open fragment file and write each fragment to cell type specific output file
with gzip.open(args.fragment_file, "rt") as infile:
  for fragment_number, fragment in enumerate(infile, start = 1):
    
    # report how many fragments have been processed
    if fragment_number % 1000000 == 0:
      frag_processed = int(fragment_number / 1000000)
      print(f"{frag_processed} million fragments processed")
    
    # split fragment line into fields and add lane to barcode
    fragment = fragment.rstrip().split('\t')
    fragment[3] = fragment[3] + lane_suffix
    
    # get output file for given cell barcode
    barcode_outfile = output_files_dict.get(fragment[3])
  
    # if no output file is found, skip to next fragment, else write fragment to output file
    if barcode_outfile is None:
      continue
    else:
      fragment = '\t'.join(fragment)
      file_connections[barcode_outfile].write(fragment + '\n')
      
# close all output file connections   
for file in list(file_connections.values()):
  file.close()
    
# compress all output files using bgzip and create indives using tabix if specified
if args.compress_and_index is True:
  print("Compressing output files and creating indices")
  for file in output_files:
    os.system(f'sortBed -i {file} | bgzip > {file}.gzip')
    os.system(f'tabix -p bed {file}.gzip')
    os.system(f'rm {file}')

print("Done")
