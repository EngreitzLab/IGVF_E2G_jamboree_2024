# Data processing code for the IGVF E2G pillar project jamboree 2024

### Install all dependencies
To install all python and R dependencies, create a conda enviroment using the provided conda
environment file:
```
conda env create -f e2g_jamboree_env.yml
```

To activate the enviroment after installation, simply use:
```
conda activate e2g_jamboree
```

### Extract fragments for each cell type from fragment files
To extract separate fragment files for all annotated cell types in a multiome fragment file, use the
`extract_cell_types_frag_file.py` python script. This takes a fragment file and a cell annotation
table as input. It also allows to specify a lane identifier to add to the cell barcodes in the
fragment file, which is important if an experiments has files from multiple lanes.

To extract fragments for all cell types, use a command similar to this:
```
python3 extract_cell_types_frag_file.py -f fragments.tsv.gz -c metadata_IGVF10.csv -l 1 -o output/fragments_lane1
```

By default, this script assumes that the cell type for each barcode is listed in the 'sample' column
in the cell annotation table. To specify another column, use the `-s` or `--sample_col` argument.

This will create a new directory called `output/fragments_lane1` containing fragment files for all
cell types that were annotated for the input fragment file. Use this script on fragment files for
all lanes of an experiment to get demultiplexed fragment files for all lanes.

To combined fragment files from multiple lanes for one cell type, we need to combine, sort, compress
and tabix index the files. This can be done using these commands:
```
cat output/fragments_lane1/fragments_cell_type.tsv output/fragments_lane2/fragments_cell_type.tsv > output/fragments_combined/fragments_cell_type.tsv
sort -k1,1 -k2,2n output/fragments_combined/fragments_cell_type.tsv > output/fragments_combined/fragments_cell_type.sorted.tsv
rm output/fragments_combined/fragments_cell_type.tsv  # delete file not longer needed
bgzip output/fragments_combined/fragments_cell_type.sorted.tsv
tabix -p bed output/fragments_combined/fragments_cell_type.sorted.tsv.gz

```

### Extract RNA counts for each cell type from RNA matrix
To extract RNA counts for all annotated cell types from a .mtx file, use the
`extract_cell_types_mtx.R` R script. This takes a RNA .mtx file and a cell annotation table as
input. It also allows to specify a lane identifier to add to the cell barcodes in the matrix file,
which is important if an experiments has files from multiple lanes.

To extract RNA counts for all cell types, use a command similar to this:
```
Rscript extract_cell_types_mtx.R -m cells_x_genes.total.mtx -c metadata_IGVF10.csv -l 1 -o output/rna_lane1
```

By default, this script assumes that the cell type for each barcode is listed in the 'sample' column
in the cell annotation table. To specify another column, use the `-s` or `--sample_col` argument.

The script will create a new directory called `output/rna_lane1` containing matrix, barcode and gene
name files for all cell types that were annotated for the input matrix. Use this script on matrix
files for all lanes of an experiment to get demultiplexed RNA matrix files for all lanes.

To combine RNA matrix files from multiple lanes for one cell type, simply use the
`combine_rna_mtx.R` R script. This takes a comma separated list of input RNA matrix files as input
and will combine them into the specified output file. This will also create combined barcode and
gene name files in the same directory.
```
Rscript combine_rna_mtx.R -i output/rna_lane1/cells_x_genes.total.cell_type.mtx,output/rna_lane2/cells_x_genes.total.cell_type.mtx -o output/rna_combined/cells_x_genes.total.cell_type.mtx
```
