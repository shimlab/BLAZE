# BLAZE (Barcode identification from Long reads for AnalyZing single cell gene Expression)
**Important Notes:** This repo is actively being updated. Please make sure you have the latest commit. To do this, you can run `git pull` before you run the script, which may potential avoid some errors.


## Keywords:
Oxford Nanopore sequencing, Demultiplexing, Single Cell, Barcode.

# Overview
Combining the single-cell RNA sequencing technology with Nanopore long sequencing enables the isoform level analysis in single cell. However, due to the relatively high error rate in Nanopore reads, the demultiplexing of cell barcode(cellBC) and Unique molecular Identifier (UMI) is challenging. This is a tool for identify cellBC solely from Nanopore reads.

# Installation

```
git clone https://github.com/youyupei/BLAZE.git
cd BLAZE
```
The scripts are in `bin`.

## Conda environment
All [dependencies](#dependencies) can be installed using "conda", the configuration file can be found at `conda_env/`:
```
conda env create -f conda_env/environment.yml
conda activate blaze
```

## <a name="dependencies"></a>Dependencies
* `Biopython`
* `pandas`
* `numpy`
* `tqdm`
* `matplotlib`



# Module:

## `get_raw_bc`: Get raw BC and BC whitelist from fastq files
This script has been tested on Chromium **Single Cell 3ʹ gene expression v3** and should be able to work on **Chromium Single Cell 3ʹ gene expression v2**, but it doesn't support any 10X 5' gene expression kit.

**Input:** 
 * *Folder of the fastq files*
 * *10X BC barcode whitelist*: list containing all the possible barcode ([more details](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)). By default, this module assumes **Single Cell 3ʹ gene expression v3** and use the file `10X_bc/3M-february-2018.zip`. Please specify a different file if you are using **Single Cell 3ʹ gene expression v2** kit.

 * *expected number of cells*: In the current version, the expected number (roughly) of cells is required input (specify `--expect-cells=xx`). 

**Example code:**
```
python3 get_raw_bc.py --expect-cells=1000 --threads=12 path/to/fastq_pass
```

**More information:**
```
python3 get_raw_bc.py -h
```

**Output:**
1. Print when running: stats of the raw barcode in reads
2. Raw barcode in each read, default filename: raw_bc.csv. It contains 3 columns
    * col1: read id
    * col2: raw barcode in read
    * col3: minimum Phred score of the bases in barcode
   
    **Note:** col2 and 3 will be empty if barcode not found. 
3. Cell-ranger style barcode whitelist, default filename: whitelist.csv
4. Knee plot using the raw barcode with high quality.

**Note:**
1. The raw barcodes are the 16nt sequence after identifed 10X adaptor within each read without correction for any basecalling error.
2. This module process individual FASTQ files in the input folder in separate CPUs to achieve multiprocessing. This means that the multiprocessing will NOT work if the input folder contains only one large FASTQ file. Splitting is recommended in this case.
