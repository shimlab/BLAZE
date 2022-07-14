# BLAZE (Barcode identification from Long reads for AnalyZing single cell gene Expression)
[![Github All Releases](https://img.shields.io/github/downloads/shimlab/BLAZE/total.svg)]()

**Important Notes:** This repo is actively being updated. Please make sure you have the latest release.

## Keywords:
Oxford Nanopore sequencing, Demultiplexing, Single Cell, Barcode.

# Overview
Combining single-cell RNA sequencing with Nanopore long-read sequencing enables isoform level analysis in single cells. However, due to the relatively high error rate in Nanopore reads, the demultiplexing of cell barcodes(cellBCs) and Unique molecular Identifiers (UMIs) is challenging. This tool enables the accurate identification of cellBCs solely from Nanopore reads. The output of BLAZE is a cellBC whitelist that can be utilised by downstream tools such as FLAMES to quantify genes and isoforms in single cells. 

# Installation

Download and unzip the [latest release](https://github.com/shimlab/BLAZE/releases/download/v1.0.0/BLAZE_v1.0.0.zip), the scripts are in `bin`. Or with command line:
```
wget https://github.com/shimlab/BLAZE/releases/download/v1.0.0/BLAZE_v1.0.0.zip
unzip BLAZE_v1.0.0.zip
cd BLAZE
```

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
This script has been tested on Chromium **Single Cell 3ʹ gene expression v3** and should be able to work on **Chromium Single Cell 3ʹ gene expression v2**, but it doesn't support any 10X 5' gene expression kits.

**Input:** 
 * *Folder of the fastq files*
 * *10X barcode whitelist*: a file containing all possible 10x barcodes ([more details](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)). The 10X barcode whitelists for both the 10X Single Cell 3ʹ gene expression v2 and v3 chemistries have been packed in `10X_bc/`. By default, this module assumes v3 chemistry and uses the file `10X_bc/3M-february-2018.zip`. You may specify `--kit-version=v2` to use the whitelist for v2 chemistry, or provide your own whitelist by specifying `--full-white-list=<filename>`. Note that barcodes outside this whitelist will never be found in the output.

 * *expected number of cells*: In the current version, the expected number of cells is a required input (specify `--expect-cells=xx`). Note that the output is NOT sensitive to the specified number, but a rough number is needed to determine the count threshold to output the whitelist. If users need to adjust the count threshold after running this script, They don't have to re-run this module but can use `bin/update_whitelist.py` (run `python3 bin/update_whitelist.py -h` for more details).

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
