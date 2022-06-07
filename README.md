# BLAZE (Barcode identification from Long reads for AnalyZing single cell gene Expression)

## Keywords:
Oxford Nanopore sequencing, Demultiplexing, Single Cell, Barcode.

# Overview
Combining the single-cell RNA sequencing technology with Nanopore long sequencing enables the isoform level analysis in single cell. However, due to the relatively high error rate  in Nanopore reads, the demultiplexing of cell barcode(cellBC) and Unique molecular Identifier (UMI) is challenging. This is a tool for identify cellBC solely from Nanopore reads.

# Installation

```
git clone https://github.com/youyupei/SC_NanoDemulti.git
unzip SC_NanoDemulti/10X_bc/3M-february-2018.zip
```
The scripts are in `bin`.

## Conda environment
All [dependencies](#dependencies) can be installed using "conda", the configuration file can be found at `conda_env/`:
```
conda env create -f conda_env/environment.yml
conda activate sc_demulti
```

## <a name="dependencies"></a>Dependencies
* `Biopython`
* `pandas`
* `numpy`
* `tqdm`
* `matplotlib`



# Module:

## `get_raw_bc`: Get raw BC and BC whitelist from fastq files
This script has been tested on Chromium **Single Cell 3ʹ gene expression v3** and should be able to work on **Chromium Single Cell 3ʹ gene expression v2**, but it doesn't support any 10X 5' gene expression kit 

**Input:** 
 * *Folder of the fastq files*
 * *10X BC barcode whitelist*: list containing all the possible barcode ([more details](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)). By default, this module assumes **Single Cell 3ʹ gene expression v3** and use the file `10X_bc/3M-february-2018.txt` (you can get this file by unzipping `10X_bc/3M-february-2018.zip`). Please specify a different file if you are using **Single Cell 3ʹ gene expression v2** kit.

 * *expected number of cells*: In the current version, the expected number (roughly) of cells is required input (specify `--expect-cells=xx`). 

**Example code:**
```
python3 get_raw_bc.py --expect-cells=1000 --threads=12 path/to/fastq_pass
```

**More information:**
```
python3 get_raw_bc.py -h
```
```
Summary: 
    Get raw BC and BC whitelist from fastq files

Usage: python3 get_raw_bc.py [Options] <fastq directory>

Options:
    -h, --help
        Print this help message.
    
    --full-bc-whitelist (required in current version)
        <path to file>: .txt file containing all the possible BCs. Users may provide
        their own whitelist. Default: 3M BC from 10X.
    
    --expect-cells (required in current version)
        <INT>:  Expected number of cells. Default: not specified
    
    --minQ:
        <INT>: Minimum phred score for all bases in a raw BC. Reads whose 
        raw BC contains one or more bases with Q<minQ is not counted 
        in the "Raw BC rank plot". Default: --minQ=15
    

    
    --out-raw-bc
        <filename_prefix>: Output a csv file for the raw BC in each read. 
                            Default: --out-raw-bc=raw_bc
    
    --out-bc-whitelist
        <filename_prefix>: Output the whitelist identified from all the reads. 
                            Default: --out-bc-whitelist=whitelist
    
    --threads
        <INT>: Number of threads used <default: # of available cpus - 1> 
```

Output:
1. Print when running: stats of the raw barcode in reads
2. Raw barcode in each read, default filename: raw_bc.csv. It contains 3 columns
    col1: read id
    col2: raw barcode in read
    col3: minimum Phred score of the bases in barcode
    **Note:** col2 and 3 will be empty if barcode not found. 
3. Cell-ranger style barcode whitelist, default filename: whitelist.csv