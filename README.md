# SC_NanoDemulti (temp name)
Demultiplexing 10X Single Cell Nanopore reads without short reads.

## Keywords:
Oxford Nanopore sequencing, Demultiplexing, Single Cell, Barcode.

# Overview
Combining the single-cell RNA sequencing technology with Nanopore long sequencing enables the isoform level analysis in single cell. However, due to the relatively high error rate  in Nanopore reads, the demultiplexing of cell barcode(cellBC) and Unique molecular Identifier (UMI) is challenging. This is a tool for identify cellBC and UMI solely from Nanopore reads.

# Installation

```
git clone https://github.com/youyupei/SC_NanoDemulti.git
```

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

## `get_raw_bc`:
```
Summary: 
    Get raw BC and BC whitelist from fastq files

Usage: python3 {argv[0]} [OPTIONS] <fastq directory>

Options:
    -h, --help
        Print this help message.
    
    --expect-cells
        <INT>:  Expected number of cells. Default: not specified
    
    --minQ:
        <INT>: Minimum phred score for all bases in a raw BC. Reads whose 
        raw BC contains one or more bases with Q<minQ is not counted 
        in the "Raw BC rank plot". Default: --minQ=15
    
    --full-bc-whitelist=
        <path to file>: .txt file containing all the possible BCs. Users may provide
        their own whitelist. Default: 3M BC from 10X.
    
    --out-raw-bc
        <filename_prefix>: Output a csv file for the raw BC in each read. 
                            Default: --out-raw-bc=raw_bc
    
    --out-bc-whitelist
        <filename_prefix>: Output the whitelist identified from all the reads. 
                            Default: --out-bc-whitelist=whitelist
    
    --threads
        <INT>: Number of threads used <default: # of available cpus - 1> 
```
