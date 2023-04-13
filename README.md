<img src="logo.png" width="300"/>

# BLAZE (Barcode identification from Long reads for AnalyZing single cell gene Expression)
[![Github All Releases](https://img.shields.io/github/downloads/shimlab/BLAZE/total.svg)](https://github.com/shimlab/BLAZE/releases/download/v1.1.0/BLAZE_v1.1.0.zip)

**Important Notes:** This repo is actively being updated. Please make sure you have the latest release.

## Keywords:
Oxford Nanopore sequencing, Demultiplexing, Single Cell, Barcode.

# Overview
Combining single-cell RNA sequencing with Nanopore long-read sequencing enables isoform level analysis in single cells. However, due to the relatively high error rate in Nanopore reads, the demultiplexing of cell barcodes and Unique molecular Identifiers (UMIs) is challenging. This tool enables the accurate identification of barcodes solely from Nanopore reads. The output of BLAZE is a barcode whitelist that can be utilised by downstream tools such as FLAMES to quantify genes and isoforms in single cells. For a detailed description of how BLAZE works and its performance across different datasets, please see our [Genome Biology paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02907-y).

# Installation

Download and unzip the [latest release](https://github.com/shimlab/BLAZE/releases/download/v1.1.0/BLAZE_v1.1.0.zip), the scripts are in `bin`. Or with the command line:
```
wget https://github.com/shimlab/BLAZE/releases/download/v1.1.0/BLAZE_v1.1.0.zip
unzip BLAZE_v1.1.0.zip
cd BLAZE
```

## Conda environment
All [dependencies](#dependencies) can be installed using "conda". The configuration file can be found at `conda_env/`:
```
conda config --set channel_priority false
conda env create -f conda_env/environment.yml
conda activate blaze
```

## <a name="dependencies"></a>Dependencies
* `Biopython`
* `pandas`
* `numpy`
* `tqdm`
* `matplotlib`
* `levenshtein`

## Test run
The following command runs BLAZE on a test dataset provided in `/test/data`. The expected output can be found [here](test/).
```
bash test/run_test.sh
```

# Module:

## `blaze.py`: Get putative barcodes and barcode whitelist from fastq files.
This script has been tested on Chromium **Single Cell 3ʹ gene expression v3** and should be able to work on **Chromium Single Cell 3ʹ gene expression v2**, but it doesn't support any 10X 5' gene expression kits.

## **Required Input:** 
 * *Folder of the fastq files*
 * *10X barcode whitelist*: a file containing all possible 10x barcodes ([more details](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)). The 10X barcode whitelists for both the 10X Single Cell 3ʹ gene expression v2 and v3 chemistries have been packed in `10X_bc/`. By default, this module assumes v3 chemistry and uses the file `10X_bc/3M-february-2018.zip`. You may specify `--kit-version=v2` to use the whitelist for v2 chemistry, or provide your own whitelist by specifying `--full-white-list=<filename>`. Note that barcodes outside this whitelist will never be found in the output.

 * *expected number of cells*: In the current version, the expected number of cells is a required input (specify `--expect-cells=xx`). Note that the output is NOT sensitive to the specified number, but a rough number is needed to determine the count threshold to output the whitelist. 

## **Output:**
1. Print when running: Overview stats of the number of reads with and without putative barcodes.
2. Putative barcode in each read, default filename: `putative_bc.csv`. It contains 3 columns:
    * col1: read id
    * col2: putative barcode (i.e. the basecalled barcode segment in each read, specifically the 16nt sequence after the identifed 10X adaptor within each read **without correction for any basecalling errors**)
    * col3: minimum Phred score of the bases in the putative barcode
   
    **Note:** col2 and 3 will be empty if no barcode is found within a read. 
3. Cell-ranger style barcode whitelist, default filename: `whitelist.csv`.
4. "Barcode rank plot"" (or "knee plot") using the high-quality putative barcodes.

## Understanding the BLAZE output
BLAZE follows the 3-step process:
1. **Locate the putative barcode in each read:**
BLAZE first searches for the putative barcode (i.e. barcode in read without error correction) by locating 10X adapter and polyT in the reads. Both the sequenced strand and reverse complement strand are considered for each read. These putative barcodes and their quality scores (minQ) are recorded in `putative_bc.csv`. The IDs of reads where no putative barcode can be located are also listed in the file but without any putative barcode or minQ score.
 
2. **Identify high-quality putative barcodes:**
BLAZE filters the putative barcodes to get a list of high-quality putative barcodes using the following criteria:
   * High-quality putative barcodes should be an exact sequence match for barcodes in the complete 10x whitelist.
   * minQ >= threshold (15 by default).
   
   Note: this step happens internally without an output file.

3. **Generate the barcode whitelist:**
BLAZE scans through the list of high-quality putative barcodes and counts the number of appearances of each unique barcode sequence.  
Finally, BLAZE picks those unique barcodes whose counts are above a quantile-based threshold and writes them into `whitelist.csv`.

Note: With the input information provided, we could assign putative barcodes to the whitelist by analysing their edit distances. In practise we have found this is not required to identify the cell-associated barcodes. In addtion, downstream tools like FLAMES already perform this assignment and so we didn't reimplement it in BLAZE. Examples of running FLAMES using BLAZE's whitelist can be find at https://github.com/youyupei/bc_whitelist_analysis/.

## Additional (optional) features

### High-sensitivity mode
By default, BLAZE is configured to minimise false-positive barcode detections and is therefore relatively conservative. BLAZE has a high-sensitivity mode for users who prioritise high recall of the barcodes (and cells) present. To use specify `--high-sensitivity-mode`. Users should be aware that high-sensitivity mode trades higher recall (i.e. more true barcodes) for potentially lower precision (i.e. more non-cell associated barcodes) and therefore we recommended running an empty drops analysis and removing cells with an ambient profile if using BLAZE HS mode (see below).

### Output a list of barcodes associated to empty droplets
Users may prioritise higher sensitivity by using [high-sensitivity mode](#high-sensitivity-mode) or a user-specified count threshold (run `python3 blaze.py -h` for more details), which may also ouput more non-cell associated barcodes. In such situations we recommended running an empty drops analysis[[1]](#1) to distinguish cell-associated barcodes and barcodes from empty droplets with an ambient RNA expression profile. BLAZE can output a list of barcodes from empty droplets (`empty_bc.csv`) by specifying `--emptydrop`.

**Criteria for selecting these empty droplet barcodes:**

1. Barcodes in 10x full whitelist, and
2. Edit distance > 4 from all selected cell-associated BCs in the barcode whitelist and
3. High-quality putative barcode count < a certain number (specified using `--emptydrop-max-count`, default: infinity)


## **Example code:**

Run BLAZE in default mode: the expected number of cells are set to be 1000 and run with 12 threads
```
python3 blaze.py --expect-cells=1000 --threads=12 path/to/fastq_pass
```
Run BLAZE in high-sensitivity mode: the expected number of cells are set to be 1000 and run with 12 threads
```
python3 blaze.py --high-sensitivity-mode --expect-cells=1000 --threads=12 path/to/fastq_pass
```

**More information:**
```
python3 blaze.py -h
```
**Note:** If you need to change any arguments after running `blaze.py`, you DO NOT have to re-run `blaze.py` but can use `bin/update_whitelist.py`. For example, after running the "Run BLAZE in default mode" code above (`putative_bc.csv` will be generated after running), if you also need output in "high-sensitivity mode" you could run: 
```
python3 update_whitelist.py --expect-cells 1000 --high-sensitivity-mode --emptydrop --out-bc-whitelist=whitelist_hs putative_bc.csv
```
To obtain an output quickly. See [more examples here](test/) or run `python3 bin/update_whitelist.py -h` for more details.


# Citing BLAZE

If you find BLAZE useful for your work, please cite our paper:

>You, Y., Prawer, Y. D., De Paoli-Iseppi, R., Hunt, C. P., Parish, C. L., Shim, H., & Clark, M. B. (2023) Identification of cell barcodes from long-read single-cell RNA-seq with BLAZE. Genome Biol 24, 66.
>[You et al. 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02907-y)


# Data availability
The data underlying the article "Identification of cell barcodes from long-read single-cell RNA-seq with BLAZE" are available from ENA under accession PRJEB54718. The processed data and scripts used in this study are available at https://github.com/youyupei/bc_whitelist_analysis/.



# References
<a id="1">[1]</a> 
Lun, A. T., Riesenfeld, S., Andrews, T., Gomes, T., & Marioni, J. C. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome biology, 20(1), 1-9.
