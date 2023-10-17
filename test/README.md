## Test run

The following command runs BLAZE on the fastq files stored in `/data`.
```
bash test_run.sh
``` 

## Expected output
Two lines of command were tested and all the expected output files can be found in `/expect_output`:

**Goal:**
Search for putative barcodes in each read,  generate a barcode list and assign reads to the barcodes . The expected number of cells is set to 500 and run with 12 threads
```
blaze --expect-cells=500 --threads=12 --output-prefix test_out/test_ data/
```
   