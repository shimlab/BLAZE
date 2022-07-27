## Test run

The following command runs BLAZE on the fastq files stored in `/data`.
```
bash test_run.sh
``` 

## Expected output
There are two modules tested and all the files can be found in `/expect_output`:
1. `blaze.py`(searches for putative barcodes in each reads, and generates whitelist):
    * putative_bc.csv
    * whitelist.csv
    * knee_plot.png
2. `update_whitelist.py` (changes the count threshold and generate a new whitelist, which is useful when users want to tune the threshold without rerunning the putative barcode searching)
    * whitelist_update.csv
    * knee_plot_update.png
    
**Note:** `update_whitelist.py` saves whitelist as "whitelist.csv" by default, which overwrites the default output "whitelist.csv" from `blaze.py`.
