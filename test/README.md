## Test run

The following command runs BLAZE on the fastq files stored in `/data`.
```
bash test_run.sh
``` 

## Expected output
Two lines ofcommand were tested and all the expected output files can be found in `/expect_output`:

### Commend 1: `blaze.py`
**Goal:**
Search for putative barcodes in each reads, and generate a whitelist using `blaze.py`. The expected number of cells are set to be 500 and run with 12 threads
   ```
   python3 bin/blaze.py --expect-cells 500 --threads 12 data/
   ```
   
**Output files:**

 * `putative_bc.csv`
 * `whitelist.csv`
 * `knee_plot.png`
   

### Commend 2: `update_whitelist.py`
**Goal:**
Generate a new whitelist using the putative barcode table (`putative_bc.csv` from previous run) without the need of rerunning the putative barcode search. The expected number of cells are set to be 500 and the new whitelist is set to be in "High-sensitivity mode", and a list of barcodes from empty droplets are requested.
   ```
   python3 ../bin/update_whitelist.py --expect-cells 500 --high-sensitivity-mode --emptydrop --out-bc-whitelist=whitelist_hs putative_bc.csv
   ```
**Output files:**

 * `whitelist_hs.csv`
 * `knee_plot_update.png`
 * `emtpy_bc.csv`

**Note:** `update_whitelist.py` saves whitelist as "whitelist.csv" by default, which overwrites the default output "whitelist.csv" from `blaze.py`.
