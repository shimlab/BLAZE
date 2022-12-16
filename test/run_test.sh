# set up and activate conda environment if you haven't
# conda env create -f conda_env/environment.yml
# conda activate blaze

# search for putative barcode in each read and obtain the whitelist
python3 ../bin/blaze.py --expect-cells=500 --threads=12 data/


# update the whitelist use a customised count threhold 
python3 ../bin/update_whitelist.py --count-threshold=50  --out-bc-whitelist=whitelist_update putative_bc.csv

# update the whitelist use a customised count threhold
python3 ../bin/update_whitelist.py --count-threshold=50  --out-bc-whitelist=whitelist_high_sensitivity --high-sensitivity-mode putative_bc.csv
