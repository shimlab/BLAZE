# set up and activate conda environment if you haven't
# conda env create -f conda_env/environment.yml
# conda activate blaze

# search for putative barcode in each read and obtain the whitelist
python3 ../bin/blaze.py --expect-cells=500 --threads=12 data/

# update the whitelist to obtain a whitelist in high-confidece mode 
python3 ../bin/update_whitelist.py --expect-cells 500 --high-sensitivity-mode --emptydrop --out-bc-whitelist=whitelist_hs putative_bc.csv
