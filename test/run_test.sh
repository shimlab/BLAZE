# set up and activate conda environment if you haven't
# conda env create -f conda_env/environment.yml
# conda activate blaze

#python3 ../bin/blaze.py --expect-cells=500 --threads=12 data/

python3 ../bin/update_whitelist.py --count-threshold=50  --out-bc-whitelist=whitelist_update putative_bc.csv
