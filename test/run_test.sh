# run the following command to test the above commands
mkdir test_out
blaze --expect-cells=500 --threads=12 --output-prefix test_out/test_ data/


blaze --expect-cells=10 --threads=1 --full-bc-whitelist /home/youyupei/github_repo/shimlab/BLAZE/blaze/10X_bc/3M-5pgex-jan-2023.txt  --kit-version 5v3 --overwrite  --output-prefix  test_  test02.fastq.gzls
