# run the following command to test the above commands
mkdir test_out
blaze --expect-cells=500 --threads=12 --output-prefix test_out/test_ data/
