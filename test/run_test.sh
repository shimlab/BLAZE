# run the following command to test the above commands
mkdir test_out
blaze --expect-cells=500 --threads=12 --output-prefix test_out/test_ data/
blaze --expect-cells=1 --threads=12 --kit-version 5v3 --output-prefix test_out/test_5prim 5prim_test.fastq.gz