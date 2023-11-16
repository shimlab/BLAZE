import unittest
import filecmp
from blaze import blaze
import blaze.helper as helper
from blaze.config import *
import os
import gzip
import shutil
import difflib

class TestBlazeMain(unittest.TestCase):

    def setUp(self):
        self.argv = \
            ' --expect-cells=500 --threads=8 --output-prefix test_ test/data/'
        self.expected_normal_files = [
            'test_' + DEFAULT_GRB_OUT_RAW_BC,
            'test_' + DEFAULT_GRB_OUT_WHITELIST,
            'test_' + DEFAULT_EMPTY_DROP_FN,
            'test_' + DEFAULT_KNEE_PLOT_FN,
            'test_' + DEFAULT_BC_STAT_FN
        ]
        self.expected_gz_files = ['test_' + DEFAULT_GRB_OUT_FASTQ]

        self.expected_dir = 'test/expect_output/'

    def decompress_gz_file(self, file_path):
        with gzip.open(file_path, 'rb') as f_in:
            with open(file_path+'.tmp', 'wb') as f_out:  # remove .gz extension
                shutil.copyfileobj(f_in, f_out)
        return file_path+'.tmp'

    def compare_gz_files(self, file1, file2):
        file1_decompressed = self.decompress_gz_file(file1)
        file2_decompressed = self.decompress_gz_file(file2)
        comparison_result = filecmp.cmp(file1_decompressed, file2_decompressed)
        
        # Clean up decompressed files
        os.remove(file1_decompressed)
        os.remove(file2_decompressed)
        if comparison_result:
            return comparison_result
        else:
            self.diff_files(file1_decompressed, file2_decompressed)
            return comparison_result
    
    def diff_files(self, file1, file2):
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            diff = difflib.unified_diff(
                f1.readlines(),  # Read lines from the first file
                f2.readlines(),  # Read lines from the second file
                fromfile=file1,
                tofile=file2,
            )
            # Output the difference to console or write to a file
            for i, line in enumerate(diff):
                if i > 10:
                    break
                print(line, end='')

    def test_main(self):
        print([str(x) for x in helper.get_files('test/data/', ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz'])])
        # Run the main function
        blaze(self.argv)

        # Check that the expected files were created and are correct
        for fn in self.expected_normal_files:
            self.assertTrue(os.path.exists(fn), f"File {fn} was not created")
            if not filecmp.cmp(fn, self.expected_dir + fn):
                with open(fn, 'r') as f1, open(self.expected_dir + fn, 'r') as f2:
                    for i in range(10):
                        print(fn, '\t' ,f1.readline())
                        print(self.expected_dir + fn, '\t', f2.readline())
                self.diff_files(fn, self.expected_dir + fn)
            self.assertTrue(filecmp.cmp(fn, self.expected_dir + fn), f"File {fn} does not match the expected output (({self.expected_dir + fn}))")
            
         
        for fn in self.expected_gz_files:
            self.assertTrue(os.path.exists(fn), f"File {fn} was not created")
            self.assertTrue(self.compare_gz_files(fn, self.expected_dir + fn) , f"File {fn} does not match the expected output ({self.expected_dir + fn})")

    def tearDown(self):
        # Remove created files after tests
        for fn in self.expected_normal_files + self.expected_gz_files:
            if os.path.exists(fn):
                os.remove(fn)

if __name__ == '__main__':
    unittest.main()