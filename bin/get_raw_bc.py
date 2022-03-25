"""
Get raw BC from either fastq file or bam file

Output:
    1. csv file a the raw BC count
    2. fastq or BAM with tagged raw BC start position
    3. Stats report (stdout):
        a. how many input read in total (only look at primary alignment if BAM provided)
        b. how many (# and %) read pass the polyT and adaptor searching process.
        c. # and % with poly T and adaptor find in both ends
        d. # adnd % with no poly T and adpator find
        e. # of BC identified
        f. # of reads exactly match the BC
        g. # of reads approximately match the BC with ED = 1
        h. # of reads approximately match the BC with ED = 2
        i. # of reads approximately match the BC with ED = 3
        j. # of reads approximately match the BC with ED > 4
    4. BC rank plot
    5. csv file of the identified cell-associated BC
"""



# step 1:
    ## get the raw BC count
    
    
import sys
import getopt
import Bio.SeqIO
from collections import defaultdict, Counter
from tqdm import tqdm
import multiprocessing as mp
import textwrap
from pathlib import Path


import helper
from polyT_adaptor_finder import Read

def parse_arg():
    
    def print_help():
        help_message=\
            f'''
            Usage: python3 {argv[0]} [OPTIONS] <fastq directory>
            
            Options:
                -h, --help
                    Print this help message.
                --threads
                    Number of threads used <default: # of available cpus - 1>
            '''
        print(textwrap.dedent(help_message))
   
    
    argv = sys.argv
    
    # Default 
    n_process = mp.cpu_count()-1
    
    # Read from options
    try: 
        opts, args = getopt.getopt(argv[1:],"h",
                    ["help","threads="])
    except getopt.GetoptError:
        helper.err_msg("Error: Invalid argument input") 
        print_help()
        sys.exit(1)    
    
    for opt, arg in opts:
        if opt in  ("-h", "--help"):
            print_help()
            sys.exit(0)
        elif opt == '--threads':
            n_process = arg
            
    # Read from args
    if not args:
        helper.err_msg("Error: Missing fastq directory.")   
        print_help()# print help doc when no command line args provided
        sys.exit(0)
    
    fastq_dir = args[0]
    
    return fastq_dir, n_process



# Parse fastq -> polyT_adaptor_finder.Read class
def get_raw_bc_from_fastq(fn):
    """
    Get raw BC from each reads from fastq file

    Parameters
    ----------
    fn : STR
        filename of the fastq.

    Returns
    -------
    None.

    """
    raw_bc = []
    raw_bc_pass = []
    # raw_bc_count = Counter([])
    # raw_bc_pass_count = Counter([])

    fastq = Bio.SeqIO.parse(test_data, "fastq")
    for i,r in enumerate(fastq):
        read = Read(sequence=str(r.seq), strand=r.letter_annotations['phred_quality'])        
        read.get_strand_and_raw_bc()
        raw_bc.append(read.raw_bc)
        raw_bc_pass.append(read.adaptor_polyT_pass)
    # raw_bc_count += Counter(raw_bc)
    # raw_bc_pass_count += Counter(raw_bc_pass)
    
    return Counter(raw_bc),  Counter(raw_bc_pass)


    


def main():
    fastq_dir, n_process = parse_arg()
    fastq_fns = list(Path(fastq_dir).rglob('*.fastq'))
    print(fastq_fns)
    exit()
    rst_futures = helper.multiprocessing_submit(get_raw_bc_from_fastq ,fastq_fns)
    
    raw_bc_count = Counter([])
    raw_bc_pass_count = Counter([])    
    
    for f in rst_futures:
        count_bc, count_pass = f.result()
        raw_bc_count += count_bc
        raw_bc_pass_count += count_pass
    print(raw_bc_pass_count)

if __name__ == '__main__':
    main()
# run multiprocessing


# Parse BAM -> polyT_adaptor_finder.Read class



# test
test_data = ['../test_data/FAP25193_pass_c3ecd52d_2789.fastq','../test_data/FAP25193_pass_c3ecd52d_2790.fastq']
raw_bc = []
raw_bc_pass = []
raw_bc_count = Counter([])

fastq = Bio.SeqIO.parse(test_data, "fastq")
for i,r in tqdm(enumerate(fastq)):
    if i < 100000:
        read = Read(sequence=str(r.seq), strand=r.letter_annotations['phred_quality'])        
        read.get_strand_and_raw_bc()
        raw_bc.append(read.raw_bc)
        raw_bc_pass.append(read.adaptor_polyT_pass)
raw_bc_count += Counter(raw_bc)


seq = 'AATCATGCTTCGTTCAATTGCACGTATACCTACTAAAGCCTTCCAGCCTACGACGCTCTTCCGATCTCCAGCGTCGGTCATAATTATAAATAATTTTTTTTTTTTTTTTTTTTTTTGTGACGGAATCTTGTCACCAGGCTGGAGATACAGTGGCATGATCTTGGCTCACTGCAACCTCAGTCTCCTGAGTTCAAGCTTTTGATTCTTCTGTGGCCTCCCGAGTAGCTGGGACTACAGGCATGCGCCACCACGCCTGGCTACTTTTTGTATTTTTAGTAGAGACAGAGTTTCACCCTTATTGGCCACTGGTCTTGAACTCCTGACCTCGTGATCCGCCCACCTCGGGCTCCCAAAGTGCTGGGATTACAGGTATGAGCCACCACCCTGTGACCTTTTCTCAGACATTCTCTCAGCTCCAGCTGACCTGGGGTAGGATTACGGCCTCAAGAAGCAGCCCCCTTTCCATCGGAAGAGCAGAACCTTGGGCCCTGGCAGAGGCGAGGGCCTGAGTGAGACAGGCATTTGGTCCTGGCTCACTGCAAGCTTATAGAGCATTGCCAGAGTCATCTTGAGACCTCTGGGGCTTCGTGGAAGTTTCCTCTGTAGAAATGGAAAAATCTTTCAGATCCGGCCCTCGAGCGCAGCATTCTCTCTACAAAATTGGAGAAAGTTTCAAAATCCTAGGCTCAGGATCTGTAAGGATGTTCTCAAATGCAAAACCCATGATTTGGAAAGTCAGTGTAGAAAAATAGTGTTGGCTGGGCGCAGTGGCTGCTCCTATAATCCCAGCATGGGAGGCCGAGGTGGGCGAGTCACGAGGTCAGGAAATCGAGACCATCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAATACAAAAATTAGCCAGCCGTGGTGGTGCATGCCTGTAGTCCCAGCTGCTCTGGAGGCTAAGGCAGGAGAATCACTTGAAACTGGGAAGGCAGAGGTTGCAGTGAGCCGAGATTACGCCACTGCACTCCAGCCTGGGCGACAGAACGAAATCTGTCTCAAAGGGAAAGAAAGGGAGGAAGAGGAGAGGAGAGGAGGAGACAAGGATGTTGGCTGAGTGTAGTGGCTCACCTTTAATCCCAACACTTTGGGAGCAAAACGGTGGATTTACCTGAGTCGGAGACCCATGTACTCTGCATTGATACCACTGTAGCCATTACGGCCTGTAAAGCAATGCGCT'
read = Read(seq)
read.get_strand_and_raw_bc()
read.raw_bc


