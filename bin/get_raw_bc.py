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
import os
import getopt
import Bio.SeqIO
from collections import defaultdict, Counter
from tqdm import tqdm
import multiprocessing as mp
import textwrap
from pathlib import Path
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import zipfile
import io
import logging


import helper
from config import *
from polyT_adaptor_finder import Read

def parse_arg():
    
    def print_help():
        help_message=\
            f'''
            Usage: python3 {argv[0]} [OPTIONS] <fastq directory>
            
            Options:
                -h, --help
                    Print this help message.
                
                --expect-cells (required in current version)
                    <INT>:  Expected number of cells. Default: not specified
                
                --minQ:
                    <INT>: Minimum phred score for all bases in a raw BC. Reads whose 
                    raw BC contains one or more bases with Q<minQ is not counted 
                    in the "Raw BC rank plot". Default: --minQ=15
                
                --kit-version:
                    Choose from v2 and v3 (for 10X Single Cell 3ʹ gene expression v2 or v3). 
                    Default: --kit_version=v3.

                --full-bc-whitelist=
                    <path to file>: .txt file containing all the possible BCs. Users may provide
                    their own whitelist. No need to specify this if users want to use the 10X whilelist. 
                    The correct version of 10X whilelist will be determined based on 10X kit version.
                
                --out-raw-bc
                    <filename_prefix>: Output a csv file for the raw BC in each read. 
                                        Default: --out-raw-bc=raw_bc
                
                --out-bc-whitelist
                    <filename_prefix>: Output the whitelist identified from all the reads. 
                                        Default: --out-bc-whitelist=whitelist
                
                --threads
                    <INT>: Number of threads used <default: # of available cpus - 1>
                                        
            '''
        print(textwrap.dedent(help_message))
   
    
    argv = sys.argv
    
    # Default 
    n_process = mp.cpu_count()-1
    exp_cells = None
    full_bc_whitelist = None
    min_phred_score = DEFAULT_GRB_MIN_SCORE
    kit = DEFAULT_GRB_KIT
    out_raw_bc = DEFAULT_GRB_OUT_RAW_BC
    out_whitelist = DEFAULT_GRB_OUT_WHITELIST
    
    # Read from options

    try: 
        opts, args = getopt.getopt(argv[1:],"h",
                    ["help","threads=","minQ=","full-bc-whitelist=",
                     "out-raw-bc=", "out-bc-whitelist=", "expect-cells=", "kit-version="])
    except getopt.GetoptError:
        helper.err_msg("Error: Invalid argument input") 
        print_help()
        sys.exit(1)    
    
    for opt, arg in opts:
        if opt in  ("-h", "--help"):
            print_help()
            sys.exit(0)
        elif opt == '--expect-cells':
            exp_cells = int(arg)
        elif opt == '--threads':
            n_process = int(arg)
        elif opt == '--minQ':
            min_phred_score = int(arg)  
        elif opt == '--full-bc-whitelist':
            full_bc_whitelist = arg  
        elif opt == "--out-raw-bc":
            out_raw_bc = arg      
        elif opt == "--out-bc-whitelist":
            out_whitelist = arg 
        elif opt == "--kit-version":
            kit = arg.lower()


    if kit not in ['v2', 'v3']:
        helper.err_msg("Error: Invalid value of --kit-version, please choose from v3 or v2") 
        sys.exit()

    if full_bc_whitelist:
        helper.warning_msg(textwrap.dedent(
            f'You are using {os.path.basename(full_bc_whitelist)} as the full barcode'\
            'whitelist. Note that the barcodes not listed in the file will never be found.'))
    else:
        if kit == 'v3':
            full_bc_whitelist = DEFAULT_GRB_WHITELIST_V3
        elif kit == 'v2':
            full_bc_whitelist = DEFAULT_GRB_WHITELIST_V2

    

    # Read from args
    if not args:
        helper.err_msg("Error: Missing fastq directory.")   
        print_help()# print help doc when no command line args provided
        sys.exit(0)
    
    # check input
    fastq_dir = args[0]
    if not os.path.isdir(fastq_dir):
        helper.err_msg("Error: Input directory doesn't exist. Note that the input should be a directory instead of file.") 
        sys.exit(1)
    if not exp_cells:
        helper.err_msg("--expect-cells is required to build the whitelist!") 
        sys.exit(1)


    # check file
    helper.check_exist([full_bc_whitelist, fastq_dir])
    return fastq_dir, n_process, exp_cells ,min_phred_score, full_bc_whitelist, out_raw_bc, out_whitelist



# Parse fastq -> polyT_adaptor_finder.Read class
def get_raw_bc_from_fastq(fn, min_q=0):
    """
    Get raw BC from each reads from fastq file

    Parameters
    ----------
    fn : STR
        filename of the fastq.
    min_q: INT
        Only count raw bc with minimum value specified
    save_raw_bc: STR
        Output filename for the raw BCs. Will not output anything by default
    Returns
    -------
    None.

    """

    
    read_list = []
    raw_bc = []
    raw_bc_pass = []
    # raw_bc_count = Counter([])
    # raw_bc_pass_count = Counter([])
    fastq = Bio.SeqIO.parse(fn, "fastq")
    for i,r in enumerate(fastq):
        read = Read(read_id = r.id, sequence=str(r.seq), 
                    phred_score=r.letter_annotations['phred_quality'], fastq_name=fn)        
        read.get_strand_and_raw_bc()
        read_list.append(read)
        if read.raw_bc_min_q and read.raw_bc_min_q >= min_q:     
            raw_bc.append(read.raw_bc)
            
        if read.raw_bc_min_q and read.raw_bc_min_q < min_q:  
            raw_bc_pass.append(100)
        else:
            raw_bc_pass.append(read.adaptor_polyT_pass)
    # raw_bc_count += Counter(raw_bc)
    # raw_bc_pass_count += Counter(raw_bc_pass)
    
    rst_df = pd.DataFrame(
        {'read_id': [r.id for r in read_list],
         'raw_bc': [r.raw_bc for r in read_list],
         'raw_bc_min_q': [r.raw_bc_min_q for r in read_list]
            }
        )
    return Counter(raw_bc), Counter(raw_bc_pass), rst_df


def qc_report(pass_count, min_phred_score):
    '''
    Print stats for 
        a. # of input read in total (only look at primary alignment if it's BAM)
        b. # and % read pass the polyT and adaptor searching process.
        c. # and % with poly T and adaptor find in both ends
        d. # adnd % with no poly T and adpator find
    Parameters
    ----------
    pass_count : Counter object
        count of different type of error
    min_phred_score :  INT
        min score used to filter out raw BC
    '''
    total_read = sum(pass_count.values())
    
    print_message=\
        f'''
        Total number of reads: 
            {total_read:,}
        Reads with unambigous polyT and adapter positions found:            
            {pass_count[0]+ pass_count[100]:,} ({(pass_count[0]+ pass_count[100])/total_read*100:.2f}% of all reads)
            {pass_count[0]:,} in which all bases in the raw BC have Q>={min_phred_score}
        Failed Reads: 
            no polyT and adapter positions found: 
                {pass_count[1]:,} ({pass_count[1]/total_read*100:.2f}% of all reads)
            polyT and adapter positions found in both end (fail to determine strand): 
                {pass_count[2]:,} ({pass_count[2]/total_read*100:.2f}% of all reads)
            multiple polyT and adapter found in one end
                {pass_count[10]:,} ({pass_count[10]/total_read*100:.2f}% of all reads)
        '''
    print(textwrap.dedent(print_message))
    
    
def get_bc_whitelist(raw_bc_count, full_bc_whitelist, exp_cells = None, count_t = None):
    '''
    Get a whitelist from all raw cell bc. If the expect number of cell is provided,
    

    Parameters
    ----------
    raw_bc_count : TYPE
        DESCRIPTION.
    exp_cell : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    list

    '''
    # use the threshold function in config.py
    percentile_count_thres = default_count_threshold_calculation
    
    
    whole_whitelist = []    
    
    if full_bc_whitelist.endswith('.zip'):
        with zipfile.ZipFile(full_bc_whitelist) as zf:
            # check if there is only 1 file
            assert len(zf.namelist()) == 1

            with io.TextIOWrapper(zf.open(zf.namelist()[0]), 
                                                    encoding="utf-8") as f:
                for line in f:
                    whole_whitelist.append(line.strip())
    else:
        with open(full_bc_whitelist, 'r') as f:
            for line in f:
                whole_whitelist.append(line.strip())
    
    whole_whitelist = set(whole_whitelist)
    
    raw_bc_count = {k:v for k,v in raw_bc_count.items() if k in whole_whitelist}

    # plot

    
    # determine real bc based on the count threshold
    if count_t:
        knee_plot(list(raw_bc_count.values()), count_t)
        return {k:v for k,v in raw_bc_count.items() if v > count_t}
    
    elif exp_cells:
        t = percentile_count_thres(list(raw_bc_count.values()), exp_cells)
        knee_plot(list(raw_bc_count.values()), t)
        return {k:v for k,v in raw_bc_count.items() if v > t}
    else: 
        raise ValueError('Invalid value of count_t and/or exp_cells.')
   
def knee_plot(counts, threshold=None):
    counts = sorted(counts)[::-1]
    plt.figure(figsize=(8, 8))
    plt.title(f'Knee plot (from high-confident raw BC)')
    plt.loglog(counts,marker = 'o', linestyle="", alpha = 1, markersize=6)
    plt.xlabel('Barcodes')
    plt.ylabel('Read counts')
    plt.axhline(y=threshold, color='r', linestyle='--', label = 'cell calling threshold')
    plt.legend()
    plt.savefig('knee_plot.png')


 
def main():
    fastq_dir, n_process, exp_cells ,min_phred_score, full_bc_whitelist, out_raw_bc, out_whitelist = parse_arg()
    
    # get raw bc
    fastq_fns = list(Path(fastq_dir).rglob('*.fastq'))
    print(f'Getting raw barcodes from {len(fastq_fns)} FASTQ files...')
    rst_futures = helper.multiprocessing_submit(get_raw_bc_from_fastq,
                                                fastq_fns, n_process=n_process, min_q=min_phred_score)
    
    raw_bc_count = Counter([])
    raw_bc_pass_count = Counter([])    
    rst_dfs = []
    for f in rst_futures:
        count_bc, count_pass, rst_df = f.result()
        raw_bc_count += count_bc
        raw_bc_pass_count += count_pass
        rst_dfs.append(rst_df)
    rst_df = pd.concat(rst_dfs)
    
    print('\nPreparing raw barcode table...')
    if out_raw_bc:
        rst_df.to_csv(out_raw_bc+'.csv', index=False)
    helper.green_msg(f'Raw barcode table saved in {out_raw_bc}.csv')
    
    # output
    print('\n----------------------stats of the raw barcode in reads--------------------------')
    qc_report(raw_bc_pass_count, min_phred_score = min_phred_score)
    print('-----------------------------------------------------\n')

    print("Getting whitelist...\n")
    
    logger = logging.getLogger()

    try:
        bc_whitelist = get_bc_whitelist(raw_bc_count,
                                full_bc_whitelist, 
                                exp_cells=exp_cells)
        with open(out_whitelist+'.csv', 'w') as f:
            for k in bc_whitelist.keys():
                f.write(k+'-1\n')
    except Exception as e:
        logger.exception(e)
        helper.err_msg(
            "Error: Failed to get whitelist. Please check the input files and settings."\
            "Note that the whilelist can be obtained"\
            f"from {out_raw_bc}.csv by using update_whitelist.py."\
            "Run \"python3 BLAZE/bin/update_whitelist.py -h\"  for more details."
            )
    else:
        helper.green_msg(f'Whitelist saved as {out_whitelist}.csv!')
        


if __name__ == '__main__':
    main()
# run multiprocessing


# Parse BAM -> polyT_adaptor_finder.Read class



# # test
# test_data = ['../test_data/FAP25193_pass_c3ecd52d_2789.fastq','../test_data/FAP25193_pass_c3ecd52d_2790.fastq']
# raw_bc = []
# raw_bc_pass = []
# raw_bc_count = Counter([])

# fastq = Bio.SeqIO.parse(test_data, "fastq")
# for i,r in tqdm(enumerate(fastq)):
#     if i < 100000:
#         read = Read(sequence=str(r.seq), strand=r.letter_annotations['phred_quality'])        
#         read.get_strand_and_raw_bc()
#         raw_bc.append(read.raw_bc)
#         raw_bc_pass.append(read.adaptor_polyT_pass)
# raw_bc_count += Counter(raw_bc)


# seq = 'AATCATGCTTCGTTCAATTGCACGTATACCTACTAAAGCCTTCCAGCCTACGACGCTCTTCCGATCTCCAGCGTCGGTCATAATTATAAATAATTTTTTTTTTTTTTTTTTTTTTTGTGACGGAATCTTGTCACCAGGCTGGAGATACAGTGGCATGATCTTGGCTCACTGCAACCTCAGTCTCCTGAGTTCAAGCTTTTGATTCTTCTGTGGCCTCCCGAGTAGCTGGGACTACAGGCATGCGCCACCACGCCTGGCTACTTTTTGTATTTTTAGTAGAGACAGAGTTTCACCCTTATTGGCCACTGGTCTTGAACTCCTGACCTCGTGATCCGCCCACCTCGGGCTCCCAAAGTGCTGGGATTACAGGTATGAGCCACCACCCTGTGACCTTTTCTCAGACATTCTCTCAGCTCCAGCTGACCTGGGGTAGGATTACGGCCTCAAGAAGCAGCCCCCTTTCCATCGGAAGAGCAGAACCTTGGGCCCTGGCAGAGGCGAGGGCCTGAGTGAGACAGGCATTTGGTCCTGGCTCACTGCAAGCTTATAGAGCATTGCCAGAGTCATCTTGAGACCTCTGGGGCTTCGTGGAAGTTTCCTCTGTAGAAATGGAAAAATCTTTCAGATCCGGCCCTCGAGCGCAGCATTCTCTCTACAAAATTGGAGAAAGTTTCAAAATCCTAGGCTCAGGATCTGTAAGGATGTTCTCAAATGCAAAACCCATGATTTGGAAAGTCAGTGTAGAAAAATAGTGTTGGCTGGGCGCAGTGGCTGCTCCTATAATCCCAGCATGGGAGGCCGAGGTGGGCGAGTCACGAGGTCAGGAAATCGAGACCATCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAATACAAAAATTAGCCAGCCGTGGTGGTGCATGCCTGTAGTCCCAGCTGCTCTGGAGGCTAAGGCAGGAGAATCACTTGAAACTGGGAAGGCAGAGGTTGCAGTGAGCCGAGATTACGCCACTGCACTCCAGCCTGGGCGACAGAACGAAATCTGTCTCAAAGGGAAAGAAAGGGAGGAAGAGGAGAGGAGAGGAGGAGACAAGGATGTTGGCTGAGTGTAGTGGCTCACCTTTAATCCCAACACTTTGGGAGCAAAACGGTGGATTTACCTGAGTCGGAGACCCATGTACTCTGCATTGATACCACTGTAGCCATTACGGCCTGTAAAGCAATGCGCT'
# read = Read(seq)
# read.get_strand_and_raw_bc()
# read.raw_bc


