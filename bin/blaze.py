"""
Get raw BC (i.e. putative barcode) from fastq file

Output:
    1. csv file a the raw BC count
    2. Stats report (stdout):
        a. how many input read in total (only look at primary alignment if BAM provided)
        b. how many (# and %) read pass the polyT and adaptor searching process.
        c. # and % with poly T and adaptor find in both ends
        d. # adnd % with no poly T and adpator find
    3. BC rank plot
    4. csv file of the BC whitelist
"""
  
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
import gzip


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
                    <INT>: Minimum phred score for all bases in a putative BC. Reads whose 
                    putative BC contains one or more bases with Q<minQ is not counted 
                    in the "Putative BC rank plot". Default: --minQ=15

                --high-sensitivity-mode:
                    Turn on the sensitivity mode, which increases the sensitivity of barcode
                    detections but potentially increase the number false/uninformative BC in
                    the whitelist.

                --kit-version:
                    Choose from v2 and v3 (for 10X Single Cell 3ʹ gene expression v2 or v3). 
                    Default: --kit_version=v3.

                --full-bc-whitelist=
                    <path to file>: .txt file containing all the possible BCs. Users may provide
                    their own whitelist. No need to specify this if users want to use the 10X whilelist. 
                    The correct version of 10X whilelist will be determined based on 10X kit version.
                
                --out-putative-bc
                    <filename_prefix>: Output a csv file for the putative BC in each read. 
                                        Default: --out-putative-bc=putative_bc
                
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
    high_sensitivity_mode = False

    
    # Read from options

    try: 
        opts, args = getopt.getopt(argv[1:],"h",
                    ["help","threads=","minQ=","full-bc-whitelist=","high-sensitivity-mode",
                     "out-putative-bc=", "out-bc-whitelist=", "expect-cells=", "kit-version="])
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
        elif opt == "--out-putative-bc":
            out_raw_bc = arg      
        elif opt == "--out-bc-whitelist":
            out_whitelist = arg 
        elif opt == "--kit-version":
            kit = arg.lower()
        elif opt == "--high-sensitivity-mode":
            high_sensitivity_mode = True


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
    return fastq_dir, n_process, exp_cells ,min_phred_score, \
            full_bc_whitelist, out_raw_bc, out_whitelist, high_sensitivity_mode

# Parse fastq -> polyT_adaptor_finder.Read class
def get_raw_bc_from_reads(reads, min_q=0):
    """
    Get putative BC from each reads from a single read batch defined in batch_iterator function

    Parameters
    ----------
    reads : LIST
        list of reads entry in Bio.SeqIO.parse
    min_q: INT
        Only count putative bc with minimum value specified
    save_putative_bc: STR
        Output filename for the putative BCs. Will not output anything by default
    Returns
    -------
    1. Counter of high-confidence putative BC
    2. Counter of high-confidence putative BC
    3. pd.DataFrame containing all putative BCs 
    """
    
    read_ids = []
    putative_bcs = []
    putative_bc_min_qs = []
    raw_bc = []
    raw_bc_pass = []

    for i,r in enumerate(reads):
        
        read = Read(read_id = r.id, sequence=str(r.seq), 
                    phred_score=r.letter_annotations['phred_quality'])    
          
        read.get_strand_and_raw_bc()
        
        read_ids.append(read.id)
        
        putative_bcs.append(read.raw_bc)
        
        putative_bc_min_qs.append(read.raw_bc_min_q)

        
        if read.raw_bc_min_q and read.raw_bc_min_q >= min_q:     
            raw_bc.append(read.raw_bc)
            
        if read.raw_bc_min_q and read.raw_bc_min_q < min_q:  
            raw_bc_pass.append(100)
        else:
            raw_bc_pass.append(read.adaptor_polyT_pass)
    # raw_bc_count += Counter(raw_bc)
    # raw_bc_pass_count += Counter(raw_bc_pass)
    
    rst_df = pd.DataFrame(
        {'read_id': read_ids,
         'putative_bc': putative_bcs,
         'putative_bc_min_q': putative_bc_min_qs
        }
        )
    return Counter(raw_bc), Counter(raw_bc_pass), rst_df

def qc_report(pass_count, min_phred_score):
    '''
    Generate report for the putative barcode detection.
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
        min score used to filter out putative BC
    '''
    total_read = sum(pass_count.values())
    
    print_message=\
        f'''
        Total number of reads: 
            {total_read:,}
        Reads with unambigous polyT and adapter positions found:            
            {pass_count[0]+ pass_count[100]:,} ({(pass_count[0]+ pass_count[100])/total_read*100:.2f}% of all reads)
            {pass_count[0]:,} in which all bases in the putative BC have Q>={min_phred_score}
        Failed Reads: 
            no polyT and adapter positions found: 
                {pass_count[1]:,} ({pass_count[1]/total_read*100:.2f}% of all reads)
            polyT and adapter positions found in both end (fail to determine strand): 
                {pass_count[2]:,} ({pass_count[2]/total_read*100:.2f}% of all reads)
            multiple polyT and adapter found in one end
                {pass_count[10]:,} ({pass_count[10]/total_read*100:.2f}% of all reads)
        '''
    print(textwrap.dedent(print_message))
    
    
def get_bc_whitelist(raw_bc_count, full_bc_whitelist, exp_cells=None, count_t=None, high_sensitivity_mode=False):
    """    
    Get a whitelist from all putative cell bc with high-confidence putative bc counts. 
    If the expect number is provided (default), a quantile-based threshold will be 
    calculated to determine the exact cells to be output. Otherwise, a user-specified 
    ount threshold will be used and the cells/Barocdes with counts above the threshold will be output.
    

    Args:
        raw_bc_count (Counter): high-confidence putative BC counts for each unique BC
        full_bc_whitelist (str): filename of the 10X whitelist
        exp_cells (int, optional): expected number of cell. Defaults to None.
        count_t (int, optional): count threshold. Defaults to None.

    Raises:
        ValueError: No valid exp_cells or count_t specified

    Returns:
        dict: 
            key: barcodes selected
            values: high-confidence putative BC counts
    """

    # use the threshold function in config.py
    if high_sensitivity_mode:
        percentile_count_thres = high_sensitivity_threshold_calculation
    else:
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
    """
    Plot knee plot using the high-confidence putative BC counts

    Args:
        counts (list): high-confidence putative BC countstion_
        threshold (int, optional): a line to show the count threshold. Defaults to None.
    """
    counts = sorted(counts)[::-1]
    plt.figure(figsize=(8, 8))
    plt.title(f'Barcode rank plot (from high-quality putative BC)')
    plt.loglog(counts,marker = 'o', linestyle="", alpha = 1, markersize=6)
    plt.xlabel('Barcodes')
    plt.ylabel('Read counts')
    plt.axhline(y=threshold, color='r', linestyle='--', label = 'cell calling threshold')
    plt.legend()
    
    out_fn = 'knee_plot'
    while os.path.exists(out_fn + '.png'):
        out_fn += '_update'
    plt.savefig(out_fn + '.png')


def read_batch_generator(fastq_fns, batch_size):
    """Generator of barches of reads from list of fastq files

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq = Bio.SeqIO.parse(handle, "fastq")
                read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            fastq = Bio.SeqIO.parse(fn, "fastq")
            read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
            for batch in read_batch:
                yield batch


def main():
    
    fastq_dir, n_process, exp_cells ,min_phred_score, full_bc_whitelist,\
        out_raw_bc, out_whitelist, high_sensitivity_mode= parse_arg()
    
    # get raw bc
    fastq_fns = helper.get_files(fastq_dir, ['*.fastq', '*.fq', '*.fastq.gz', '*.fg.gz'])
    print(f'Getting putative barcodes from {len(fastq_fns)} FASTQ files...')

    read_batchs = read_batch_generator(fastq_fns, batch_size=500)

    rst_futures = helper.multiprocessing_submit(get_raw_bc_from_reads,
                                            read_batchs, n_process=n_process, min_q=min_phred_score)

    raw_bc_count = Counter([])
    raw_bc_pass_count = Counter([])    
    rst_dfs = []
    for idx, f in enumerate(rst_futures):
        count_bc, count_pass, rst_df = f.result() #write rst_df out
        raw_bc_count += count_bc
        raw_bc_pass_count += count_pass
        if idx == 0:
            rst_df.to_csv(out_raw_bc+'.csv', index=False)
        else:
            rst_df.to_csv(out_raw_bc+'.csv', mode='a', index=False, header=False)

    helper.green_msg(f'Putative barcode table saved in {out_raw_bc}.csv')
    
    # output
    print('\n----------------------stats of the putative barcodes--------------------------')
    qc_report(raw_bc_pass_count, min_phred_score = min_phred_score)
    print('-----------------------------------------------------\n')

    print("Getting whitelist...\n")
    
    logger = logging.getLogger()

    try:
        bc_whitelist = get_bc_whitelist(raw_bc_count,
                                full_bc_whitelist, 
                                exp_cells=exp_cells,
                                high_sensitivity_mode=high_sensitivity_mode)
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
