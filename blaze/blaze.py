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
import shlex
import os
import getopt
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
from fast_edit_distance import edit_distance
import logging

import blaze.helper as helper
from blaze.config import *
import blaze.polyT_adaptor_finder as polyT_adaptor_finder
import blaze.read_assignment as read_assignment
from datetime import datetime

# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def parse_arg(argv):
    
    def print_help():
        help_message=\
            f'''
            Description:
                BLAZE2 is a tool for demultiplexing 10X single cell long-read RNA-seq data.
                It takes fastq files as input and output a whitelist of barcodes and a fastq 
                with demultiplexed reads.

            Usage: blaze  --expect-cells <INT> [OPTIONS] <fastq directory>

            Required argument:
                --expect-cells <INT> (required)
                            Expected number of cells.

            Options:
                -h, --help
                    Print this help message.

                --output-prefix <prefix>
                    Filename of output files. Default: --output-prefix {DEFAULT_PREFIX}
                    Note that the output can be directed to a different directory by specifying
                    the path in the prefix. E.g., --output-prefix /path/to/output/prefix
                
                --output-fastq <fastq filename>
                    Filename of output fastq file name. Default: --output-fastq {DEFAULT_GRB_OUT_FASTQ}
                    Note that if the filename has to end with .fastq, .fq, .fastq.gz or .fq.gz.

                --no-demultiplexing:
                    Only output the whitelist and do not perform the demultiplexing step.
                
                --max-edit-distance <INT>
                    Maximum edit distance allowed between a putative barcode and a barcode 
                    for a read/putative barcdoe to be assigned to the barcode. Default: --max-edit-distance {DEFAULT_ASSIGNMENT_ED}
                
                --kit-version <v2 or v3>:
                    Choose from 10X Single Cell 3ʹ gene expression v2 or v3. 
                    Default: --kit_version v3.

                --minQ <INT>:
                    Putative BC contains one or more bases with Q<minQ is not counted 
                    in the "Putative BC rank plot". Default: --minQ=15
                
                --full-bc-whitelist <path to file>
                    txt file containing all the possible BCs. You may provide your own whitelist. 
                    No need to specify this if users want to use the 10X whilelist. The correct 
                    version of 10X whilelist will be determined based on 10X kit version.

                --overwrite
                    Overwrite the old file when the output file(s) exist. If not specified, 
                    the steps generating the existing file(s) will be skipped.

                --threads <INT>
                    <INT>: Number of threads used <default: # of available cpus - 1>

                --batch-size <INT>
                    <INT>: Number of reads this program process together as a batch. Not that if 
                    the specified number larger than the number of reads in each fastq files, the 
                    batch size will be forced to be the number of reads in the file. <Default: 1000>

                --minimal_stdout
                    Minimise the command-line printing

            High sensitivity mode:
                --high-sensitivity-mode:
                    Turn on the sensitivity mode, which increases the sensitivity of barcode
                    detections but potentially increase the number false/uninformative BC in
                    the whitelist. 
                    Note that --emptydrop is recommanded specified with this mode (See details below).

            Empty droplet BCs
                --emptydrop-max-count <INT>
                    Only select empty droplet barcodes supported by a maximum number of high-confidence
                    putative barcode count. (Default: Inf, i.e. no maximum number is set
                    and any barcodes with ED>= {DEFAULT_EMPTY_DROP_MIN_ED} to the barcodes
                    in whitelist can be selected as empty droplets)

            '''
        print(textwrap.dedent(help_message))
   
    if not argv:
        argv = sys.argv
    else:
        argv = ['blaze.py'] + shlex.split(argv)
    
    
    # Default 
    out_fastq_fn = DEFAULT_GRB_OUT_FASTQ
    prefix = DEFAULT_PREFIX
    overwrite = False
    n_process = mp.cpu_count()-1
    exp_cells = None
    full_bc_whitelist = None
    min_phred_score = DEFAULT_GRB_MIN_SCORE
    kit = DEFAULT_GRB_KIT
    batch_size=1000
    high_sensitivity_mode = False   
    emptydrop_max_count = np.inf
    do_demultiplexing = True
    max_edit_distance = DEFAULT_ASSIGNMENT_ED
    minimal_out = False

    # Read from options
    try: 
        opts, args = getopt.getopt(argv[1:],"h",
                    ["help","threads=","minQ=","full-bc-whitelist=","high-sensitivity-mode",
                     "output-prefix=", "expect-cells=", "overwrite",
                     "kit-version=", "batch-size=", "emptydrop-max-count=", "output-fastq=",
                     "no-demultiplexing", "max-edit-distance=", "minimal_stdout"])
    except getopt.GetoptError:
        helper.err_msg("Error: Invalid argument input") 
        print_help()
        sys.exit(1)    
    
    for opt, arg in opts:
        if opt in  ("-h", "--help"):
            print_help()
            sys.exit(0)
        elif opt == '--output-prefix':
            prefix = arg
        elif opt == '--expect-cells':
            exp_cells = int(arg)
        elif opt == '--overwrite':
            overwrite = True
        elif opt == '--output-fastq':
            out_fastq_fn = arg
        elif opt == '--threads':
            n_process = int(arg)
        elif opt == '--minQ':
            min_phred_score = int(arg)  
        elif opt == '--full-bc-whitelist':
            full_bc_whitelist = arg  
        elif opt == "--out-putative-bc":
            out_raw_bc_fn = arg      
        elif opt == "--kit-version":
            kit = arg.lower()
        elif opt == "--high-sensitivity-mode":
            high_sensitivity_mode = True
            #emptydrop = True
        elif opt == "--batch-size":
            batch_size = int(arg)
        elif opt == "--emptydrop-max-count":
            emptydrop_max_count = int(arg)
        elif opt == "--no-demultiplexing":
            do_demultiplexing = False
        elif opt == "--max-edit-distance":
            max_edit_distance = int(arg)
        elif opt == '--minimal_stdout':
            minimal_out = True
    # output filename:
    out_fastq_fn = prefix + out_fastq_fn
    out_raw_bc_fn = prefix + DEFAULT_GRB_OUT_RAW_BC
    out_whitelist_fn = prefix + DEFAULT_GRB_OUT_WHITELIST
    out_emptydrop_fn = prefix + DEFAULT_EMPTY_DROP_FN
    out_plot_fn = prefix + DEFAULT_KNEE_PLOT_FN
    summary_fn = prefix + DEFAULT_BC_STAT_FN

    if kit not in ['v2', 'v3']:
        helper.err_msg("Error: Invalid value of --kit-version, please choose from v3 or v2") 
        sys.exit()

    if full_bc_whitelist:
        helper.warning_msg(textwrap.dedent(
            f'You are using {os.path.basename(full_bc_whitelist)} as the full barcode '\
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

    # check required argument
    if not exp_cells:
        helper.err_msg("--expect-cells is required to build the whitelist!") 
        sys.exit(1)
    
    # check input
    fastq_dir = args[0]
    helper.check_exist([full_bc_whitelist, fastq_dir])
    
    if os.path.isdir(fastq_dir):
        fastq_fns = helper.get_files(fastq_dir, ['*.fastq', '*.fq', '*.fastq.gz', '*.fg.gz'])
    elif os.path.isfile(fastq_dir):
        fastq_fns = [fastq_dir]
    else:
        helper.err_msg(f"File type of input file/dir {fastq_fns} is not supported.")
        sys.exit(1)

    # check output filenames
    # files_status = {
    #     out_raw_bc_fn: os.path.exists(out_raw_bc_fn), # False for not existing
    #     out_whitelist_fn: os.path.exists(out_whitelist_fn),
    #     out_fastq_fn: os.path.exists(out_fastq_fn)
    # }

    if not helper.check_suffix(out_fastq_fn, ['.fastq', '.fq', '.fastq.gz', '.fg.gz']):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{out_fastq_fn}' should end with '.fastq', '.fq', '.fastq.gz' or '.fg.gz'. Please check the config.py file in BLAZE.")
        sys.exit(1)

    if not helper.check_suffix(out_raw_bc_fn, '.csv'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{out_raw_bc_fn}' should end with '.csv'. Please check the config.py file in BLAZE.")
        sys.exit(1)

    if not helper.check_suffix(out_whitelist_fn, '.csv'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{out_whitelist_fn}' should end with '.csv'. Please check the config.py file in BLAZE.")
        sys.exit(1)

    if not helper.check_suffix(out_emptydrop_fn, '.csv'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{out_emptydrop_fn}' should end with '.csv'. Please check the config.py file in BLAZE.")
        sys.exit(1)
    if not helper.check_suffix(out_plot_fn, '.png'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{out_plot_fn}' should end with '.png'. Please check the config.py file in BLAZE.")
        sys.exit(1)
    
    # helper.warning_msg(f"Filename {out_fastq_fn} exists, will output demultiplexed reads into {out_fastq_fn}.")

    return fastq_fns, out_fastq_fn, n_process, exp_cells ,min_phred_score, \
            full_bc_whitelist, out_raw_bc_fn, out_whitelist_fn, \
            high_sensitivity_mode, batch_size, out_emptydrop_fn, emptydrop_max_count, \
            overwrite, out_plot_fn, do_demultiplexing, max_edit_distance, summary_fn,\
            minimal_out

# Parse fastq -> polyT_adaptor_finder.Read class
def get_raw_bc_from_reads(reads, min_q=0):
    """
    Get putative BC from each reads from a batch of read (can be defined by batch_iterator function)


    Parameters
    ----------
    reads : LIST
        list of reads named tuple with names: id, seq, q_letter
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
    # init
    read_ids = []
    putative_bcs = []
    putative_bc_min_qs = []
    raw_bc = []
    raw_bc_pass = []
    umis = []
    trim_idxs = []
    pre_bc_flankings = []
    post_umi_flankings = []

    for i,r in enumerate(reads):
        
        # create read object 
        read = polyT_adaptor_finder.Read(read_id = r.id, sequence=str(r.seq), 
                    phred_score=r.q_letter)    
        

        read.get_strand_and_raw_bc()
        read_ids.append(read.id)
        putative_bcs.append(read.raw_bc)
        putative_bc_min_qs.append(read.raw_bc_min_q)
        umis.append(read.putative_UMI)
        trim_idxs.append(read.adator_trimming_idx)
        pre_bc_flankings.append(read.pre_bc_flanking)
        post_umi_flankings.append(read.post_umi_flanking)
        
        if read.raw_bc_min_q and read.raw_bc_min_q >= min_q:     
            raw_bc.append(read.raw_bc)
            
        if read.raw_bc_min_q and read.raw_bc_min_q < min_q:  
            raw_bc_pass.append(100) #tag for low quality putative bc
        else:
            raw_bc_pass.append(read.adaptor_polyT_pass)

    rst_df = pd.DataFrame(
        {'read_id': read_ids,
         'putative_bc': putative_bcs,
         'putative_bc_min_q': putative_bc_min_qs,
        'putative_umi': umis,
        'umi_end': trim_idxs,
        'pre_bc_flanking': pre_bc_flankings,
        'post_umi_flanking': post_umi_flankings
        }
        )
    return Counter(raw_bc), Counter(raw_bc_pass), rst_df

def qc_report(pass_count, min_phred_score, stdout=True, out_fn=None):
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
    
    print_message=textwrap.dedent(
        f'''
        \n----------------------stats of the putative barcodes--------------------------

        Total number of reads: 
            {total_read:,}
        Reads with unambiguous polyT and adapter positions found:            
            {pass_count[0]+ pass_count[100]:,} ({(pass_count[0]+ pass_count[100])/total_read*100:.2f}% of all reads)
            {pass_count[0]:,} in which all bases in the putative BC have Q>={min_phred_score}
        Failed Reads: 
            no polyT and adapter positions found: 
                {pass_count[1]:,} ({pass_count[1]/total_read*100:.2f}% of all reads)
            polyT and adapter positions found in both end (fail to determine strand): 
                {pass_count[2]:,} ({pass_count[2]/total_read*100:.2f}% of all reads)
            multiple polyT and adapter found in one end
                {pass_count[10]:,} ({pass_count[10]/total_read*100:.2f}% of all reads)
        -------------------------------------------------------------------------------\n'
        ''')
    if stdout:
        print(print_message)
    if out_fn:
        with open(out_fn, 'w') as f:
            f.write(print_message + '\n')

    
def get_bc_whitelist(raw_bc_count, full_bc_whitelist, exp_cells=None, 
                    count_t=None, high_sensitivity_mode=False, 
                    output_empty = False, empty_max_count = np.inf, 
                    out_plot_fn = DEFAULT_KNEE_PLOT_FN):
    f"""    
    Get a whitelist from all putative cell bc with high-confidence putative bc counts. 
    If the expect number is provided (default), a quantile-based threshold will be 
    calculated to determine the exact cells to be output. Otherwise, a user-specified 
    ount threshold will be used and the cells/Barocdes with counts above the threshold will be output.
    
    If high_sensitivity_mode = True, the high sensitivity (HS) mode is turned on which uses
    more relaxed threshold  

    If in output_empty=True, a list of BCs that are most likely corresponding to 
    empty droplets will also be produced autimatically , which might be useful in 
    downstream analysis.
        Criteria of selecting these BC:
            1. BC in 10x full whitelist, and
            2. At least {DEFAULT_EMPTY_DROP_MIN_ED} away from all selected cell-associated BCs in whitelist. and
            3. BC selected as empty droplets should below a certain high-confidence putative barcode count (empty_max_count)

    Args:
        raw_bc_count (Counter): high-confidence putative BC counts for each unique BC
        full_bc_whitelist (str): filename of the 10X whitelist
        exp_cells (int, optional): expected number of cell. Defaults to None.
        count_t (int, optional): count threshold. Defaults to None.

    Raises:
        ValueError: No valid exp_cells or count_t specified

    Returns:
        dict 1: 
            key: barcodes selected
            values: high-confidence putative BC counts
       
        list  (high_sensitivity_mode only):
            barcodes most likely associated to empty droplet
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
        knee_plot(list(raw_bc_count.values()), count_t, out_plot_fn)
        cells_bc = {k:v for k,v in raw_bc_count.items() if v > count_t}

        if not output_empty:
            return cells_bc, []
        else:
            # create a confident list of empty drops in high sensitivity mode
            logger.info("Creating emtpy droplets barocde list...")
            ept_bc = []
            ept_bc_max_count = min(cells_bc.values())
            ept_bc_max_count = min(ept_bc_max_count, empty_max_count)
            ept_bc_candidate = \
                [k for k,v in raw_bc_count.items() if v < ept_bc_max_count]
            for k in ept_bc_candidate:
                if min([edit_distance(k, x, max_ed = DEFAULT_EMPTY_DROP_MIN_ED) for x in cells_bc.keys()]) >= DEFAULT_EMPTY_DROP_MIN_ED:
                    ept_bc.append(k)

                # we don't need too much BC in this list
                if len(ept_bc) >  DEFAULT_EMPTY_DROP_NUM:
                    break      
            return cells_bc, ept_bc

    elif exp_cells:
        t = percentile_count_thres(list(raw_bc_count.values()), exp_cells)
        knee_plot(list(raw_bc_count.values()), t, out_plot_fn)

        cells_bc = {k:v for k,v in raw_bc_count.items() if v > t}
        
        if not output_empty:
            return cells_bc, []
        else:
            # create a confident list of empty drops in high sensitivity mode
            logger.info("Creating emtpy droplets barocde list...")
            ept_bc = []
            ept_bc_max_count = min(cells_bc.values())
            ept_bc_max_count = min(ept_bc_max_count, empty_max_count)
            ept_bc_candidate = \
                [k for k,v in raw_bc_count.items() if v < ept_bc_max_count]
            for k in ept_bc_candidate:
                if min([edit_distance(k, x, max_ed = DEFAULT_EMPTY_DROP_MIN_ED) for x in cells_bc.keys()]) >= DEFAULT_EMPTY_DROP_MIN_ED:
                    ept_bc.append(k)

                # we don't need too much BC in this list
                if len(ept_bc) >  DEFAULT_EMPTY_DROP_NUM:
                    break
            return cells_bc, ept_bc

    else: 
        raise ValueError('Invalid value of count_t and/or exp_cells.')

def knee_plot(counts, threshold=None, out_fn = 'knee_plot.png'):
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
    plt.savefig(out_fn)


def read_batch_generator(fastq_fns, batch_size):
    """Generator of barches of reads from list of fastq files

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq = helper.fastq_parser(handle)
                read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            with open(fn, "r") as handle:
                fastq = helper.fastq_parser(handle)
                read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch



def main(argv=None):
    
    fastq_fns, out_fastq_fn, n_process, exp_cells ,min_phred_score, \
        full_bc_whitelist, out_raw_bc_fn, out_whitelist_fn, \
        high_sensitivity_mode, batch_size, out_emptydrop_fn, \
        emptydrop_max_count, overwrite, out_plot_fn, do_demultiplexing, \
        max_edit_distance, summary_fn,minimal_out  = parse_arg(argv)
    
    # Start running: Welcome logo
    if not minimal_out:
        print(textwrap.dedent(
            f'''\n\nWelcome to 
                {BLAZE_LOGO}
        '''))

    ######################
    ###### Getting putative barcodes
    ######################
    if not os.path.exists(out_raw_bc_fn) or overwrite:
        if os.path.exists(out_raw_bc_fn) and overwrite:
            logger.info(helper.warning_msg(
                f"The output putative barcodes table `{out_raw_bc_fn}` exist. It will be overwritten...",
                 printit = False))
        logger.info(f'Getting putative barcodes from {len(fastq_fns)} FASTQ files...')
        read_batchs = read_batch_generator(fastq_fns, batch_size=batch_size)
        
        # # temp---------------------------------------------
        # import cProfile
        #     # Create a cProfile object
        # profiler = cProfile.Profile()

        # # Start profiling
        # profiler.enable()
        # get_raw_bc_from_reads(next(read_batchs), min_q=min_phred_score)
        # # Stop profiling
        # profiler.disable()
        # profiler.dump_stats('profiling_results.cprof')
        # profiler.print_stats(sort='cumulative')

        # exit()
        # # ---------------------------------------------------------------------------

        rst_futures = helper.multiprocessing_submit(get_raw_bc_from_reads,
                                                read_batchs, n_process=n_process, 
                                                min_q=min_phred_score)

        raw_bc_count = Counter([])
        raw_bc_pass_count = Counter([])    
        rst_dfs = []
        for idx, f in enumerate(rst_futures):
            count_bc, count_pass, rst_df = f.result() #write rst_df out

            raw_bc_count += count_bc
            raw_bc_pass_count += count_pass
            if idx == 0:
                rst_df.to_csv(out_raw_bc_fn, index=False)
            else:
                rst_df.to_csv(out_raw_bc_fn, mode='a', index=False, header=False)
        
        helper.green_msg(f'Putative barcode table saved in {out_raw_bc_fn}')
        
        # ----------------------stats of the putative barcodes--------------------------
        qc_report(raw_bc_pass_count, min_phred_score=min_phred_score, 
                  out_fn=summary_fn, stdout= not minimal_out)

    
    else:
        logger.info(helper.warning_msg(textwrap.dedent(
            f"""
            Warning: `{out_raw_bc_fn}` exists. BLAZE would NOT re-generate the file and the existing file 
            will be directly used for downstream steps. If you believe it needs to be updated, please
            change the --output_prefix or remove/rename the existing file. 
            
            Note: there is no need to update this file if the input data remain the same and the previous
            run that generated this file finished successfully. It wouldn't change with other specified 
            arguments . However if you are running using a modified config.py file or the existing  `{out_raw_bc_fn}`
            was generated by a different version of BLAZE, updating the file is suggested.
            """
        ), printit = False))
            # read table
        dfs = pd.read_csv(out_raw_bc_fn, chunksize=1_000_000)
        
        # get bc count dict (filtered by minQ)
        
        raw_bc_count = Counter()
        for df in tqdm(dfs, desc = 'Counting high-quality putative BC'):
            raw_bc_count += Counter(df[
                df.putative_bc_min_q >=min_phred_score].putative_bc.value_counts().to_dict())
            
        
    ######################
    ###### Whitelisting
    ######################
    if not os.path.exists(out_whitelist_fn) or overwrite:
        if overwrite or not os.path.exists(out_emptydrop_fn):
            logger.info("Getting barcode whitelist and empty droplet barcode list...\n")

        if overwrite:
            logger.info(helper.warning_msg(
                f"Warning: `{out_whitelist_fn}` and `{out_emptydrop_fn}` will be overwritten if exist...",
                printit = False
            ))

        elif os.path.exists(out_emptydrop_fn):
            logger.info(helper.warning_msg(
                f"Warning: `{out_whitelist_fn}` doesn't exist, `{out_emptydrop_fn}` and `{out_plot_fn}` will be overwritten if exist...",
                printit = False
            ))

        try:
            bc_whitelist, ept_bc = get_bc_whitelist(raw_bc_count,
                                    full_bc_whitelist, 
                                    exp_cells=exp_cells,
                                    high_sensitivity_mode=high_sensitivity_mode,
                                    output_empty=True,
                                    empty_max_count=emptydrop_max_count,
                                    out_plot_fn = out_plot_fn)
            with open(out_whitelist_fn, 'w') as f:
                for k in bc_whitelist.keys():
                    f.write(k+'-1\n')

            with open(out_emptydrop_fn, 'w') as f:
                for k in ept_bc:
                    f.write(k+'-1\n')
                helper.green_msg(f'Empty droplet barcode list saved as `{out_emptydrop_fn}`.')    

        except Exception as e:
            logger.exception(e)
            helper.err_msg(
                "Error: Failed to get whitelist. Please check the input files and settings."
                )
        else:
            helper.green_msg(f'Whitelist saved as `{out_whitelist_fn}`!')

    elif os.path.exists(out_whitelist_fn) and not overwrite:
        logger.info(helper.warning_msg(
                f"Warning: `{out_whitelist_fn}` exist, the whitelisting step will be skipped."
                , printit = False))
        if not os.path.exists(out_emptydrop_fn):
            logger.info(helper.warning_msg(
                f"Warning: BLAZE will use existing `{out_whitelist_fn}` for the downstread steps and will not re-generate the {out_emptydrop_fn}."
                f"If the file is required, please remove/rename the existing `{out_whitelist_fn}` and rerun."
            , printit = False))
        if not os.path.exists(out_plot_fn):
            logger.info(helper.warning_msg(
                f"Warning: BLAZE will use existing `{out_whitelist_fn}` for the downstread steps and will not re-generate the {out_plot_fn}."
                f"If the file is required, please remove/rename the existing `{out_whitelist_fn}` and rerun."
            , printit = False))
    
    ######################
    ###### Demultiplexing
    ######################
    if do_demultiplexing and (not os.path.exists(out_fastq_fn) or overwrite):
        if overwrite:
            logger.info(helper.warning_msg(
                f"Warning:  {out_fastq_fn} will be overwritten if exist...",
                printit = False
            ))
        logger.info("Assigning reads to whitelist.\n")
        read_assignment.main_multi_thread(fastq_fns, 
                                        out_fastq_fn, 
                                        out_raw_bc_fn, 
                                        out_whitelist_fn,
                                        max_edit_distance,
                                        n_process,
                                        out_fastq_fn.endswith('.gz'), 
                                        batch_size)
    elif os.path.exists(out_fastq_fn) and \
        os.path.getmtime(out_fastq_fn) < os.path.getmtime(out_whitelist_fn):
        logger.info(helper.warning_msg(
            f"Warning: The `{out_fastq_fn}` exists and has NOT been updated. However,"
            f"the existing `{out_fastq_fn}` is older than the upstream output {out_whitelist_fn}."
            f"If it needs to be re-generated. Please remove/rename the existing `{out_fastq_fn}` and re-run BLAZE "
        , printit = False))
if __name__ == '__main__':
    main()
