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
from collections import defaultdict, Counter
from tqdm import tqdm
import textwrap
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import zipfile
import io
import logging
import gzip
from fast_edit_distance import edit_distance
import logging


from blaze.parser import parse_arg
import blaze.helper as helper
from blaze.config import *
import blaze.polyT_adaptor_finder as polyT_adaptor_finder
import blaze.read_assignment as read_assignment

# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Parse fastq -> polyT_adaptor_finder.Read class
def get_raw_bc_from_reads(reads, min_q=0, kit=None, **kwargs):
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
                    phred_score=r.q_letter, kit=kit, **kwargs)    
        

        read.get_strand_and_raw_bc()
        read_ids.append(read.id)
        putative_bcs.append(read.raw_bc)
        putative_bc_min_qs.append(read.raw_bc_min_q)
        umis.append(read.putative_UMI)
        trim_idxs.append(read.polyT_trimming_idx)
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
        'polyT_end': trim_idxs,
        'pre_bc_flanking': pre_bc_flankings,
        'post_umi_flanking': post_umi_flankings
        }
        )
    return Counter(raw_bc_pass), rst_df

def add_summary(content, args, write_mode='w'):
    with open(args.summary_fn, write_mode) as f:
        f.write(content)
    if not args.minimal_stdout:
        print(content)

def bc_search_qc_report(pass_count, args):
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
            {pass_count[0]:,} in which all bases in the putative BC have Q>={args.minQ}
        Failed Reads: 
            no polyT and adapter positions found: 
                {pass_count[1]:,} ({pass_count[1]/total_read*100:.2f}% of all reads)
            polyT and adapter positions found in both end (fail to determine strand): 
                {pass_count[2]:,} ({pass_count[2]/total_read*100:.2f}% of all reads)
            multiple polyT and adapter found in one end
                {pass_count[10]:,} ({pass_count[10]/total_read*100:.2f}% of all reads)
        -------------------------------------------------------------------------------\n
        ''')
    return print_message

def get_bc_whitelist(raw_bc_count, full_bc_whitelist=None, exp_cells=None, 
                    count_t=None,force_cell_n=None, high_sensitivity_mode=False, 
                    output_empty = True, empty_max_count = np.inf, 
                    out_plot_fn = DEFAULT_KNEE_PLOT_FN, args=None):
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
            2. At least {out_plot_fn} away from all selected cell-associated BCs in whitelist. and
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

    # Check args:
    if args:
        full_bc_whitelist = args.full_bc_whitelist
        exp_cells = args.expect_cells
        count_t = args.count_threshold
        high_sensitivity_mode = args.high_sensitivity_mode
        out_plot_fn = args.out_plot_fn
        output_empty=args.out_emptydrop_fn
        empty_max_count = args.empty_max_count
        force_cell_n = args.force_cells
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
    if force_cell_n:
        if force_cell_n > len(raw_bc_count):
            logger.warning(helper.warning_msg(
                f"force_cells ({force_cell_n}) is larger than the number of unique barcodes found in the data ({len(raw_bc_count)})."))
            count_t = 0
        # count threshold is the minimum count of the top N cells
        else:
            count_t = sorted(list(raw_bc_count.values()))[-force_cell_n]
    if count_t is not None:
        knee_plot(list(raw_bc_count.values()), count_t, out_plot_fn)
        cells_bc = {k:v for k,v in raw_bc_count.items() if v >= count_t}
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

def print_logo(args):
    if not args.minimal_stdout:
        print(textwrap.dedent(
            f'''\n\nWelcome to 
                {BLAZE_LOGO}
        '''))

def main(args_string=None):
    # Start running: Welcome logo
    args = parse_arg(args_string)
    print_logo(args)

    # TMP: print all the arguments in args
    #for arg in vars(args):
    #    print(f"{arg}: {getattr(args, arg)}")
    
    ######################
    ###### Getting putative barcodes
    ######################
    if args.do_bc_search:
        logger.info(f'Getting putative barcodes from {len(args.fastq_fns)} FASTQ files...')
        read_batchs = read_batch_generator(args.fastq_fns, batch_size=args.batch_size)

            
        rst_futures = helper.multiprocessing_submit(get_raw_bc_from_reads,
                                                read_batchs, n_process=args.threads, 
                                                min_q=args.minQ, kit=args.kit_version,
                                                umi_len=args.umi_len)
    

        raw_bc_pass_count = defaultdict(int)

        for idx, f in enumerate(rst_futures):
            count_pass, rst_df = f.result() #write rst_df out
            for k,v in count_pass.items():
                raw_bc_pass_count[k] += v
            if idx == 0:
                rst_df.to_csv(args.out_raw_bc_fn, index=False)
            else:
                rst_df.to_csv(args.out_raw_bc_fn, mode='a', index=False, header=False)
        
        helper.green_msg(f'Putative barcode table saved in {args.out_raw_bc_fn}', printit=True)
        
        # ----------------------stats of the putative barcodes--------------------------

        add_summary(bc_search_qc_report(raw_bc_pass_count, args), 
                    args=args, write_mode='w')

    ######################
    ###### Whitelisting
    ######################
    if args.do_whitelisting:
        # get bc count dict (filtered by minQ)
        dfs = pd.read_csv(args.out_raw_bc_fn, chunksize=1_000_000)
        raw_bc_count = Counter()
        for df in tqdm(dfs, desc = 'Counting high-quality putative BC', unit='M reads'):
            raw_bc_count += Counter(df[
                df.putative_bc_min_q >=args.minQ].putative_bc.value_counts().to_dict())

        try:
            bc_whitelist, ept_bc = get_bc_whitelist(raw_bc_count, args=args)

            # output
            with open(args.out_whitelist_fn, 'w') as f:
                for k in bc_whitelist.keys():
                    f.write(k+'\n')
            with open(args.out_emptydrop_fn, 'w') as f:
                for k in ept_bc:
                    f.write(k+'\n')
                helper.green_msg(f'Empty droplet barcode list saved as `{args.out_emptydrop_fn}`.', printit=True)    
            # write to summary
            add_summary(f'\nIdentified # of cells: {len(bc_whitelist)}\n', args, write_mode='a')

        except Exception as e:
            logger.exception(e)
            helper.err_msg(
                "Error: Failed to get whitelist. Please check the input files and settings." , printit=True
                )
        else:
            helper.green_msg(f'Whitelist saved as `{args.out_whitelist_fn}`!', printit=True)


    ######################
    ###### Demultiplexing
    ######################
    if args.do_demultiplexing:
        logger.info("Assigning reads to whitelist.\n")
        # write to fastq
        demul_count_tot, count_tot = read_assignment.assign_read(args=args) 
        # write to summary
        add_summary(f'\nTotal reads: {count_tot}'
                    f'\nTotal reads in cells: {demul_count_tot} ({demul_count_tot/count_tot*100:.2f}%)', args, write_mode='a')

    
if __name__ == '__main__':
    main()
