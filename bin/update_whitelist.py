# Generate new whitelist from baw_bc table

import argparse
import textwrap
import pandas as pd
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm
import numpy as np

from get_raw_bc import get_bc_whitelist
from config import *
import helper


def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        This script can be used to generate a new whitelist from the raw_bc table
        output from 'get_raw_bc.py'. Users may specify different argment used in 
        'get_raw_bc.py' to obtain a different whitelist.
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Required positional argument
    parser.add_argument('raw_bc_csv', type=str,
                        help='Filename of the raw_bc_csv')

    # required name argment
    requiredNamed = parser.add_argument_group('Either one of these argument is required')
    requiredNamed.add_argument('--expect-cells',type=int, help='<INT>:  Expected number of cells.')
    requiredNamed.add_argument('--count-threshold', type=int,
                        help='Use a user-specified count threshold for BCs in '
                        'the output whitelist. BCs with high-confidence raw '
                        'BC count larger than the threshold will be included in'
                        ' the output whitelist. Note that specifying a count '
                        'threshold means any specified value for --expect-cells' 
                        ' and --high-sensitivity-mode will be ignored')
    requiredNamed.add_argument('--high-sensitivity-mode', type=bool, nargs='?',
                    const=True, default=False,
                    help='Use a high-sensitivity way of choosing whitelist'
                    '(using the knee point of cumulative count curve). '
                    'With --high-sensitivity-mode True, the output would include '
                    'more cell-associated BC with risks of also including false '
                    'positive BCs （i.e. BCs in empty droplets and/or BCs do '
                    'not exist）')

    # Optional positional argument
    #parser.add_argument('opt_pos_arg', type=int, nargs='?',help)

    # Optional argument
    parser.add_argument('--kit-version', type=str, default='v3',
                        help= textwrap.dedent(
                            '''
                            Choose from v2 and v3 (for 10X Single Cell 3ʹ gene expression v2 or v3). 
                            '''))
    parser.add_argument('--minQ', type=int, default=15,
                        help= textwrap.dedent(
                            '''
                            <INT>: Minimum phred score for all bases in a raw BC. 
                            Reads whose raw BC contains one or more bases with 
                            Q<minQ is not counted in the "Raw BC rank plot".'''))

    parser.add_argument('--full-bc-whitelist', type=str, default=None,
                        help='''<path to file>: .txt file containing all the possible BCs. Users may provide
        their own whitelist. No need to specify this if users want to use the 10X whilelist. The correct version
        of 10X whilelist will be determined based on 10X kit version''')
    parser.add_argument('--out-bc-whitelist', type=str, default=DEFAULT_GRB_OUT_WHITELIST,
                        help='''<filename_prefix>: Output the whitelist identified from all the reads.''')
    parser.add_argument('--cr-style', type=bool, nargs='?',const=True, default=True,
                        help='Output the whitelist in Cellranger style')
    parser.add_argument('--chunk-size', type=int, default=1_000_000,
                        help='Chunksize when reading the input file. Please use'
                        'smaller number if memory is not sufficient.')

    
    args = parser.parse_args()

    if not args.expect_cells and not args.count_threshold \
                             and not args.high_sensitivity_mode:
        helper.err_msg("Missing argument --expect-cells or --count-threshold"
                        "or --high-sensitivity-mode.") 
        sys.exit(1)
    if args.expect_cells and args.count_threshold:
        helper.warning_msg(textwrap.dedent(
                f'''
                Warning: You have specified both '--expect-cells' and '--count-threshold'. \
'--expect-cells' will be ignored.                
                '''))
        args.expect_cells = None
    
    if args.high_sensitivity_mode and args.count_threshold:
        helper.warning_msg(textwrap.dedent(
                f'''
                Warning: You have specified both '--high-sensitivity-mode' and '--count-threshold'. \
'--high-sensitivity-mode' will be ignored.                
                '''))
        args.high_sensitivity_mode = False
    
    if args.high_sensitivity_mode and args.expect_cells and not args.count_threshold:
        helper.warning_msg(textwrap.dedent(
                f'''
                Warning: You have specified both '--high-sensitivity-mode' and '--expect-cells'. \
the output will be in high sensitivity mode and '--expect-cells' will be ignored.                
                '''))
        args.expect_cells = None
    
    
    args.kit_version = args.kit_version.lower()
    if args.kit_version not in ['v2', 'v3']:
        helper.err_msg("Error: Invalid value of --kit-version, please choose from v3 or v2") 
        sys.exit()

    if args.full_bc_whitelist:
        helper.warning_msg(textwrap.dedent(
                f'You are using {os.path.basename(args.full_bc_whitelist)} as the full barcode'\
                'whitelist. Note that the barcodes not listed in the file will never be found.'))
    else:
        if args.kit_version == 'v3':
            args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_V3
        elif args.kit_version == 'v2':
            args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_V2

    # check file 
    helper.check_exist([args.full_bc_whitelist, args.raw_bc_csv])
    return args



def main(args):
    # read table
    dfs = pd.read_csv(args.raw_bc_csv, chunksize=args.chunk_size)
    
    # get bc count dict (filtered by minQ)
    
    raw_bc_count = Counter()
    for df in tqdm(dfs, desc = 'Counting high-confidence raw BC'):
        raw_bc_count += Counter(df[
            df.raw_bc_min_q >=args.minQ].raw_bc.value_counts().to_dict())

    print(f'{sum(raw_bc_count.values())} raw BCs have'
    f' minimum base quality >= {args.minQ}')
    print('Preparing whitelist...')
    bc_whitelist = get_bc_whitelist(raw_bc_count,
                                    args.full_bc_whitelist, 
                                    args.expect_cells,
                                    args.count_threshold,
                                    args.high_sensitivity_mode)

    if args.cr_style:
        with open(args.out_bc_whitelist+'.csv', 'w') as f:
            for k in bc_whitelist.keys():
                f.write(k+'-1\n')
    else:
        with open(args.out_bc_whitelist+'.csv', 'w') as f:
            for k in bc_whitelist.keys():
                f.write(k+'\n')
    helper.green_msg(f'Whitelist saved as {args.out_bc_whitelist}.csv!')
if __name__ == '__main__':
    args = parse_arg()
    #print(args)
    main(args)