# Generate new whitelist from baw_bc table

import argparse
import textwrap
import pandas as pd
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm

from blaze.blaze import get_bc_whitelist
from blaze.config import *
import blaze.helper as helper


def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        This script can be used to generate a new whitelist from the putative bc table
        output from 'blaze.py'. Users may specify different argment used in 
        'blaze.py' to obtain a different whitelist.
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Required positional argument
    parser.add_argument('putative_bc_csv', type=str,
                        help='Filename of the putative_bc csv file output from blaze.py')

    # required name argment
    requiredNamed = parser.add_argument_group('Either one of these argument is required')
    requiredNamed.add_argument('--expect-cells',type=int, help='<INT>:  Expected number of cells.')
    requiredNamed.add_argument('--count-threshold', type=int,
                        help='Output the whitelist in Cellranger style')

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
                            <INT>: Minimum phred score for all bases in a putative BC. 
                            Reads whose putative BC contains one or more bases with 
                            Q<minQ is not counted in the "Putative BC rank plot".'''))
    parser.add_argument('--full-bc-whitelist', type=str, default=None,
                        help='''<path to file>: .txt file containing all the possible BCs. Users may provide
        their own whitelist. No need to specify this if users want to use the 10X whilelist. The correct version
        of 10X whilelist will be determined based on 10X kit version''')
    parser.add_argument('--out-bc-whitelist', type=str, default=DEFAULT_GRB_OUT_WHITELIST,
                        help='''<filename_prefix>: Output the whitelist identified from all the reads.''')
    parser.add_argument('--cr-style', type=bool, nargs='?',const=True, default=True,
                        help='Output the whitelist in Cellranger style')
    parser.add_argument('--chunk-size', type=int, default=1_000_000,
                        help='Chunk size when reading the input file. Please use'
                        'smaller number if memory is not sufficient.')
    parser.add_argument('--high-sensitivity-mode', action='store_true',
                        help='''Turn on the sensitivity mode, which increases the sensitivity of barcode
                    detections but potentially increase the number false/uninformative BC in
                    the whitelist.''')
    parser.add_argument('--emptydrop', action='store_true',
                        help='''Output list of BCs corresponding to empty droplets (filename: {DEFAULT_EMPTY_DROP_FN}), 
                    which could be used to estimate ambiant RNA expressionprofile.''')
    parser.add_argument('--emptydrop-max-count', type=float, default=np.inf,
                        help=textwrap.dedent(
                            '''
                                Only select barcodes supported by a maximum number of high-confidence
                                putative barcode count. (Default: Inf, i.e. no maximum number is set
                                and any barcodes with ED>= {DEFAULT_EMPTY_DROP_MIN_ED} to the barcodes
                                in whitelist can be selected as empty droplets)'''))
    args = parser.parse_args()

    if not args.expect_cells and not args.count_threshold:
        helper.err_msg("Missing argument --expect-cells or --count-threshold.") 
        sys.exit(1)
    if (args.expect_cells or args.high_sensitivity_mode) and args.count_threshold:
        helper.warning_msg(textwrap.dedent(
                f'''
                Warning: You have specified'--count-threshold'. Options
                "--high_sensitivity_mode" and "--expect-cells" would be ignored if
                specified.        
                '''))
        args.high_sensitivity_mode = False
    
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
    helper.check_exist([args.full_bc_whitelist, args.putative_bc_csv])
    return args



def main(args):
    # read table
    dfs = pd.read_csv(args.putative_bc_csv, chunksize=args.chunk_size)
    
    # get bc count dict (filtered by minQ)
    
    raw_bc_count = Counter()
    for df in tqdm(dfs, desc = 'Counting high-quality putative BC'):
        raw_bc_count += Counter(df[
            df.putative_bc_min_q >=args.minQ].putative_bc.value_counts().to_dict())

    if args.high_sensitivity_mode:
        print('Preparing whitelist...(high-sensitivity-mode)')
    else: 
        print('Preparing whitelist...')
    bc_whitelist, ept_bc = get_bc_whitelist(raw_bc_count,
                                    args.full_bc_whitelist, 
                                    args.expect_cells,
                                    args.count_threshold,
                                    args.high_sensitivity_mode,
                                    args.emptydrop,
                                    args.emptydrop_max_count)
    if args.cr_style:
        with open(args.out_bc_whitelist+'.csv', 'w') as f:
            for k in bc_whitelist.keys():
                f.write(k+'-1\n')
        if ept_bc:
            with open(DEFAULT_EMPTY_DROP_FN, 'w') as f:
                for k in ept_bc:
                    f.write(k+'-1\n')
    else:
        with open(args.out_bc_whitelist+'.csv', 'w') as f:
            for k in bc_whitelist.keys():
                f.write(k+'\n')
        if ept_bc:
            with open(DEFAULT_EMPTY_DROP_FN, 'w') as f:
                for k in ept_bc:
                    f.write(k+'\n')
    helper.green_msg(f'Whitelist saved as {args.out_bc_whitelist}.csv!')
if __name__ == '__main__':
    args = parse_arg()
    #print(args)
    main(args)
