# Generate new whitelist from baw_bc table

import argparse
import textwrap
import pandas as pd

from get_raw_bc import get_bc_whitelist
from config import *
import helper
import sys

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
                        help='Output the whitelist in Cellranger style')

    # Optional positional argument
    #parser.add_argument('opt_pos_arg', type=int, nargs='?',help)

    # Optional argument
    parser.add_argument('--minQ', type=int, default=15,
                        help= textwrap.dedent(
                            '''
                            <INT>: Minimum phred score for all bases in a raw BC. 
                            Reads whose raw BC contains one or more bases with 
                            Q<minQ is not counted in the "Raw BC rank plot".'''))

    parser.add_argument('--full-bc-whitelist', type=str, default=DEFAULT_GRB_WHITELIST,
                        help='''<path to file>: .txt file containing all the possible BCs. Users may provide
        their own whitelist.''')
    parser.add_argument('--out-bc-whitelist', type=str, default=DEFAULT_GRB_OUT_WHITELIST,
                        help='''<filename_prefix>: Output the whitelist identified from all the reads.''')
    parser.add_argument('--cr_style', type=bool, nargs='?',const=True, default=True,
                        help='Output the whitelist in Cellranger style')
    
    args = parser.parse_args()

    if not args.expect_cells and not args.count_threshold:
        helper.err_msg("Missing argument --expect-cells or --count-threshold.") 
        sys.exit(1)
    if args.expect_cells and args.count_threshold:
        helper.warning_msg(textwrap.dedent(
                f'''
                You have specified both '--expect-cells' and '--count-threshold'. \
'--expect-cells' will be ignored.                
                '''))
    return args



def main(args):
    # read table
    df = pd.read_csv(args.raw_bc_csv)
    # get bc count dict (filtered by minQ)
    raw_bc_count = df[df.raw_bc_min_q >=args.minQ].raw_bc.value_counts().to_dict()
    
    bc_whitelist = get_bc_whitelist(raw_bc_count,
                                    args.full_bc_whitelist, 
                                    args.expect_cells,
                                    args.count_threshold)

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