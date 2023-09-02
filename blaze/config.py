# This file stores the parameter used in this repo

import os
import numpy as np

## Output prefix
DEFAULT_PREFIX = ''

####################################################
############# polyT and adaptor finding#############
####################################################
## adaptor finding
ADPT_SEQ='CTTCCGATCT' #searched adaptor sequence
ADPT_WIN=200 #search adaptor in subsequence from both end of the reads with this size
ADPT_MAC_MATCH_ED=2 #minimum proportion of match required when searching

## poly T searching
PLY_T_LEN=4 #length of searched poly T
PLY_T_WIN=200 #search poly T in subsequence from both end of the reads with this size
PLY_T_MIN_MATCH_PROP=1#minimum proportion of match required when searching
PLY_T_NT_AFT_ADPT=(20,50)#a poly T should locate within this range downstream an adaptor



####################################################
#######    DEFAULT in getting putative bc     ######
####################################################
# input
DEFAULT_GRB_MIN_SCORE=15
DEFAULT_GRB_KIT='v3'
DEFAULT_UMI_SIZE= 12 if DEFAULT_GRB_KIT=='v3' else 10

# The 10X barcode whitelists has been packed in the package
DEFAULT_GRB_WHITELIST_V3=os.path.join(os.path.dirname(__file__), '10X_bc', '3M-february-2018.zip')
DEFAULT_GRB_WHITELIST_V2=os.path.join(os.path.dirname(__file__), '10X_bc', '737K-august-2016.txt')

#output
DEFAULT_GRB_OUT_RAW_BC='putative_bc.csv'
DEFAULT_GRB_OUT_WHITELIST = 'whitelist.csv'
DEFAULT_GRB_OUT_FASTQ = "matched_reads.fastq.gz"
DEFAULT_GRB_FLANKING_SIZE = 5

####################################################
#####    DEFAULT in generating  whitelist     ######
####################################################
# quantile based threshold
def default_count_threshold_calculation(count_array, exp_cells):
    top_count = np.sort(count_array)[::-1][:exp_cells]
    return np.quantile(top_count, 0.95)/20

def high_sensitivity_threshold_calculation(count_array, exp_cells):
    top_count = np.sort(count_array)[::-1][:exp_cells]
    return np.quantile(top_count, 0.95)/200

# list for empty drops (output in high-sensitivity mode)
DEFAULT_EMPTY_DROP_FN = 'emtpy_bc_list.csv'
DEFAULT_KNEE_PLOT_FN = 'knee_plot.png'
DEFAULT_BC_STAT_FN = "summary.txt"
DEFAULT_EMPTY_DROP_MIN_ED = 5 # minimum edit distance from emtpy drop BC to selected BC
DEFAULT_EMPTY_DROP_NUM = 2000 # number of BC in the output
    

####################################################
#####    DEFAULT in Demultiplexing            ######
####################################################

DEFAULT_ASSIGNMENT_ED = 2 
# Make sure this is smaller than DEFAULT_GRB_FLANKING_SIZE
assert DEFAULT_GRB_FLANKING_SIZE >= DEFAULT_ASSIGNMENT_ED
DEFAULT_ED_FLANKING = DEFAULT_ASSIGNMENT_ED



BLAZE_LOGO = \
"""
BBBBBBBBBBBBBBB  LLLLL       AAAAAA     ZZZZZZZZZZZZZZEEEEEEEEEEEEEEE
BBBBB&&&&&&&BBBB LLLLL      AAAAAAAA    ZZZZZZZZZZZZZZEEEEEEEEEEEEEEE
BBBBB^^^^^^!BBBB LLLLL     AAAAAAAAAA.         ZZZZZ. EEEE
BBBBB       BBBB LLLLL    AAAA    AAAA.       ZZZZZ.  EEEEEEEEEEEEEEE
BBBBBBBBBBBBBBB  LLLLL   AAAAAAAAAAAAAA      ZZZZZ.   EEEEEEEEEEEEEEE
BBBBBBBBBBBBBBB  LLLLL  AAAAAAAAAAAAAAAA.   ZZZZZ.    EEEE
BBBBB       BBBB LLLLL AAAAAA       AAAAA ZZZZZZZZZZZZEEEEEEEEEEEEEEE
BBBBB       BBBB LLLLLAAAAAA         AAAAAZZZZZZZZZZZZEEEEEEEEEEEEEEE
BBBBBBBBBBBBBBBB LLLLLLLLLLLLLLLLLLLL . ^PPPPPPY.   ^PPPPPPY.    7PPP
BBBBBBBBBBBBBBB  LLLLLLLLLLLLLLLLLLLL...!BBBBBBP:...!BBBBBBP:.::.?BBB
"""

