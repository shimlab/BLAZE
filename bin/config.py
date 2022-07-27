# This file stores the parameter used in this repo

import os
import numpy as np

####################################################
############# polyT and adaptor finding#############
####################################################
## adaptor finding

ADPT_SEQ='CTTCCGATCT' #searched adaptor sequence
ADPT_WIN=200 #search adaptor in subsequence from both end of the reads with this size
ADPT_MIN_MATCH_PROP=0.8 #minimum proportion of match required when searching

## poly T searching
PLY_T_LEN=4 #length of searched poly T
PLY_T_WIN=200 #search poly T in subsequence from both end of the reads with this size
PLY_T_MIN_MATCH_PROP=1#minimum proportion of match required when searching
PLY_T_NT_AFT_ADPT=(20,50)#a poly T should locate within this range downstream an adaptor



####################################################
############    DEFAULT in get_raw_bc          #####
####################################################
# input
DEFAULT_GRB_MIN_SCORE=15
DEFAULT_GRB_KIT='v3'

DEFAULT_GRB_WHITELIST_V3=\
        os.path.join(os.path.dirname(__file__), '../10X_bc/3M-february-2018.zip')
DEFAULT_GRB_WHITELIST_V2=\
        os.path.join(os.path.dirname(__file__), '../10X_bc/737K-august-2016.txt')

#output
DEFAULT_GRB_OUT_RAW_BC='putative_bc'
DEFAULT_GRB_OUT_WHITELIST = 'whitelist'

# defualt count threshold
def default_count_threshold_calculation(count_array, exp_cells):
    top_count = np.sort(count_array)[::-1][:exp_cells]
    return np.quantile(top_count, 0.95)/20
