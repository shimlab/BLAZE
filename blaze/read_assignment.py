from fast_edit_distance import edit_distance, sub_edit_distance
import multiprocessing as mp
import tempfile
import pandas as pd
import numpy as np
import gzip
import os
from tqdm import tqdm
import logging
import sys

import blaze.helper as helper
from blaze.config import *


# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def match_bc(bc, whitelist, max_ed):
    """
    Args:
        bc (str): putative bc to be assigned
        whitelist (set): set of bc to assign to
        max_ed (int): maximum edit distance allowed for a assignment
    Returns:
        assigned barcode <str>:  if a unambiguous barcode was found in the whitelist
        '': if no barcode was found in the whitelist
        'ambiguous': if multiple barcode was found in the whitelist
    """
    if not bc:
        return ''
    if bc in whitelist:
        return bc
    best_ed = max_ed
    bc_hit = ''
    for i in whitelist:
        ed = edit_distance(bc, i, best_ed) 
        if ed < best_ed:
            best_ed = ed
            bc_hit = i
        elif ed == best_ed:
            if not bc_hit:
                bc_hit = i
            else: 
                bc_hit = 'ambiguous'
                best_ed -= 1
                if best_ed < 0:
                    return ''
    return bc_hit

def test_match_bc_row(row, whitelist, max_ed):
    rst = match_bc_row(row, whitelist, max_ed)
    if len(rst) != 3:
        print(len(rst))
    return rst

def match_bc_row(row, whitelist, max_ed):
    """
    Args:
        row (row in pandas.DataFrame)
        whitelist (set): list of bc to assign to
        max_ed (int): maximum edit distance allowed for a assignment
    Returns:
        1. assigned barcode <str>:  if a unambiguous barcode was found in the whitelist
            '': if no barcode was found in the whitelist
            'ambiguous': if multiple barcode was found in the whitelist
        2. adjusted putative_umi
        3. strand '+' for read with positive strand (start with BC,umi,polyT...), 
                  '-' for read with negative strand

    """
    if not row.umi_end:
        strand = ''
    elif row.umi_end > 0:
        strand = '+'
    else: 
        strand = '-'

    if not row.putative_bc or row.putative_bc in whitelist:
        return pd.Series([row.putative_bc, row.putative_umi, strand])
    else:
        bc = row.putative_bc
    
    best_ed = max_ed
    bc_hit = ''
    bc_hit_end_idx = -1
    # extending the putative barcode from both sides for potential indels
    bc = row.pre_bc_flanking[-DEFAULT_ED_FLANKING:] + bc + row.putative_umi[:DEFAULT_ED_FLANKING]
    for i in whitelist:
        ed, end_idx = sub_edit_distance(i, bc, best_ed) 
        if ed < best_ed:
            best_ed = ed
            bc_hit = i
            bc_hit_end_idx = end_idx
        elif ed == best_ed:
            if not bc_hit:
                bc_hit = i
                bc_hit_end_idx = end_idx
            else: 
                bc_hit = 'ambiguous'
                best_ed -= 1
                if best_ed < 0:
                    return pd.Series(['', row.putative_umi, strand])
    
    if bc_hit == 'ambiguous' or bc_hit == '':
        return pd.Series(['', row.putative_umi, strand])
    else:
        pass

    # adjust the umi start position
    umi_adj = bc_hit_end_idx - (len(bc) - 1 -DEFAULT_ED_FLANKING )
    out_umi = row.putative_umi
    if umi_adj > 0:
        out_umi = row.putative_umi[umi_adj:] + row.post_umi_flanking[:umi_adj]
    elif umi_adj < 0:
        out_umi =  row.putative_bc[umi_adj:] + row.putative_umi[:umi_adj]
    return pd.Series([bc_hit, out_umi, strand])
            
# Function to modify the header
def fastq_modify_header(record, barcode, UMI):
    record.id = f"{barcode}_{UMI}#{record.id}"
    record.description = ""
    return record

def fastq_trim_seq(record, start = 0, end = None, return_str = False):
    """trim reads using 0-based index
    note end can be either positive or negative, following python idexing rules
    """
    if not end:
        qscore = \
            record.letter_annotations['phred_quality'][start:]
        seq = record.seq[start:]
        if return_str:
            return seq, qscore
        else:
            record.letter_annotations = {}
            record.seq = seq
            record.letter_annotations['phred_quality']=qscore
            return record
    else:
        qscore = \
            record.letter_annotations['phred_quality'][start:end]
        seq = record.seq[start:end]
        if return_str:
            return seq, qscore
        else:
            record.letter_annotations = {}
            record.seq = seq
            record.letter_annotations['phred_quality']=qscore
            return record


def batch_barcode_to_fastq(read_batches_with_idx, assignment_df ,gz = True):
    """Take a read batch, write fastq/fastq.gz to a tmp file
    """
    read_batch, start_df_idx = read_batches_with_idx
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_file_path = temp_file.name
    output_handle = gzip.open(temp_file_path, 'wt') if gz else open(temp_file_path, 'w')
    out_buffer = ""
    read_idx = 0
    for r in read_batch:
        row = assignment_df.iloc[start_df_idx+read_idx]#the row in putative bc table
        read_idx += 1 
        try:
            assert row.read_id == r.id            
        except AssertionError:
            helper.err_msg("Different order in putative bc file and input fastq!")
            sys.exit()

        if not row.BC_corrected:
            continue

        if row.umi_end < 0:
            seq = r.seq[:int(row.umi_end)]
            qscore = r.qscore[:int(row.umi_end)]
        else:
            seq = r.seq[int(row.umi_end):]
            qscore = r.qscore[int(row.umi_end):]
        
        out_buffer += f"@{row.BC_corrected}_{row.putative_umi}#{row.read_id}_{row.strand}\n"
        out_buffer += str(seq) + '\n'
        out_buffer += '+\n'
        out_buffer += qscore + '\n'
    
    output_handle.write(out_buffer)
    output_handle.close()
    return temp_file.name


def assign_barcodes(putative_bc_csv, whitelsit_csv, n_process, max_ed):
    """Read the putative barcodes as pd.dataframe,
    Assign all putative barcode to whiteliested barcode, 

    Args:
        putative_bc_csv (str): filename
        whitelsit_csv (str): filname
        n_process (int): number of process 

    Returns:
        pd.dataframe: Same as the input csv, with an additional column called BC_corrected
    """
    
    # read putative barcode
    df = pd.read_csv(putative_bc_csv)
    df = df.fillna('') # replace nan with empty string
    # df['putative_umi'] = df['putative_umi'].fillna('')
    # df['pre_bc_flanking'] = df['pre_bc_flanking'].fillna('')
    # df['post_bc_flanking'] = df['post_bc_flanking'].fillna('')   

    # read whitelist
    whitelist = [] 
    with open(whitelsit_csv, 'r') as f:
        for line in f:
            whitelist.append(line.split('-')[0])
    whitelist = set(whitelist)
    
    new_cols = pd.DataFrame(columns = ['BC_corrected','putative_umi', 'strand'])
    new_cols = helper.df_multiproceccing_apply(df, 
                                        match_bc_row,
                                        n_process = n_process,
                                        max_ed = max_ed,
                                        whitelist=whitelist
                                        )

    df[['BC_corrected','putative_umi', 'strand']] = new_cols


    logger.info(helper.green_msg(f"Demultiplexing finshied: ", printit = False))
    logger.info(helper.green_msg(f"Successfully demultiplexed reads / Total reads: {sum(df.BC_corrected!='')} / {len(df.BC_corrected)}. ", printit = False))
    return df

def main_multi_thread(fastq_fns, fastq_out, putative_bc_csv, 
                      whitelsit_csv, max_ed, n_process, gz, batchsize):
    """Main function: Demultiplex fastq files using putative barcode csv and whitelist csv
    """
    class read_fastq:
        """This class is for mimic the SeqIO fastq record. The SeqIO is not directly used because it's slow.
        """
        def __init__(self,title, sequence, qscore, quality_map = False):
            self.id = title.split()[0].strip("@")
            self.seq = sequence
            self.qscore = qscore

    def read_batch_generator_with_idx(fastq_fns, batch_size):
        """Generator of barches of reads from list of fastq files with the idx of the first read
        in each batch

        Args:
            fastq_fns (list): fastq filenames
            batch_size (int, optional):  Defaults to 100.
        """
        read_idx = 0
        for fn in fastq_fns:
            if str(fn).endswith('.gz'):
                with gzip.open(fn, "rt") as handle:
                    fastq =\
                        (read_fastq(title, sequence, qscore) for title, sequence, qscore in helper.fastq_parser(handle))

                    batch_iter = helper.batch_iterator(fastq, batch_size=batch_size)
                    
                    for batch in batch_iter:
                        batch_len = len(batch)
                        yield batch, read_idx
                        read_idx += batch_len
            else:
                with open(fn) as handle:
                    fastq =\
                        (read_fastq(title, sequence, qscore) for title, sequence, qscore in helper.fastq_parser(handle))
                    read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
                    for batch in read_batch:
                        batch_len = len(batch)
                        yield batch, read_idx
                        read_idx += batch_len
    
    assignment_df = assign_barcodes(putative_bc_csv, whitelsit_csv, n_process, max_ed)
    r_batches_with_idx = \
        read_batch_generator_with_idx(fastq_fns, batchsize)

    logger.info("Reads assignment completed.")
    logger.info(f"Writing to tmp fastq files...")
    rst_futures = helper.multiprocessing_submit(batch_barcode_to_fastq, 
                           r_batches_with_idx, 
                           n_process=n_process*4,
                           schduler = "thread",
                           pbar_func=lambda x: len(x[0]),
                            assignment_df = assignment_df,
                            gz = gz)
    
    tmp_files = []
    for f in  rst_futures:
        tmp_fn = f.result()
        tmp_files.append(f.result())
    logger.info(f"Concatenating tmp fastq files to {fastq_out}")
    helper.concatenate_files(fastq_out, *tmp_files)
    logger.info(helper.green_msg(f"Demultiplexed read saved in {fastq_out}!", printit = False))

if __name__ == '__main__':
    #main_multi_thread()
    pass