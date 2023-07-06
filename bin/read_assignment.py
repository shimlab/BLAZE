from fast_edit_distance import edit_distance, sub_edit_distance
import multiprocessing as mp
import tempfile
import pandas as pd
import numpy as np
import swifter
import gzip
import Bio
from Bio import SeqIO
import os
from tqdm import tqdm
import logging
import sys

import helper
from config import *


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

    """
    bc = row.putative_bc
    if not bc:
        return pd.Series(['', row.putative_umi])
    if bc in whitelist:
        return pd.Series([bc, row.putative_umi])
    best_ed = max_ed
    bc_hit = ''
    bc = row.pre_bc_flanking[-DEFAULT_ED_FLANKING:] + bc + row.putative_umi[:DEFAULT_ED_FLANKING]
    for i in whitelist:
        ed, end_idx = sub_edit_distance(i, bc, best_ed) 
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
                    return pd.Series(['', row.putative_umi])
    if bc_hit == 'ambiguous' or not bc_hit:
        return pd.Series(['', row.putative_umi])
    
    umi_adj = end_idx - (len(bc) - 1 -DEFAULT_ED_FLANKING )
    out_umi = row.putative_umi
    if umi_adj > 0:
        out_umi = row.putative_umi[umi_adj:] + row.post_umi_flanking[:umi_adj]
    elif umi_adj < 0:
        out_umi =  row.putative_bc[umi_adj:] + row.putative_umi[:umi_adj]
    

    return pd.Series([bc_hit, out_umi])
            
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

        if not row.putative_bc:
            continue

        if row.umi_end < 0:
            seq = r.seq[:int(row.umi_end)]
            qscore = r.qscore[:int(row.umi_end)]
        else:
            seq = r.seq[int(row.umi_end):]
            qscore = r.qscore[int(row.umi_end):]
        
        out_buffer += f"@{row.BC_corrected}_{row.putative_umi}#{row.read_id}\n"
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

    # df['BC_corrected'] =\
    #   df['putative_bc'].\
    #     swifter.allow_dask_on_strings(enable=True).\
    #         set_dask_scheduler('processes').\
    #             set_npartitions(n_process).\
    #                 apply(match_bc, whitelist=whitelist, max_ed=max_ed)
    

    df[['BC_corrected','putative_umi']] =\
      df.swifter.allow_dask_on_strings(enable=True).\
            set_dask_scheduler('processes').\
                set_npartitions(n_process).\
                    apply(match_bc_row, axis=1, whitelist=whitelist, max_ed=max_ed)
    

    #df = df.set_index('read_id')# need UMI
    print(df.head())
    return df

def main_multi_thread(fastq_fns, fastq_out, putative_bc_csv, whitelsit_csv, max_ed, n_process, gz, batchsize):
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
                         (read_fastq(title, sequence, qscore) for title, sequence, qscore in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle))

                    batch_iter = helper.batch_iterator(fastq, batch_size=batch_size)
                    
                    for batch in batch_iter:
                        batch_len = len(batch)
                        yield batch, read_idx
                        read_idx += batch_len
            else:
                with open(fn) as handle:
                    fastq =\
                        (read_fastq(title, sequence, qscore) for title, sequence, qscore in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle))
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
    rst_futures = helper.multithreading_submit(batch_barcode_to_fastq, 
                           r_batches_with_idx, 
                           threads=n_process*4,
                           pbar_func=lambda x: len(x[0]),
                            assignment_df = assignment_df,
                            gz = gz)
    
    tmp_files = []
    for f in  rst_futures:
        tmp_fn = f.result()
        tmp_files.append(f.result())
    logger.info(f"Concatenating tmp fastq files to {fastq_out}")
    helper.concatenate_files(fastq_out, *tmp_files)
    logger.info(helper.green_msg(f"Demultiplexing completed!", printit = False))


# def main_single_thread_wrt(fastq_fns, fastq_out, putative_bc_csv, whitelsit_csv, n_process, gz, batchsize):
#     # check file exist. Avoid overwrite
#     i = 1
#     while os.path.exists(fastq_out):
#         fastq_out = f'{fastq_out.split("_")[0]}_{i}'
#         i+=1

#     output_handle = gzip.open(fastq_out, 'wt') if gz else open(fastq_out, 'w')
#     assignment_df = assign_barcodes(putative_bc_csv, whitelsit_csv, n_process)


#     r_batches = blaze.read_batch_generator(fastq_fns, batchsize)
    
#     read_idx = 0
#     for read_batch in tqdm(r_batches):
#         records = []
#         for r in read_batch:
#             row = assignment_df.iloc[read_idx]#the row in putative bc table
#             read_idx += 1 
            
#             try:
#                 assert row.read_id == r.id
#             except AssertionError:
#                 helper.err_msg("Different order in putative bc file and input fastq!")

#             if not row.putative_bc:
#                 continue
#             r = fastq_modify_header(r, row.BC_corrected, row.putative_umi)

#             if row.umi_end < 0:
#                 r = fastq_trim_seq(r,end = int(row.umi_end))
#             else:
#                 r = fastq_trim_seq(r,start = int(row.umi_end))
#             records.append(r)
#         SeqIO.write(records, output_handle, "fastq")
#     output_handle.close()

if __name__ == '__main__':
    #main_multi_thread()
    pass