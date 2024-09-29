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
from io import StringIO
import io
import gzip


import blaze.helper as helper
from blaze.config import *


# setup logging
LOG_FORMAT = '(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

class read_fastq:
    """This class is for mimic the Bio.SeqIO fastq record. The SeqIO is not directly used because it's slow.
    """
    def __init__(self,title, sequence, qscore, quality_map = False):
        self.id = title.split()[0].strip("@")
        self.seq = sequence
        self.qscore = qscore

def _match_bc_row(row, whitelist, max_ed, minQ):
    """
    Args:
        row (row in pandas.DataFrame)
        whitelist (set): list of bc to assign to
        max_ed (int): maximum edit distance allowed for a assignment
        minQ (int): minimum quality score of putative BC for read assignment
    Returns:
        1. assigned barcode <str>:  if a unambiguous barcode was found in the whitelist
            '': if no barcode was found in the whitelist
            'ambiguous': if multiple barcode was found in the whitelist
        2. adjusted putative_umi
        3. strand '+' for read with transcript strand (End with polyA)
                  '-' for read with negative strand (start with BC,umi,polyT...)

    """
    if not row.polyT_end:
        strand = ''
    elif row.polyT_end > 0:
        strand = '-'
    else: 
        strand = '+'

    if minQ and row.putative_bc_qscore < minQ:
        return ['', '', '']

    if not row.putative_bc or row.putative_bc in whitelist:
        return [row.putative_bc, row.putative_umi, strand]
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
                    return ['', row.putative_umi, strand]
    
    if bc_hit == 'ambiguous' or bc_hit == '':
        return ['', row.putative_umi, strand]
    else:
        pass

    # adjust the umi start position
    umi_adj = bc_hit_end_idx - (len(bc) - 1 -DEFAULT_ED_FLANKING )
    out_umi = row.putative_umi
    if umi_adj > 0:
        out_umi = row.putative_umi[umi_adj:] + row.post_umi_flanking[:umi_adj]
    elif umi_adj < 0:
        out_umi =  row.putative_bc[umi_adj:] + row.putative_umi[:umi_adj]
    return [bc_hit, out_umi, strand]
            
def batch_barcode_to_fastq(r_batches_with_idx, assignment_df ,gz = True):
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
            helper.err_msg("Different order in putative bc file and input fastq!", printit = True)
            sys.exit()

        if not row.BC_corrected:
            continue

        if row.polyT_end < 0:
            seq = r.seq[:int(row.polyT_end)]
            qscore = r.qscore[:int(row.polyT_end)]
        else:
            seq = r.seq[int(row.polyT_end):]
            qscore = r.qscore[int(row.polyT_end):]
        
        out_buffer += f"@{row.BC_corrected}_{row.putative_umi}#{row.read_id}_{row.strand}\tCB:Z:{row.BC_corrected}\tUB:Z:{row.putative_umi}\n"
        out_buffer += str(seq) + '\n'
        out_buffer += '+\n'
        out_buffer += qscore + '\n'
    
    output_handle.write(out_buffer)
    output_handle.close()
    return temp_file.name

def _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batch_size):
    """Generator of barches of reads from list of fastq files with the idx of the first read
    in each batch

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 100.
    """
    read_idx = 0
    putative_bc_f = open(putative_bc_csv, 'r')
    putative_bc_header = next(putative_bc_f)

    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in helper.fastq_parser(handle))

                batch_iter = helper.batch_iterator(fastq, batch_size=batch_size)
                
                for batch in batch_iter:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO(
                            putative_bc_header + \
                            ''.join([next(putative_bc_f) for x in range(batch_len)])
                        ))
                    yield batch, read_idx, batch_bc_df
                    read_idx += batch_len
        else:
            with open(fn) as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in helper.fastq_parser(handle))
                read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO( putative_bc_header + \
                        ''.join([next(putative_bc_f) for x in range(batch_len)]))
                    )

                    yield batch, read_idx, batch_bc_df
                    read_idx += batch_len
    putative_bc_f.close()

    
def _assign_read_batches(r_batch, whitelist, max_ed, gz, restrand ,minQ=0):
    """Single-thread function:
        Assign all putative barcode to whiteliested barcode

    Args:
        r_batch (tuple): yield from read_and_bc_batch_generator_with_idx 
                                            (read_batch, start_idx, batch_bc_df)
        whitelist (list): list of barcode 
        n_process (int): number of process 
        max_ed (int): maximum edit distance allowed for a assignment
        gz (bool): output fastq is gzipped or not
        minQ (int): minimum quality score for read assignment

    Returns:
        1. pd.dataframe: Same as the input csv, with an additional column called 'BC_corrected', 'putative_umi', 'strand'
        2. binary str: output fastq string
        3. int: number of reads successfully assigned to a barcode
        4. total number of reads in the batch
    """
    # read putative barcode
    read_batch, start_df_idx, df = r_batch
    df = df.fillna('')
    whitelist = set(whitelist)
    out_buffer = ''
    
    new_cols = []
    for row in df.itertuples():
        new_cols.append(_match_bc_row(row, whitelist, max_ed, minQ))

    df[['BC_corrected','putative_umi', 'strand']] = new_cols
    demul_read_count = sum(df.BC_corrected!='')

    # create read list

    for r, bc in zip(read_batch, df.itertuples()):
        try:
            assert bc.read_id == r.id            
        except AssertionError:
            helper.err_msg("Different order in putative bc file and input fastq!", printit = True)
            sys.exit()

        if not bc.BC_corrected:
            continue
        if bc.polyT_end < 0:
            seq = r.seq[:int(bc.polyT_end)]
            qscore = r.qscore[:int(bc.polyT_end)]
        else:
            seq = r.seq[int(bc.polyT_end):]
            qscore = r.qscore[int(bc.polyT_end):]
        
        if restrand and bc.strand == '-':
            seq = helper.reverse_complement(seq)
            qscore = qscore[::-1]
        
        # write to fastq
        out_buffer += f"@{bc.BC_corrected}_{bc.putative_umi}#{bc.read_id}_{bc.strand}\tCB:Z:{bc.BC_corrected}\tUB:Z:{bc.putative_umi}\n"
        out_buffer += str(seq) + '\n' 
        out_buffer += '+\n'
        out_buffer += qscore + '\n'
    
    if gz:
        b_out_buffer = gzip.compress(out_buffer.encode('utf-8'))
    else:
        b_out_buffer = out_buffer.encode('utf-8')

    return df, b_out_buffer, demul_read_count, len(read_batch)
    # logger.info(helper.green_msg(f"Demultiplexing finshied: ", printit = False))
    # logger.info(helper.green_msg(f"Successfully demultiplexed reads / Total reads: {sum(df.BC_corrected!='')} / {len(df.BC_corrected)}. ", printit = False))

def assign_read(fastq_fns=None, fastq_out=None, putative_bc_csv=None, 
                    whitelsit_csv=None, max_ed=None, n_process=None,
                    batchsize=None, minQ=0, restrand=True, args=None):
    """Main function: Demultiplex fastq files using putative barcode csv and whitelist csv
        Input:
            fastq_fns (list): list of fastq filenames
            fastq_out (str): output fastq filename
            putative_bc_csv (str): putative barcode csv filename
            whitelsit_csv (str): whitelist csv filename
            max_ed (int): maximum edit distance allowed for a assignment
            n_process (int): number of process 
            gz (bool): output fastq is gzipped or not
            batchsize (int): batch size for read assignment
            minQ (int): minimum quality score for read assignment (NOTE: the minQ specified in args is not used as that is for whitelisting)
    """
    # check input
    if args:
        fastq_fns = args.fastq_fns
        fastq_out = args.output_fastq
        putative_bc_csv = args.out_raw_bc_fn
        whitelsit_csv = args.out_whitelist_fn
        max_ed = args.max_edit_distance
        n_process = args.threads
        batchsize = args.batch_size
        restrand = args.restrand

    gz = fastq_out.endswith('.gz')

    # greating generator for read and putative barcode batches
    r_batches = \
        _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batchsize)
    
    # read the whitelist
    whitelist = [] 
    with open(whitelsit_csv, 'r') as f:
        for line in f:
            whitelist.append(line.split('-')[0].strip())
    
    # assign putative barcode to whitelist
    # single thread version
    if n_process == 1:
        demul_count_tot = 0
        count_tot = 0
        with open(fastq_out, 'wb') as output_handle:
            pbar = tqdm(unit="Reads", desc='Processed')
            for r_batch in r_batches:
                _, b_fast_str, demul_count, read_count = _assign_read_batches(r_batch, whitelist, max_ed,  gz, restrand)
                demul_count_tot += demul_count
                count_tot += read_count
                output_handle.write(b_fast_str)
                pbar.update(read_count)

        logger.info(helper.green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!", printit = False))

    # multi-thread version
    else:
        rst_futures = helper.multiprocessing_submit(_assign_read_batches, 
                            r_batches, 
                            n_process=n_process,
                            schduler = "process",
                            pbar_func=lambda x: len(x[0]),
                            whitelist = whitelist,
                            max_ed = max_ed,
                            gz = gz,
                            restrand = restrand)
        
        # collect results
        demul_count_tot = 0
        count_tot = 0
        with open(fastq_out, 'wb') as output_handle:
            for f in rst_futures:
                df, b_fast_str, demul_count, read_count = f.result()
                demul_count_tot += demul_count
                count_tot += read_count
                output_handle.write(b_fast_str)
        logger.info(helper.green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!", printit = False))

    return demul_count_tot, count_tot

if __name__ == '__main__':
    #assign_read()
    pass