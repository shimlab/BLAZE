from collections import namedtuple
import argparse
from tqdm import tqdm
import os
import gzip
import numpy as np
import logging

import helper

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
class read_fastq:
    """This class is for mimic the Bio.SeqIO fastq record. The SeqIO is not directly used because it's slow.
    """
    def __init__(self,title, sequence, qscore, quality_map = False):
        self.id = title.split()[0].strip("@")
        self.seq = sequence
        self.qscore = qscore

def parse_command_line():
    parser = argparse.ArgumentParser(description='Fastq processing script')
    parser.add_argument('fastqs', nargs='+', help='Input fastq filenames or directory')
    parser.add_argument('fastq_out', help='Output fastq filename')
    # add optional arguments for whether to output gzipped fastq files
    parser.add_argument('--gzip_out', action='store_true', help='Output gzipped fastq files if set')
    parser.add_argument('--n_process', type=int, default=1, help='Number of processes to use (default: 1)')
    args = parser.parse_args()
    return args.fastqs, args.fastq_out, args.n_process, args.gzip_out


def polyT_trimming_idx(seq, reverse=False, umi_end_idx=0, polyT_end_minT=7, polyT_end_check_win=10):
    """sequencing after trimming the adaptor and UMI
    Args:
        seq (str): read sequence
        reverse (bool, optional): read is from transcript strand, Defaults to False.
        umi_end_idx (int, optional): end of UMI, Defaults to 0. Note that the index is related to polyT strand (i.e. '+')
        polyT_end_minT (int, optional): minimum number of T to call polyT, Defaults to 7.
        polyT_end_check_win (int, optional): window size to check polyT, Defaults to 10.
    Returns:
        index of the end of the polyT, negative if it's polyA strand
    """

    # take reverse complement if read is coming from transcript strand (with ployA instead ployT)

    assert polyT_end_minT < polyT_end_check_win, "polyT_end_minT should be smaller than polyT_end_check_win"
    if reverse:
        seq = helper.reverse_complement(seq)

    # find the first appearance of TTTT in seq
    polyT_start = seq.find('TTTT', umi_end_idx)
    if polyT_start == -1:
        return umi_end_idx
    
    read_code = np.array([int(x == 'T') for x in seq])
    for idx, nt in enumerate(read_code[polyT_start:]):
        if nt == 1:
            polyT_start += 1
        elif sum(read_code[polyT_start:polyT_start+polyT_end_check_win]) >= polyT_end_minT: 
            polyT_start += 1
        else:
            break 
    return int(polyT_start) if not reverse else int(-polyT_start)


def _read_batch_generator(fastq_fns, batch_size):
    """Generator of barches of reads from list of fastq files with the idx of the first read
    in each batch

    Args:
        fastq_fns (list): fastq filenames/or directory
        batch_size (int, optional):  Defaults to 100.
    """

    # if any element in fastq_fns is a directory, get all fastq files inside
    fastq_fns_processed = []
    for fn in fastq_fns:
        if os.path.isdir(fn):
            fastq_fns_processed += [os.path.join(fn, x) for x in os.listdir(fn) if x.endswith('.fastq') or x.endswith('.fastq.gz')]
        else:
            fastq_fns_processed.append(fn)

    fastq_fns = fastq_fns_processed

    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in helper.fastq_parser(handle))

                batch_iter = helper.batch_iterator(fastq, batch_size=batch_size)
                for batch in batch_iter:
                    yield batch
                
        else:
            with open(fn) as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in helper.fastq_parser(handle))
                batch_iter = helper.batch_iterator(fastq, batch_size=batch_size)

                for batch in batch_iter:
                    yield batch


def _proc_read_batches(read_batch, gz):
    out_buffer = ''
    # create read list
    for r in read_batch:
        strand = r.id.split('_')[-1].strip()
        polyT_end = polyT_trimming_idx(r.seq, reverse=strand=='-')
        if polyT_end < 0:
            seq = r.seq[:int(polyT_end)]
            qscore = r.qscore[:int(polyT_end)]
        else:
            seq = r.seq[int(polyT_end):]
            qscore = r.qscore[int(polyT_end):]
        
        out_buffer += '@' + r.id + '\n'
        out_buffer += str(seq) + '\n'
        out_buffer += '+\n'
        out_buffer += qscore + '\n'

    
    if gz:
        b_out_buffer = gzip.compress(out_buffer.encode('utf-8'))
    else:
        b_out_buffer = out_buffer.encode('utf-8')

    return b_out_buffer, len(read_batch)
    # logger.info(helper.green_msg(f"Demultiplexing finshied: ", printit = False))
    # logger.info(helper.green_msg(f"Successfully demultiplexed reads / Total reads: {sum(df.BC_corrected!='')} / {len(df.BC_corrected)}. ", printit = False))



def polyT_trimmer(fastq_fns, fastq_out, n_process, gz, batchsize):
    """Main function: Demultiplex fastq files using putative barcode csv and whitelist csv
    """
    # greating generator for read and putative barcode batches
    r_batches = \
        _read_batch_generator(fastq_fns, batchsize)
    # assign putative barcode to whitelist
    # single thread version
    if n_process == 1:
        with open(fastq_out, 'wb') as output_handle:
            pbar = tqdm(unit="Reads", desc='Processed')
            for r_batch in r_batches:
                b_fast_str, read_count = _proc_read_batches(r_batch,  gz)
                output_handle.write(b_fast_str)
                pbar.update(read_count)

        logger.info(helper.green_msg(f"PolyT trimming completed! saved in {fastq_out}!", printit = False))

    # multi-thread version
    else:
        rst_futures = helper.multiprocessing_submit(_proc_read_batches, 
                            r_batches, 
                            n_process=n_process,
                            schduler = "process",
                            pbar_func=len,
                            gz = gz)
        
        # collect results
        with open(fastq_out, 'wb') as output_handle:
            for f in rst_futures:
                b_fast_str, read_count = f.result()
                output_handle.write(b_fast_str)
        logger.info(helper.green_msg(f"PolyT trimming completed! saved in {fastq_out}!", printit = False))

if __name__ == '__main__':
    fastq_fns, fastq_out, n_process, gzip_out = parse_command_line()
    polyT_trimmer(fastq_fns, fastq_out, n_process, gz=gzip_out, batchsize=4000)
        

