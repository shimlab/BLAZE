import numpy as np
import pandas as pd
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
from pathlib import Path
from tqdm import tqdm
import os
import sys
import shutil
from collections import namedtuple


def reverse_complement(seq):
	'''
	Args: <str>
		queried seq
	Returns: <str>
		reverse_complement seq
	'''
	comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
					'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
	letters = \
		[comp[base] if base in comp.keys() else base for base in seq]
	return ''.join(letters)[::-1]

def err_msg(msg, printit = False):
    CRED = '\033[91m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def warning_msg(msg, printit = False):
    CRED = '\033[93m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def green_msg(msg, printit = False):
    CRED = '\033[92m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def bold_text(text, printit = False):
    if printit:
        print(f"\033[1m{text}\033[0m")
    else:
        return f"\033[1m{text}\033[0m"

def sliding_window_sum(array, window) :
    cum = np.cumsum(array)  
    return cum[window:] - cum[:-window]

def sliding_window_mean(array, window) :
    cum = np.cumsum(array)  
    return (cum[window:] - cum[:-window]) / window


class param:
    def __init__(self, **kwargs):
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
    def add(self, attr_name, attr_val, overwrite = True):
        if attr_name in __dict__.keys() and not overwrite:
            pass
        else:    
            self.__dict__[attr_name] = attr_val
    def rm(self, attr_name):
        if attr_name in self.__dict__.keys():
            del self.__dict__[attr_name]
    def __str__(self):
        return str(self.__dict__)
    
    def check(self, attr_list, add_none = True, silent = True):
        """
        Check whether the attributes in a given list are present. 

        Parameters
        ----------
        attr_list : LIST
            list of strings of the attributes to check 
        add_none : BOOL
            whether or not to create the missed attributes with value None
        silent : BOOL
            always return True if silent = True
        Returns
        -------
        True/False if all attributes is present
        """
        try:
            assert isinstance(add_none, bool)
            assert isinstance(silent, bool)
        except (AttributeError, TypeError):
            raise AssertionError("add_none and silent must be bool variable.")
        
        check_res = True
        for attr in attr_list:
            if attr not in self.__dict__.keys():
                check_res = False
                self.__dict__[attr] = None
        return check_res if not silent else True
    

def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,
                           pbar=True, pbar_unit='Read',pbar_func=len, 
                           schduler = 'process', *arg, **kwargs):
    """multiple processing or threading, 

    Args:
        func: function to be run parallely
        iterator: input to the function in each process/thread
        n_process (int, optional): number of cores or threads. Defaults to mp.cpu_count()-1.
        pbar (bool, optional): Whether or not to output a progres bar. Defaults to True.
        pbar_unit (str, optional): Unit shown on the progress bar. Defaults to 'Read'.
        pbar_func (function, optional): Function to calculate the total length of the progress bar. Defaults to len.
        schduler (str, optional): 'process' or 'thread'. Defaults to 'process'.

    Yields:
        return type of the func: the yield the result in the order of submit
    """
    class fake_future:
        # a fake future class to be used in single processing
        def __init__(self, rst):
            self.rst = rst
        def result(self):
            return self.rst

    if schduler == 'process':
        # make sure the number of process is not larger than the number of cores
        n_process = min(n_process-1, mp.cpu_count()-1)
        if n_process > 1:
            executor = concurrent.futures.ProcessPoolExecutor(n_process)
    elif schduler == 'thread':
        if n_process > 1:
            executor = concurrent.futures.ThreadPoolExecutor(n_process)
    else:
        green_msg('Error in multiprocessing_submit: schduler should be either process or thread', printit=True)
        sys.exit(1)

    if pbar:
        _pbar = tqdm(unit=pbar_unit, desc='Processed')
        
    # run in single process/thread if n_process < 1
    if n_process <= 1:
        for it in iterator:
            yield fake_future(func(it, *arg, **kwargs))
            if pbar:
                _pbar.update(pbar_func(it))
        return

    # A dictionary which will contain the future object
    max_queue = n_process
    futures = {}
    n_job_in_queue = 0
    
    # make sure the result is yield in the order of submit.
    job_idx = 0
    job_completed = {}

    # submit the first batch of jobs
    while n_job_in_queue < max_queue:
        i = next(iterator, None)
        if i is None:
            break
        futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
        job_idx += 1
        n_job_in_queue += 1
        job_to_yield = 0
    # yield the result in the order of submit and submit new jobs
    while True:
        # will wait until as least one job finished
        # batch size as value, release the cpu as soon as one job finished
        job = next(as_completed(futures), None)

        # yield the completed job in the order of submit  
        if job is not None:
            job_completed[futures[job][1]] = job, futures[job][0]
            del futures[job]

        # 
        if job is None and i is None and len(job_completed)==0:
            break

        # check order
        while job_to_yield in job_completed.keys():
            # update pregress bar based on batch size
            if pbar:
                _pbar.update(job_completed[job_to_yield][1])
            yield job_completed[job_to_yield][0]
            del job_completed[job_to_yield]
            
            # submit new job
            i = next(iterator, None)
            if i is not None:
                futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
                job_idx += 1
                
            job_to_yield += 1

# multiproces panda data frame  
def procee_batch(df,  row_func, *arg, **kwargs):
        return df.apply(row_func, axis=1, *arg, **kwargs)

def df_multiproceccing_apply(df, func, n_process, aggr_func=pd.concat, pbar = True, pbar_unit='Read',pbar_func=len,*arg, **kwargs):
    """This is a re-implementation which is very similar to dask apply. 
    """
    def sub_df_generator(df, n_part):
        """ 
        Generator: split a large df by row into multiple 
        """
        sub_size = int(np.ceil(len(df) / n_part))
        batch_start, batch_end = 0, sub_size
        while batch_start < len(df):
            yield df.iloc[batch_start:batch_end,].copy()
            batch_start += sub_size
            batch_end += sub_size

    #run 100 batch per process on average but ensure at least 1000 reads per batch
    num_batch = min(int(len(df)/1000)+1, n_process*100) 
    df_iter = sub_df_generator(df, n_part=num_batch)
    
    rst_futures = multiprocessing_submit(procee_batch, df_iter, n_process=n_process ,
                           pbar = pbar, pbar_unit=pbar_unit,pbar_func=pbar_func, 
                           schduler = 'process', row_func = func, *arg, **kwargs)

    for i in rst_futures:
        yield i.result()
    # rsts = [x.result() for x in rst_futures]
    # return aggr_func(rsts)
    
# concatenate multiple files
def concatenate_files(output_file, *input_files):
    with open(output_file, 'wb') as outfile:
        for input_file in input_files:
            with open(input_file, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
            os.remove(input_file)    

# get file with a certian extensions
def get_files_by_suffix(search_dir, suffix, recursive=True):
    files = []
    if isinstance(suffix, str):
        suffix = [suffix]
    if recursive:
        for i in suffix:
            files.extend(Path(search_dir).rglob(i))
        return sorted(files)
    else:
        for i in suffix:
            files.extend(Path(search_dir).glob(i))
        return sorted(files)

# check file exist. Exit if not
def check_files_exist(file_list):
    if isinstance(file_list, str):
        file_list = [file_list]
    exit_code = 0
    for fn in file_list:
        if not os.path.exists(fn):
            exit_code = 1
            err_msg(f'Error: can not find {fn}', printit=True)
    if exit_code == 1:
        sys.exit(1)
    else:
        return True

# split any iterator in to batches  
def batch_iterator(iterator, batch_size):
    """generateor of batches of items in a iterator with batch_size.
    """
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(entry)
        
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):
        yield batch

# a light class for a read in fastq file
read_tuple = namedtuple('read_tuple', ['id', 'seq', 'q_letter'])
def fastq_parser(file_handle):
    while True:
        id = next(file_handle, None)
        if id is None:
            break
        seq = next(file_handle)
        next(file_handle) # skip  '+'
        q_letter = next(file_handle)
        yield read_tuple(id[1:].split()[0], seq.strip(), q_letter.strip())
        


# validate filename (check suffix):
def check_suffix(filename, suffix_lst):
    """check the suffix of a filename
    Args:
        filename (str)
        suffix_lst (list/str)
    """
    # check output filenames
    if isinstance(suffix_lst, str):
        return filename.endswith(suffix_lst)
    if any([filename.endswith(x) for x in suffix_lst]):
        return True
    else: 
        return False

