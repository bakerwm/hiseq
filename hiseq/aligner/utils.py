#!/usr/bin/env python3

"""Functions for aligner
check fastq/index, arguments, ...

file existence 
arguments 

"""

import pyfastx
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import *
# from aligner_index import *


def check_fx_args(fq1, fq2=None, **kwargs):
    """Check the fastx in arguments
    
    Parameters
    ----------
    fq1 : str or list
    fq2 : None, str or list

    check:
    1. file exists
    2. file type
    3. fq paired
    4. check_empty
    """
    if isinstance(fq1, str):
        k1 = check_fx(fq1, **kwargs)
        if isinstance(fq2, str):
            k2 = check_fx(fq2, **kwargs)
            p1 = check_fx_paired(fq1, fq2, **kwargs)
            out = all([k1, k2, p1])
        elif fq2 is None:
            out = k1
        else:
            log.error('fq2 not valid, expect str or NoneType, got {}'.format(
                type(fq2).__name__))
            out = False
    elif isinstance(fq1, list):
        k1 = all([check_fx(i, **kwargs) for i in fq1])
        if isinstance(fq2, list):
            k2 = all([check_fx(i, **kwargs) for i in fq2])
            if len(fq1) == len(fq2):
                k3 = all([check_fx_paired(a, b) for a,b in zip(fq1, fq2)])
                out = all([k1, k2, k3])
            else:
                log.error('fq1, fq2 not in same length')
                out = False
        elif fq2 is None:
            out = k1
        else:
            log.error('fq2 not valid, expect list or NoneType, got {}'.format(
                type(fq2).__name__))
            out = False
    else:
        log.error('fq1 expect str or list, got {}'.format(
            type(fq1).__name__))
        out = False
    return out


# fq files
def check_fx(fx, **kwargs):
    """Check the fastq/a files
    1. file exist
    2. fq1 required
    3. fq type: fasta/q
    """
    kwargs['check_empty'] = kwargs.get('check_empty', True)
    kwargs['show_error'] = kwargs.get('show_error', False)
    if isinstance(fx, str):
        if check_file(fx, **kwargs):
            try:
                fx_type = Fastx(fx).format # fasta/q
                out = fx_type in ['fasta', 'fastq']
            except ValueError as err:
                log.info('Failed to read file, with error: {}'.format(err))
                out = False
#             finally:
#                 log.info('check_fx() ...')            
        else:
            if kwargs['show_error']:
                log.error('fx failed, {}'.format(fx))
            out = False
    else:
        out = False
    return out


def check_fx_paired(fq1, fq2, **kwargs):
    """Check the fq1 and fq2 are paired or not
    the name of fq1, fq2 are the same

    fq1: @the-name/1
    fq2: @the-name/2

    scan the first read from the two files
    """
    if isinstance(fq1, str) and isinstance(fq2, str):
        if check_fx(fq1, **kwargs) and check_fx(fq2, **kwargs):
            # check fq name:
            fx1 = pyfastx.Fastx(fq1)
            fx2 = pyfastx.Fastx(fq2)
            for a,b in zip(fx1, fx2):
                out = a[0][:-1] == b[0][:-1]
                break
        else:
            out = False
    else:
        out = False
    return out


def check_file(x, **kwargs):
    """Check the x file
    1. file exists
    
    Parameters
    ----------
    x : str
        Path to a file
    
    Keyword Parameters
    ------------------
    show_error : bool
        Show the error messages
        
    show_log : bool
        Show the log messages 
        
    check_empty : bool
        Check if the file is empty or not,  gzipped empty file, size=20
        
    emptycheck : bool
        see check_empty
    """
    show_error = kwargs.get('show_error', False)
    show_log = kwargs.get('show_log', False)
    check_empty = kwargs.get('check_empty', False)
    emptycheck = kwargs.get('emptycheck', False) # for old version
    if isinstance(x, str):
        if file_exists(x):
            x_size = os.stat(x).st_size
            # empty gzipped file, size=20
            q_size = 20 if x.endswith('.gz') else 0
            out = x_size > q_size if check_empty or emptycheck else True
            if show_log:
                flag = 'ok' if out else 'failed'
                log.info('{:<6s} : {}'.format(flag, x))
        else:
            if show_error:
                log.error('file not exists: {}'.format(x))
            out = False # failed
    elif isinstance(x, list):
        out = all([check_file(i, **kwargs) for i in x])
    else:
        log.error('x expect str or list, got {}'.format(type(x).__name__))
        out = False
    return out

