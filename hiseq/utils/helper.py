
# -*- coding: utf-8 -*-

"""
Common functions for pipeline construction
file modification,
...
"""


from .bam import *
from .bed import *
from .download import *
from .fastx import *
from .featurecounts import *
from .file import *
from .utils import *
## tmp file
from .file import list_file as listfile
from .file import list_dir as listdir
from .file import fx_name as fq_name
from .file import symlink_file as file_symlink
from .file import copy_file as file_copy
from .file import remove_file as file_remove
from .file import file_is_gzipped as is_gz
from .file import file_nrows as file_row_counter 
from .file import gzip_file as gzip_cmd
from .bam import is_sam_flag as sam_flag_check
from .utils import update_obj


# import os
# import sys
# import re
# import gzip
# import shutil
# import signal
# import json
# import yaml
# import toml
# import numbers
# import pickle
# import fnmatch
# import tempfile
# import logging
# import functools
# import subprocess
# import collections
# import numpy
# import pysam
# import pybedtools
# import pathlib
# import binascii
# import datetime
# import pandas as pd
# import Levenshtein
# from itertools import combinations
# # from .args import args_init, ArgumentsInit
# ### local test ###
# # from args import args_init, ArgumentsInit # for local test
# 
# 
# def index_checker(seq_list, mm=0):
#     """
#     Check if the index compatiable, with specific mimatch
#     could be differenct in length
#     only consider mutations
#     """
#     item_pairs = list(combinations(seq_list, 2))
# 
#     tag = 0
#     for (a, b) in item_pairs:
#         a = a[:len(b)]
#         b = b[:len(a)] # the same length
#         m = Levenshtein.distance(a, b)
#         if m > mm:
#             print(m, a, b)
#             tag += 1
# 
#     return tag == 0
# 
# logging.basicConfig(
#     format='[%(asctime)s %(levelname)s] %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     stream=sys.stdout)
# log = logging.getLogger(__name__)
# log.setLevel('INFO')
# 
# 
# # Decorator:
# class Logger(object):
#     def __init__(self, level='INFO'):
#         logging.basicConfig(
#             format='[%(asctime)s %(levelname)s] %(message)s',
#             datefmt='%Y-%m-%d %H:%M:%S',
#             stream=sys.stdout)
#         self.logger = logging.getLogger(__name__)
#         self.logger.setLevel(level)
# 
#     def __call__(self, fn):
#         @functools.wraps(fn)
#         def decorated(*args, **kwargs):
#             try:
#                 self.logger.info('{0} - {1} - {2}'.format(
#                     fn.__name__,
#                     args,
#                     kwargs))
#                 result = fn(*args, **kwargs)
#                 self.logger.info(result)
#                 # return result
#             except Exception as ex:
#                 self.logger.info('Exception {0}'.format(ex))
#                 raise ex
#             return result
#         return decorated
# 
# 
# 
# def check_fq(fq):
#     """
#     Make sure
#     fq: str or list, or None
#     """
#     if fq is None:
#         # raise ValueError('fq1 required, got None')
#         pass
#     elif isinstance(fq, list):
#         fq = file_abspath(fq)
#     elif isinstance(fq, str):
#         fq = [file_abspath(fq)]
#     else:
#         log.error('fq failed, Nont, str, list expected, got {}'.format(type(fq).__name__))
# 
#     return fq
# 
# 
# def fq_paired(fq1, fq2):
#     """
#     Make sure fq1 and fq2, proper paired, by name: R1, R2
#     """
#     fq1 = check_fq(fq1)
#     fq2 = check_fq(fq2)
# 
#     if isinstance(fq1, str) and isinstance(fq2, str):
#         return distance(fq1, fq2) == 1
#     elif isinstance(fq1, list) and isinstance(fq2, list):
#         return [distance(i, j) == 1 for i, j in zip(fq1, fq2)]
#     else:
#         log.warning('fq not paired: {}, {}'.format(fq1, fq2))
#         return False
# 
# 
# ## 1. files and path ##
# def listdir(path, full_name=True, recursive=False, include_dir=False):
#     """
#     List all the files within the path
#     """
#     out = []
#     for root, dirs, files in os.walk(path):
#         if full_name:
#             dirs = [os.path.join(root, d) for d in dirs]
#             files = [os.path.join(root, f) for f in files]
#         out += files
# 
#         if include_dir:
#             out += dirs
# 
#         if recursive is False:
#             break
# 
#     return sorted(out)
# 
# 
# def listfile(path='.', pattern='*', full_name=True, recursive=False,
#     include_dir=False):
#     """
#     Search files by the pattern, within directory
#     fnmatch.fnmatch()
# 
#     pattern:
# 
#     *       matches everything
#     ?       matches any single character
#     [seq]   matches any character in seq
#     [!seq]  matches any char not in seq
# 
#     An initial period in FILENAME is not special.
#     Both FILENAME and PATTERN are first case-normalized
#     if the operating system requires it.
#     If you don't want this, use fnmatchcase(FILENAME, PATTERN).
# 
#     example:
#     listfile('./', '*.fq')
#     """
#     fn_list = listdir(path, full_name, recursive, include_dir=include_dir)
#     fn_list = [f for f in fn_list if fnmatch.fnmatch(os.path.basename(f), pattern)]
#     return sorted(fn_list)
# 
# 
# def list_fx(path, recursive=False):
#     """List the fastq, fasta files
#     *.fq
#     *.fastq
#     *.fa
#     *.fasta
#     """
#     f_list = listfile(path, recursive=recursive)
#     # re for fa/fq files
#     p = re.compile('\.f(ast)?(a|q)(\.gz)?', flags=re.IGNORECASE)
#     fx_list = [i for i in f_list if p.search(i)]
# 
#     return fx_list
# 
# 
# def list_fq_files(path, pattern='*'):
#     """
#     Parse fastq files within path, using the prefix
#     PE reads
#     SE reads
# 
#     *.fastq
#     *.fq
#     *.fastq.gz
#     *.fq.gz
# 
#     _1.
#     _2.
#     """
#     # all fastq files: *f[astq]+(.gz)?
#     fq_list = listfile(path, '*q.gz') # *fastq.gz, *fq.gz
#     fq_list.extend(listfile(path, '*q')) # *fastq, *fq
# 
#     # filter
#     if pattern == '*':
#         hit_list = fq_list
#     else:
#         p = re.compile(r'(_[12])?.f(ast)?q(.gz)?$')
#         hit_list = [f for f in fq_list if p.search(os.path.basename(f))]
# 
#     # chk1
#     if len(hit_list) == 0:
#         log.error('no fastq files found: {}'.format(path))
# 
#     # determine SE or PE
#     r0 = r1 = r2 = []
#     for i in hit_list:
#         p1 = re.compile(r'_[rR]?1.f(ast)?q(.gz)?$') # read1
#         p2 = re.compile(r'_[rR]?2.f(ast)?q(.gz)?$') # read2
#         if p1.search(i):
#             r1.append(i)
#         elif p2.search(i):
#             r2.append(i)
#         else:
#             r0.append(i)
# 
#     # chk2
#     if len(r2) > 0 and not len(r1) == len(r2):
#         log.error('read1 and read2 not equal: \nread1: {}\nread2: {}'.format(r1, r2))
# 
#     # organize [[r1, r2], [r1, r2], ...]
#     out_list = []
#     for i, j in zip(r1, r2):
#         out_list.append([i, j])
# 
#     for i in r0:
#         out_list.append([i, None])
# 
#     # if len(r1) > 6:
#     #     log.warning('too many records matched : {} \n{}'.format(pattern, r1))
# 
#     return out_list
# 
# 
# def is_gz(filepath):
#     if os.path.exists(filepath):
#         with open(filepath, 'rb') as test_f:
#             return binascii.hexlify(test_f.read(2)) == b'1f8b'
#     else:
#         if filepath.endswith('.gz'):
#             return True
#         else:
#             return False
# 
# 
# def file_prefix(fn, with_path=False):
#     """
#     extract the prefix of a file
#     remove extensions
#     .gz, .fq.gz
#     """
# #     assert isinstance(fn, str)
# #     p1 = os.path.splitext(fn)[0]
# #     px = os.path.splitext(fn)[1]
# #     if px.endswith('gz') or px.endswith('.bz2'):
# #         px = os.path.splitext(p1)[1] + px
# #         p1 = os.path.splitext(p1)[0]
# #     if not with_path:
# #         p1 = os.path.basename(p1)
# #     return [p1, px]
#     if isinstance(fn, str):
#         if fn.endswith('.gz') or fn.endswith('.bz2'):
#             out = os.path.splitext(fn)[0]
#         out = os.path.splitext(out)[0]
#         if not with_path:
#             out = os.path.basename(out)
#     elif isinstance(fn, list):
#         out = [file_prefix(i, with_path) for i in fn]
#     elif fn is None:
#         out = None
#     else:
#         log.error('unknown fn, str,list,None expected, got {}'.format(
#             type(fn).__name__))
#     return out
# 
# 
# def check_file(x, show_log=False, emptycheck=False):
#     """
#     if x (file, list) exists or not, empty
# 
#     Whether the file is empty
#     see: https://stackoverflow.com/a/2507871/2530783
#     os.stat("file").st_size == 0    
#     """
#     if isinstance(x, str):
#         flag = 'ok' if os.path.exists(x) else 'failed'
#         if flag == 'ok' and emptycheck:
#             flag = 'ok' if os.stat(x).st_size > 20 else 'failed' # in case of gzipped empty file (size=20)
#         if show_log is True:
#             log.info('{:<6s} : {}'.format(flag, x))
#         return flag == 'ok' # os.path.exists(x)
#     elif isinstance(x, list):
#         return all([check_file(i, show_log, emptycheck) for i in x])
#     else:
#         log.warning('x, str and list expected, {} got'.format(type(x).__name__))
#         return None
# 
# 
# def check_path(x, show_log=False, create_dirs=True):
#     """
#     Check if x is path, Create path
#     """
#     if isinstance(x, str):
#         if os.path.isdir(x):
#             tag = True
#         else:
#             if create_dirs is True:
#                 try:
#                     os.makedirs(x)
#                     tag = True
#                 except:
#                     tag = False
#             else:
#                 tag = False
#         # show log
#         flag = 'ok' if tag is True else 'failed'
#         if show_log is True:
#             log.info('{:<6s} : {}'.format(flag, x))
#         return tag
#     elif isinstance(x, list):
#         return all([check_path(i, show_log, create_dirs) for i in x])
#     else:
#         log.warning('expect str and list, not {}'.format(type(x)))
#         return None
# 
# 
# def fq_name(fq, include_path=False, pe_fix=False):
#     """
#     parse the name of fastq file:
#     .fq.gz
#     .fastq.gz
#     .fq
#     .fastq
#     (also for fasta, fa)
#     """
#     # if isinstance(x, str):
#     #     fname = file_prefix(x)[0]
#     #     fname = re.sub('[._][rR]?1$', '', fname)
#     #     fname = re.sub('_\d$', '', fname)
#     # elif isinstance(x, list):
#     #     fname = [fq_name(f) for f in x]
#     # elif x is None:
#     #     log.warning('None type detected')
#     #     fname = None
#     # else:
#     #     log.warning('unknown input detected')
#     #     fname = None
#     # return fname
#     ###############
#     # fastq.gz, fastq, fasta.gz, fa.gz
#     # p1 = re.compile('(_[12])?[.](fast|f)[aq](.gz)?$', re.IGNORECASE)
#     p1 = re.compile('[.](fast|f)[aq](.gz)?$', re.IGNORECASE)
#     p2 = re.compile('[._][12]$', re.IGNORECASE)
#     if isinstance(fq, str):
#         fq = fq if include_path is True else os.path.basename(fq)
#         fq_tmp = re.sub(p1, '', fq) # r1_1.fq.gz : r1_1
#         if pe_fix is True:
#             fq_tmp = re.sub(p2, '', fq_tmp) # r1_1.fq.gz: r1
#         return fq_tmp
#     elif isinstance(fq, list):
#         return [fq_name(x, include_path=include_path, pe_fix=pe_fix) for x in fq]
#     else:
#         log.warning('unknown type found: {}'.format(type(fq)))
#         return fq
# 
# 
# def fq_name_rmrep(fq, include_path=False, pe_fix=True):
#     """
#     x, filename, or list
# 
#     Remove the *.rep[123], *.REP[123] from tail
#     """
#     # if isinstance(x, str):
#     #     return fq_name(x).rstrip('rep|REP|r|R||_|.|1|2')
#     # elif isinstance(x, list):
#     #     return [fq_name_rmrep(i) for i in x]
#     # else:
#     #     pass
#     ##################
#     p1 = re.compile('[._](rep|r)[0-9]+$', re.IGNORECASE) # _rep1, _r1
#     if isinstance(fq, str):
#         return re.sub(p1, '', fq_name(fq, include_path, pe_fix))
#     elif isinstance(fq, list):
#         return [fq_name_rmrep(x, include_path, pe_fix) for x in fq]
#     else:
#         log.warning('unknown type found: {}'.format(type(fq)))
#         return fq
# 
# 
# def file_abspath(file):
#     """
#     Create os.path.abspath() for files, directories
#     """
#     if file is None:
#         return None
#     elif isinstance(file, str):
#         return os.path.abspath(file)
#     elif isinstance(file, list):
#         return [file_abspath(x) for x in file]
#     else:
#         log.warning('unknown type found: {}'.format(type(file)))
#         return file
# 
# 
# def file_exists(file, isfile=True):
#     """
#     os.path.exists() for files, os.path.isfile()
#     """
#     if file is None:
#         return None
#     elif isinstance(file, str):
#         chk1 = os.path.exists(file) 
#         chk2 = os.path.isfile(file)
#         chk3 = os.path.isdir(file) # candidate
#         return chk1 and chk2 if isfile else chk1
#     elif isinstance(file, list):
#         return [file_exists(x, isfile) for x in file]
#     else:
#         log.warning('unknown type found: {}'.format(type(file)))
#         return file
# 
# 
# def file_symlink(src, dest, absolute_path=False):
#     """
#     Create symlink for dest files
#     
#     file, directories
#     """
#     if isinstance(src, str) and isinstance(dest, str):
#         src = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
#         dest = os.path.abspath(os.path.expanduser(os.path.expandvars(dest)))
# 
#         if os.path.isdir(src):
#             # target exists
#             if os.path.exists(dest):
#                 log.info('symlink() skipped, dest exists: {}'.format(dest))
#             else:
#                 if absolute_path:
#                     os.symlink(src, dest)
#                 else:
#                     rel_src = os.path.relpath(src, dest)
#                     os.symlink(rel_src, dest)
# 
#         elif os.path.isfile(src):
#             src_filename = os.path.basename(src)
#             src_dirname = os.path.dirname(src)
# 
#             if os.path.isdir(dest):
#                 dest_dirname = dest
#                 dest = os.path.join(dest_dirname, src_filename)
#             else:
#                 dest_dirname = os.path.dirname(dest)
# 
#             if os.path.exists(dest):
#                 log.info('symlink() skipped, dest exists: {}'.format(dest))
#             elif absolute_path:
#                 os.symlink(src, dest)
#             else:
#                 rel_src_dirname = os.path.relpath(src_dirname, dest_dirname)
#                 rel_src = os.path.join(rel_src_dirname, src_filename)
#                 os.symlink(rel_src, dest)
# 
#         else:
#             log.error('src path should be file or directory')
# 
#     else:
#         log.error('src path should be string, got: {}'.format(
#             type(src).__name__))
# 
# 
# def file_remove(x, ask=True):
#     """
#     Remove files, directories
#     """
#     del_list = []
#     undel_list = []
#     if isinstance(x, str):
#         if os.path.isfile(x) and file_exists(x):
#             del_list.append(x)
#         else:
#             undel_list.append(x)
#     elif isinstance(x, list):
#         for f in x:
#             if isinstance(f, str):
#                 if os.path.isfile(f) and file_exists(f):
#                     del_list.append(f)
#                 else:
#                     undel_list.append(f)
#             else:
#                 undel_list.append(f)
#     elif isinstance(x, dict):
#         for f in list(x.values()):
#             if isinstance(f, str):
#                 if os.path.isfile(f) and file_exists(f):
#                     del_list.append(f)
#                 else:
#                     undel_list.append(f)
#             else:
#                 undel_list.append(f)
#     else:
#         log.info('Nothing removed, str, list, dict expected, got {}'.format(
#             type(x).__name__))
# 
#     # remove files
#     if len(del_list) > 0:
#         del_msg = ['{:>6s}: {}'.format('remove', i) for i in del_list]
#         undel_msg = ['{:>6s}: {}'.format('skip', i) for i in undel_list]
#         msg = '\n'.join(del_msg + undel_msg)
#         log.info('Removing files: \n' + msg)
#         
#         if ask:
#             ask_msg = input('Removing the files? [Y|n]： ')
#         else:
#             ask_msg = 'Y'
# 
#         # remove
#         if ask_msg.lower() in ['y', 'yes']:
#             for f in del_list:
#                 os.remove(f)
#             log.info('{} files removed'.format(len(del_list)))
#         else:
#             log.info('Nothing removed, skipped')
# 
# 
# def file_copy(src, dest, force=False):
#     """
#     copy files to dest
#     """
#     if src is None:
#         log.error('file_copy() skipped, src is NoneType')
#     elif os.path.isfile(src):
#         if isinstance(dest, str):
#             if os.path.isdir(dest):
#                 dest_file = os.path.join(dest, os.path.basename(src))
#                 shutil.copy(src, dest_file)
#             elif os.path.isfile(dest):
#                 if force:
#                     shutil.copy(src, dest)
#                 else:
#                     log.warning('file_copy() skipped, dest exists: {}'.format(dest))
#             elif os.path.exists(os.path.dirname(dest)):
#                 shutil.copy(src, dest)
#             else:
#                 log.error('file_copy() skipped, dest is not file or dir')
#     else:
#         log.error('file_copy() skipped, src is not file')
# 
# 
# def path_remove(x, ask=True):
#     """
#     Remove directories
#     """
#     del_list = []
#     undel_list = []
#     if isinstance(x, str):
#         if os.path.isdir(x) and file_exists(x):
#             del_list.append(x)
#         else:
#             undel_list.append(x)
#     elif isinstance(x, list):
#         for f in x:
#             if os.path.isdir(f) and file_exists(x):
#                 del_list.append(f)
#             else:
#                 undel_list.append(f)
#     elif isinstance(x, dict):
#         for f in list(x.values()):
#             if os.path.isdir(f) and file_exists(x):
#                 del_list.append(f)
#             else:
#                 undel_list.append(f)
#     else:
#         log.info('Nothing removed, str, list, dict expected, got {}'.format(
#             type(x).__name__))
# 
#     # remove files
#     if len(del_list) > 0:
#         del_msg = ['{:>6s}: {}'.format('remove', i) for i in del_list]
#         undel_msg = ['{:>6s}: {}'.format('skip', i) for i in undel_list]
#         msg = '\n'.join(del_msg + undel_msg)
#         log.info('Removing files: \n' + msg)
#         
#         if ask:
#             ask_msg = input('Removing the files? [Y|n]： ')
#         else:
#             ask_msg = 'Y'
# 
#         # remove
#         if ask_msg.lower() in ['y', 'yes']:
#             for f in del_list:
#                 # os.remove(f)
#                 shutil.rmtree(f)
#             log.info('{} files removed'.format(len(del_list)))
#         else:
#             log.info('Nothing removed, skipped')
# 
# 
# def update_obj(obj, d, force=True, remove=False):
#     """
#     d: dict
#     force: bool, update exists attributes
#     remove: bool, remove exists attributes
#     Update attributes from dict
#     force exists attr
#     """
#     # fresh start
#     if remove is True:
#         for k in obj.__dict__:
#             delattr(obj, k)
#     # add attributes
#     if isinstance(d, dict):
#         for k, v in d.items():
#             if not hasattr(obj, k) or force:
#                 setattr(obj, k, v)
# 
#     return obj
# 
# 
# def file_row_counter(fn):
#     """
#     count the file rows
#     count '\n' 
#     from @glglgl on stackoverflow, modified
#     https://stackoverflow.com/a/9631635/2530783
#     """
#     def blocks(files, size = 1024 * 1024):
#         while True:
#             b = files.read(size)
#             if not b: break
#             yield b
#     freader = gzip.open if is_gz(fn) else open
#     with freader(fn, 'rt', encoding="utf-8", errors='ignore') as fi:
#         return sum(bl.count('\n') for bl in blocks(fi))
# 
# 
# def merge_names(x):
#     """
#     Get the name of replictes
#     common in left-most
#     """
#     assert isinstance(x, list)
#     name_list = [os.path.basename(i) for i in x]
#     name_list = [re.sub('.rep[0-9].*$', '', i) for i in name_list]
#     return list(set(name_list))[0]
# 
# 
# ## 2. commandline ##
# def run_shell_cmd(cmd):
#     """This command is from 'ENCODE-DCC/atac-seq-pipeline'
#     https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py
# 
#     save log to file
#     """
#     p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
#         stdin=subprocess.PIPE,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE,
#         universal_newlines=True,
#         preexec_fn=os.setsid) # to make a new process with a new PGID
#     pid = p.pid
#     pgid = os.getpgid(pid)
#     log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
#     stdout, stderr = p.communicate(cmd)
#     rc = p.returncode
#     err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
#         pid,
#         pgid,
#         rc,
#         stderr.strip(),
#         stdout.strip())
#     if rc:
#         # kill all child processes
#         try:
#             os.killpg(pgid, signal.SIGKILL)
#         except:
#             log.error(err_str)
# 
#     return (rc, stdout.strip('\n'), stderr.strip('\n'))
# 
# 
# def gzip_cmd(src, dest, decompress=True, rm=True):
#     """
#     Gzip Compress or Decompress files using gzip module in python
#     rm, True/False, whether remove old file
# 
#     # check the src file by extension: .gz
#     """
#     if os.path.exists(dest):
#         log.warning('file exists, skipped - {}'.format(dest))
#     else:
#         if decompress:
#             if is_gz(src):
#                 with gzip.open(src, 'rb') as r, open(dest, 'wb') as w:
#                     shutil.copyfileobj(r, w)
#             else:
#                 log.warning('not a gzipped file: {}'.format(src))
#                 shutil.copy(src, dest)
#         else:
#             if is_gz(src):
#                 log.warning('input is gzipped file, no need gzip')
#                 shutil.copy(src, dest)
#             else:
#                 with open(src, 'rb') as r, gzip.open(dest, 'wb', compresslevel=1) as w:
#                     shutil.copyfileobj(r, w)
# 
#     # output
#     if rm is True:
#         os.remove(src)
# 
#     return dest
# 
# 
# def list_uniquer(seq, sorted=True, idfun=None):
#     """
#     seq: a list with items
#     sorted: whether sort the output(unique)
# 
#     get the unique of inlist
# 
#     see1: Markus
#     remove duplicates from a list while perserving order
#     https://stackoverflow.com/a/480227/2530783
# 
#     def f7(seq):
#         seen = set()
#         seen_add = seen.add
#         return [x for x in seq if not (x in seen or seen_add(x))]
# 
#     see2: ctcherry
#     https://stackoverflow.com/a/89202/2530783
# 
#     def f5(seq, idfun=None):
#         # order preserving
#         if idfun is None:
#             def idfun(x): return x
#         seen = {}
#         result = []
#         for item in seq:
#             marker = idfun(item)
#             # in old Python versions:
#             # if seen.has_key(marker)
#             # but in new ones:
#             if marker in seen: continue
#             seen[marker] = 1
#             result.append(item)
#         return result
#     """
#     if idfun is None:
#         def idfun(x): return x # for None
# 
#     if not isinstance(seq, list):
#         log.error('list required, but get {}'.format(type(seq)))
#         return [] # blank
#     elif sorted is True:
#         return list(set(seq))
#     else:
#         seen = set()
#         return [x for x in seq if x not in seen and not seen.add(x)]
# 
# 
# def args_checker(d, x, update=False):
#     """Check if dict and x are consitent
#     d is dict
#     x is pickle file
#     """
#     assert isinstance(d, dict)
#     flag = None
#     if os.path.exists(x):
#         # read file to dict
#         with open(x, 'rb') as fh:
#             d_checker = pickle.load(fh)
#         if d == d_checker:
#             flag = True
#         else:
#             if update:
#                 with open(x, 'wb') as fo:
#                     pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
#     elif isinstance(x, str):
#         # save dict to new file
#         with open(x, 'wb') as fo:
#             pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
#     else:
#         log.error('illegal x= argument: %s' % x)
# 
#     return flag
# 
# 
# def args_logger(d, x, overwrite=False):
#     """Format dict, save to file
#         key: value
#     """
#     assert isinstance(d, dict)
#     n = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
#     if os.path.exists(x) and overwrite is False:
#         return True
#     else:
#         with open(x, 'wt') as fo:
#             fo.write('\n'.join(n) + '\n')
#         return '\n'.join(n)
# 
# 
# ## 3. utils ##
# def sam_flag_check(query, subject):
#     """
#     Check two numbers, query (for filtering) in subject or not
#     convert to binary mode
#     q: 0000010  (2)
#     s: 1011011  (91)
#     q in s
#     range: 0 - 2048 (SAM flag)
#     """
#     def to_bin(n):
#         return '{0:012b}'.format(n)
# 
#     # convert to binary mode
#     q = to_bin(eval(query))
#     s = to_bin(eval(subject))
# 
#     # check q, s
#     flag = True
#     for j, k in zip(q[::-1], s[::-1]):
#         if not j == '1':
#             continue
#         if eval(j) - eval(k) > 0:
#             flag = False
#             break
# 
#     return flag
# 
# 
# ################################################################################
# ## functions for pipeline
# class Toml(object):
#     """
#     Processing TOML files
# 
#     {toml, yaml, json, pickle} -> {dict} -> {toml, yaml, json, pickle}
# 
#     guess input format
# 
#     Known issue:
#     1. do not support Date and time 
#     2. do not support non{str, list} in to_text() 
#     3. 
# 
#     !!!! Do not support Date and time in TOML format
#     """
#     def __init__(self, x=None, **kwargs):
#         self = update_obj(self, kwargs, force=True)
#         self.x = x
# 
#         if x is None:
#             # self.x_dict = {}
#             pass
#         else:
#             x_fmt = self.guess_fmt(x)
#             if x_fmt:
#                 if x_fmt == 'dict':
#                     self.x_dict = collections.OrderedDict(sorted(x.items()))
#                 elif x_fmt == 'JSON':
#                     self.x_dict = self.from_json(x)
#                 elif x_fmt == 'YAML':
#                     self.x_dict = self.from_yam(x)
#                 elif x_fmt == 'TOML':
#                     self.x_dict = self.from_toml(x)
#                 elif x_fmt == 'pickle':
#                     self.x_dict = self.from_pickle(x)
#                 else:
#                     self.x_dict = {} # empty
#             else:
#                 raise ValueError('x, str or dict expected, got {}'.format(
#                     type(x).__name__))
# 
# 
#     def guess_fmt(self, x):
#         """
#         Guess the format of input x
#     
#         file format:
#         - yaml
#         - json
#         - pickle
# 
#         data format:
#         - dict
#         """
#         # supported
#         fmt = {
#             'json': 'JSON',
#             'yaml': 'YAML',
#             'yml': "YAML",
#             'toml': 'TOML',
#             'pickle': 'pickle'
#         }
# 
#         if isinstance(x, str):
#             x_ext = os.path.splitext(x)[1]
#             x_ext = x_ext.lstrip('.').lower()
#             x_fmt = fmt.get(x_ext, None)
# 
#         elif isinstance(x, dict):
#             x_fmt = 'dict'
# 
#         else:
#             x_fmt = None
# 
#         return x_fmt
# 
# 
#     def _tmp(self, suffix='.txt'):
#         """
#         Create a tmp file to save json object
#         """
#         tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
#             delete=False)
#         return tmp.name
# 
# 
#     def to_json(self, x):
#         """
#         Save dict to file in Json format
#         """
#         x = file_abspath(x) # convert to absolute path
#         if isinstance(self.x_dict, dict):
#             if isinstance(x, str):
#                 if file_exists(os.path.dirname(x), isfile=False):
#                     with open(x, 'wt') as w:
#                         json.dump(self.x_dict, w, indent=4, sort_keys=True)
# 
#                     return x
#                 else:
#                     log.warning('to_json() failed, could not write to file: \
#                         {}'.format(x))
#             else: 
#                 log.warning('to_json() failed, x, str expected, got {}'.format(
#                     type(x).__name__))
#         else:
#             log.warning('Expect input for Toml(x=)')
# 
# 
#     def from_json(self, x):
#         """
#         Parsing data from JSON file
#         """
#         try:
#             with open(x, 'r') as r:
#                 if os.path.getsize(x) > 0:
#                     d = json.load(r)
#                     return collections.OrderedDict(sorted(d.items()))
#         except:
#             log.warning('from_json() failed')
# 
# 
#     def to_dict(self):
#         """
#         Convert file to dict
#         """
#         if isinstance(self.x_dict, dict):
#             return self.x_dict
#         else:
#             log.warning('Expect input for Toml(x=)')
# 
# 
#     def from_dict(self, x):
#         """
#         Parsing data from dict
#         """
#         return collections.OrderedDict(sorted(x.items())) \
#             if isinstance(x, dict) else None
# 
# 
#     def to_yaml(self, x):
#         """
#         Saving data to file in YAML format
#         """
#         x = file_abspath(x) # convert to absolute path
#         if isinstance(self.x_dict, dict):
#             if isinstance(x, str):
#                 if file_exists(os.path.dirname(x), isfile=False):
#                     with open(x, 'wt') as w:
#                         json.dump(self.x_dict, w, indent=4, sort_keys=True)
# 
#                     return x
#                 else:
#                     log.warning('to_yaml() failed, could not write to file: \
#                         {}'.format(x))
#             else: 
#                 log.warning('to_yaml() failed, x, str expected, got {}'.format(
#                     type(x).__name__))
#         else:
#             log.warning('Expect input for Toml(x=)')
#         
# 
#     def from_yaml(self, x):
#         """
#         Parsing data from YAML file
#         """
#         with open(x, 'r') as r:
#             try:
#                 d = yaml.safe_load(r)
#                 return collections.OrderedDict(sorted(d.items()))
#             except yaml.YAMLError as exc:
#                 log.warning(exc)
# 
# 
#     def to_toml(self, x):
#         """
#         Saving config to TOML file
#         """
#         x = file_abspath(x) # convert to absolute path
#         if isinstance(self.x_dict, dict):
#             if isinstance(x, str):
#                 if file_exists(os.path.dirname(x), isfile=False):
#                     with open(x, 'w') as w:
#                         x_new = toml.dump(self.x_dict, w)
# 
#                     return x
#                 else:
#                     log.warning('to_toml() failed, could not write to file: \
#                         {}'.format(x))
#             else: 
#                 log.warning('to_toml() failed, x, str expected, got {}'.format(
#                     type(x).__name__))
#         else:
#             log.warning('Expect input for Toml(x=)')
# 
# 
#     def from_toml(self, x):
#         """
#         Parsing TOML file
#         """
#         try:
#             d = toml.load(x)
#             return collections.OrderedDict(sorted(d.items()))
#         except:
#             log.warning('failed to parsing file: {}'.format(x))
# 
# 
#     def from_pickle(self, x):
#         """
#         Parsing data from pickle file
#         """
#         try:
#             with open(x, 'rb') as r:
#                 d = pickle.load(r)
#                 return collections.OrderedDict(sorted(d.items()))
#         except:
#             log.warning('failed to parsing file: {}'.format(x))
# 
# 
#     def to_pickle(self, x):
#         """
#         Saving data to pickle file
#         """
#         x = file_abspath(x) # convert to absolute path
#         if isinstance(self.x_dict, dict):
#             if isinstance(x, str):
#                 if file_exists(os.path.dirname(x), isfile=False):
#                     with open(x, 'wb') as w:
#                         pickle.dump(self.x_dict, w, 
#                                     protocol=pickle.HIGHEST_PROTOCOL)
# 
#                     return x
#                 else:
#                     log.warning('to_pickle() failed, could not write to file: \
#                         {}'.format(x))
#             else: 
#                 log.warning('to_pickle() failed, x, str expected, got {}'.format(
#                     type(x).__name__))
#         else:
#             log.warning('Expect input for Toml(x=)')
# 
# 
#     def to_text(self, x):
#         """
#         Saving dict in text format
#         """
#         x = file_abspath(x) # convert to absolute path
#         if isinstance(self.x_dict, dict):
#             if isinstance(x, str):
#                 if file_exists(os.path.dirname(x), isfile=False):
#                     msg = '\n'.join([
#                         '{:30s} | {:<40s}'.format(k, str(v)) for k, v in self.x_dict.items()
#                         ])
#                     with open(x, 'wt') as w:
#                         w.write(msg + '\n')
# 
#                     return x
#                 else:
#                     log.warning('to_text() failed, could not write to file: \
#                         {}'.format(x))
#             else: 
#                 log.warning('to_text() failed, x, str expected, got {}'.format(
#                     type(x).__name__))
#         else:
#             log.warning('Expect input for Toml(x=)')
# 
# 
# class Config(object):
#     """Working with config, in dict/yaml/toml/pickle formats
#     load/dump
#     
#     Example:
#     1. write to file
#     >>> Config(d).dump('out.json')
#     >>> Config(d).dump('out.toml')
#     >>> Config(d).dump('out.pickle')
#     
#     2. load from file
#     >>> d = Config().load('in.yaml')
# 
#     read/write data
#     """
#     def __init__(self, x=None, **kwargs):
#         self = update_obj(self, kwargs, force=True)
#         self.x = x
# 
# 
#     def load(self, x=None):
#         """Read data from x, auto-recognize the file-type
#         toml
#         json
#         pickle
#         txt
#         ...
#         """
#         if x == None:
#             x = self.x # dict or str
# 
#         if x is None:
#             x_dict = None # {} ?
#         elif isinstance(x, dict):
#             x_dict = collections.OrderedDict(sorted(x.items()))
#         elif isinstance(x, str):
#             reader = self.get_reader(x)
#             if reader is None:
#                 x_dict = None
#                 log.error('unknown x, {}'.format(x))
#             else:
#                 x_dict = reader(x)
#         else:
#             x_dict = None
#             log.warning('dump(x=) dict,str expect, got {}'.format(
#                 type(x).__name__))
# 
#         return x_dict
# 
# 
#     def dump(self, d=None, x=None):
#         """Write data to file x, auto-recognize the file-type
#         d str or dict, data
#         x str file to save data(dict)
# 
#         toml
#         json
#         pickle
#         txt
#         ...
#         """
#         if d is None:
#             d = self.load(self.x)
#         # make sure: dict
#         if isinstance(x, str):
#             writer = self.get_writer(x)
#             if writer is None:
#                 log.error('unknown x, {}'.format(x))
#             else:
#                 writer(d, x)
#         else:
#             log.warning('dump(x=) expect str, got {}'.format(
#                 type(x).__name__))
# 
# 
#     def guess_format(self, x):
#         """Guess the file format, by file extension
#     
#         file format:
#         - toml
#         - yaml
#         - json
#         - pickle
# 
#         data format:
#         - dict
#         """
#         formats = {
#             'json': 'json',
#             'yaml': 'yaml',
#             'yml': "yaml",
#             'toml': 'toml',
#             'pickle': 'pickle'
#         }
# 
#         if isinstance(x, str):
#             x_ext = os.path.splitext(x)[1]
#             x_ext = x_ext.lstrip('.').lower()
#             x_format = formats.get(x_ext, None)
#         elif isinstance(x, dict):
#             x_format = 'dict'
#         else:
#             x_format = None
# 
#         return x_format
# 
# 
#     def get_reader(self, x):
#         """Get the reader for file x, based on the file extension
#     
#         could be: json/yaml/toml/pickle
#         """
#         x_format = self.guess_format(x)
#         readers = {
#             'json': self.from_json,
#             'yaml': self.from_yaml,
#             'toml': self.from_toml,
#             'pickle': self.from_pickle
#         }
#         return readers.get(x_format, None)
# 
# 
#     def get_writer(self, x):
#         """Get the reader for file x, based on the file extension
# 
#         could be: json/yaml/toml/pickle
#         """
#         x_format = self.guess_format(x)
#         writers = {
#             'json': self.to_json,
#             'yaml': self.to_yaml,
#             'toml': self.to_toml,
#             'pickle': self.to_pickle
#         }
#         return writers.get(x_format, None)
# 
# 
#     def from_json(self, x):
#         """Loding data from JSON file
#         x should be file
#         """
#         d = None
#         if file_exists(x):
#             try:
#                 with open(x, 'r') as r:
#                     if os.path.getsize(x) > 0:
#                         d = json.load(r)
#                         d = collections.OrderedDict(sorted(d.items()))
#             except Exception as exc:
#                 log.error('from_json() failed, {}'.format(exc))
#             finally:
#                 return d
#         else:
#             log.error('from_json() failed, file not exists: {}'.format(x))
# 
# 
#     def from_yaml(self, x):
#         """Loding data from YAML file
#         x should be file
#         """
#         d = None
#         if file_exists(x):
#             try:
#                 with open(x, 'r') as r:
#                     if os.path.getsize(x) > 0:
#                         d = yaml.load(r, Loader=yaml.FullLoader)
#                         d = collections.OrderedDict(sorted(d.items()))
#             except Exception as exc:
#                 log.error('from_yaml() failed, {}'.format(exc))
#             finally:
#                 return d
#         else:
#             log.error('from_yaml() failed, file not exists: {}'.format(x))
#         # with open(x, 'r') as r:
#         #     try:
#         #         d = yaml.safe_load(r)
#         #         return collections.OrderedDict(sorted(d.items()))
#         #     except yaml.YAMLError as exc:
#         #         log.warning(exc)
# 
# 
#     def from_toml(self, x):
#         """Loding data from TOML file
#         x should be file
#         """
#         d = None
#         if file_exists(x):
#             try:
#                 with open(x, 'r') as r:
#                     if os.path.getsize(x) > 0:
#                         d = toml.load(x)
#                         d = collections.OrderedDict(sorted(d.items()))
#             except Exception as exc:
#                 log.error('from_toml() failed, {}'.format(exc))
#             finally:
#                 return d
#         else:
#             log.error('from_toml() failed, file not exists: {}'.format(x))
# 
# 
#     def from_pickle(self, x):
#         """Loding data from pickle file
#         x should be file
#         """
#         d = None
#         if file_exists(x):
#             try:
#                 with open(x, 'rb') as r:
#                     if os.path.getsize(x) > 0:
#                         d = pickle.load(r)
#                         d = collections.OrderedDict(sorted(d.items()))
#             except Exception as exc:
#                 log.error('from_pickle() failed, {}'.format(exc))
#             finally:
#                 return d
#         else:
#             log.error('from_pickle() failed, file not exists: {}'.format(x))
# 
# 
#     def to_json(self, d, x):
#         """Writing data to JSON file
#         d dict, data to file
#         x None or str, path to JSON file, or return string
#         """
#         x = file_abspath(x)
#         if not isinstance(d, dict):
#             log.error('to_json(d=) failed, dict expect, got {}'.format(
#                 type(d).__name__))
#         elif not isinstance(x, str):
#             log.error('to_json(d=) failed, str expect, got {}'.format(
#                 type(x).__name__))
#         elif not file_exists(os.path.dirname(x), isfile=False):
#             log.error('to_json(x=) failed, file not exists: {}'.format(x))
#         else:
#             try:
#                 with open(x, 'wt') as w:
#                     json.dump(d, w, indent=4, sort_keys=True)
#                 # return x
#             except Exception as exc:
#                 log.error('to_json() failed, {}'.format(exc))
# 
# 
#     def to_yaml(self, d, x):
#         """Writing data to YAML file
#         d dict, data to file
#         x str, path to YAML file
# 
#         yaml.dump(), does not support OrderedDict
#         Solution: OrderedDict -> json -> dict
#         """
#         x_yaml = x
#         x = os.path.splitext(x_yaml)[0] + '.toml'
#         log.warning('OrderedDict is not supported in YAML, save as TOML instead: {}'.format(x))
#         # check
#         x = file_abspath(x)
#         if not isinstance(d, dict):
#             log.error('to_yaml(d=) failed, dict expect, got {}'.format(
#                 type(d).__name__))
#         elif not isinstance(x, str):
#             log.error('to_yaml(d=) failed, str expect, got {}'.format(
#                 type(x).__name__))
#         elif not file_exists(os.path.dirname(x), isfile=False):
#             log.error('to_yaml(x=) failed, file not exists: {}'.format(x))
#         else:
#             try:
#                 with open(x, 'wt') as w:
#                     toml.dump(d, w)
#                 # return x
#             except Exception as exc:
#                 log.error('to_yaml() failed, {}'.format(exc))
# 
# 
#     def to_toml(self, d, x):
#         """Writing data to TOML file
#         d dict, data to file
#         x str, path to TOML file
#         """        
#         x = file_abspath(x)
#         if not isinstance(d, dict):
#             log.error('to_toml(d=) failed, dict expect, got {}'.format(
#                 type(d).__name__))
#         elif not isinstance(x, str):
#             log.error('to_toml(d=) failed, str expect, got {}'.format(
#                 type(x).__name__))
#         elif not file_exists(os.path.dirname(x), isfile=False):
#             log.error('to_toml(d=) failed, file not exists: {}'.format(x))
#         else:
#             try:
#                 with open(x, 'wt') as w:
#                     toml.dump(d, w)
#                 # return x
#             except Exception as exc:
#                 log.error('to_toml() failed, {}'.format(exc))
# 
# 
#     def to_pickle(self, d, x):
#         """Writing data to pickle file
#         d dict, data to file
#         x str, path to pickle file
#         """        
#         x = file_abspath(x)
#         if not isinstance(d, dict):
#             log.error('to_pickle(d=) failed, dict expect, got {}'.format(
#                 type(d).__name__))
#         elif not isinstance(x, str):
#             log.error('to_pickle(x=) failed, str expect, got {}'.format(
#                 type(x).__name__))
#         elif not file_exists(os.path.dirname(x), isfile=False):
#             log.error('to_pickle(x=) failed, file not exists: {}'.format(x))
#         else:
#             try:
#                 with open(x, 'wb') as w:
#                     pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)
#                 # return x
#             except Exception as exc:
#                 log.error('to_pickle() failed, {}'.format(exc))
# 
# 
#     def to_log(self, d, x, stdout=False):
#         """Writing data to log file: key: value format
#         d dict, data to file
#         x str, path to pickle file
#         """
#         x = file_abspath(x)
#         if not isinstance(d, dict):
#             log.error('to_log(d=) failed, dict expect, got {}'.format(
#                 type(d).__name__))
#         elif not isinstance(x, str):
#             log.error('to_log(x=) failed, str expect, got {}'.format(
#                 type(x).__name__))
#         elif not file_exists(os.path.dirname(x), isfile=False):
#             log.error('to_log(x=) failed, file not exists: {}'.format(x))
#         else:
#             try:
#                 # organize msg
#                 msg = []
#                 for k, v in d.items():
#                     if isinstance(v, str) or isinstance(v, numbers.Number) or isinstance(v, bool):
#                         v = str(v)
#                     elif isinstance(v, list):
#                         v = ', '.join(map(str, v))
#                     else:
#                         v = '...' # skip
#                     msg.append('{:30s} | {:<40s}'.format(k, v))
#                 # save
#                 with open(x, 'wt') as w:
#                     w.write('\n'.join(msg) + '\n')
#                 if stdout:
#                     print('\n'.join(msg))
#                 # return x
#             except Exception as exc:
#                 log.error('to_log() failed, {}'.format(exc))
# 
# 
#     def _tmp(self, suffix='.txt'):
#         """
#         Create a tmp file to save json object
#         """
#         tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
#             delete=False)
#         return tmp.name
#     
#     
# 
# def in_dict(d, k):
#     """
#     Check the keys in dict or not
#     """
#     assert isinstance(d, dict)
#     if isinstance(k, str):
#         k_list = [k]
#     elif isinstance(k, list):
#         k_list = list(map(str, k))
#     else:
#         log.warning('expect str and list, not {}'.format(type(k)))
#         return False
# 
#     return all([i in d for i in k_list])
# 
# 
# def in_attr(x, a, return_values=True):
#     """
#     Check a (attributes) in object a or not
#     return the values or not
#     """
#     if isinstance(a, str):
#         a_list = [a]
#     elif isinstance(a, list):
#         a_list = list(map(str, a))
#     else:
#         log.warning('expect str and list, not {}'.format(type(a)))
#         return False
# 
#     # status
#     status = all([hasattr(x, i) for i in a_list])
# 
#     if status and return_values:
#         # values
#         return [getattr(x, i) for i in a_list]
#     else:
#         return status
# 
# 
# def dict_to_log(d, x, overwrite=False):
#     """
#     Convert dict to log style
#         key | value
#     """
#     assert isinstance(d, dict)
#     logout = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
#     if overwrite is True or not os.path.exists(x):
#         with open(x, 'wt') as w:
#             w.write('\n'.join(logout) + '\n')
# 
#     return '\n'.join(logout)
# 
# 
# def dict_to_pickle(d, x):
#     """
#     Convert dict to pickle
#     """
#     assert isinstance(d, dict)
#     with open(x, 'wb') as w:
#         pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)
# 
# 
# def pickle_to_dict(x):
#     """
#     Convert pickle file to dict
#     """
#     with open(x, 'rb') as r:
#         return pickle.load(r)
# 
# 
# class Dict2Obj(object):
#     """
#     >>> d = {'a': 1, 'b': 2}
#     >>> b = Dict2obj(**d)
#     >>> print(b.a)
#     1
# 
#     class AttributeDict(object):
#     A class to convert a nested Dictionary into an object with key-values
#     accessibly using attribute notation (AttributeDict.attribute) instead of
#     key notation (Dict["key"]). This class recursively sets Dicts to objects,
#     allowing you to recurse down nested dicts (like: AttributeDict.attr.attr)
#     from: http://databio.org/posts/python_AttributeDict.html
#     by Nathan Sheffield
#     """
#     def __init__(self, **entries):
#         self.add_entries(**entries)
# 
#     def add_entries(self, **entries):
#         for key, value in entries.items():
#             if type(value) is dict:
#                 self.__dict__[key] = AttributeDict(**value)
#             else:
#                 self.__dict__[key] = value
# 
#     def __getitem__(self, key):
#         """
#         Provides dict-style access to attributes
#         """
#         return getattr(self, key)
# 
# 
# class Dict2Class(object):
#     """
#     Turns a dictionary into a class
#     from: https://www.blog.pythonlibrary.org/2014/02/14/python-101-how-to-change-a-dict-into-a-class/
#     by Mike
#     """
#     #----------------------------------------------------------------------
#     def __init__(self, d):
#         """Constructor"""
#         for k, v in d.items():
#             setattr(self, k, v)
# 
#         # for key in d:
#         #     setattr(self, key, d[key])
# 
# 
# def bed2gtf(infile, outfile):
#     """Convert BED to GTF
#     chrom chromStart chromEnd name score strand
#     """
#     with open(infile) as r, open(outfile, 'wt') as w:
#         for line in r:
#             fields = line.strip().split('\t')
#             start = int(fields[1]) + 1
#             w.write('\t'.join([
#                 fields[0],
#                 'BED_file',
#                 'gene',
#                 str(start),
#                 fields[2],
#                 '.',
#                 fields[5],
#                 '.',
#                 'gene_id "{}"; gene_name "{}"'.format(fields[3], fields[3])
#                 ]) + '\n')
#     return outfile
# 
# 
# ### deprecated - BEGIN ###
# def listfiles(path, full_name=True, recursive=False, include_dir=False):
#     """
#     List all the files within the path
#     """
#     out = []
#     for root, dirs, files in os.walk(path):
#         if full_name:
#             dirs = [os.path.join(root, d) for d in dirs]
#             files = [os.path.join(root, f) for f in files]
#         out += files
# 
#         if include_dir:
#             out += dirs
# 
#         if recursive is False:
#             break
#     return out
# 
# 
# def listfiles2(pattern, path='.', full_name=True, recursive=False):
#     """
#     List all the files in specific directory
#     fnmatch.fnmatch()
# 
#     pattern:
# 
#     *       matches everything
#     ?       matches any single character
#     [seq]   matches any character in seq
#     [!seq]  matches any char not in seq
# 
#     An initial period in FILENAME is not special.
#     Both FILENAME and PATTERN are first case-normalized
#     if the operating system requires it.
#     If you don't want this, use fnmatchcase(FILENAME, PATTERN).
# 
#     example:
#     listfiles('*.fq', './')
#     """
#     fn_list = listfiles(path, full_name, recursive, include_dir=False)
#     fn_list = [f for f in fn_list if fnmatch.fnmatch(f, pattern)]
#     return fn_list
# 
# 
# def is_path(path, create = True):
#     """
#     Check path, whether a directory or not
#     if not, create it
#     """
#     assert isinstance(path, str)
#     if os.path.exists(path):
#         return True
#     else:
#         if create:
#             try:
#                 os.makedirs(path)
#                 return True
#             except IOError:
#                 log.error('failed to create directories: %s' % path)
#         else:
#             return False
# ### deprecated - END ###
# 
# 
# class Json(object):
#     """
#     Wrapper for *.json file
#     1.json to dict
#     2.dict to json -> file
#     """
#     def __init__(self, input, **kwargs):
#         self.input = input
#         self.mission()
# 
# 
#     def mission(self):
#         """
#         Check what to do, based on the input args
#         """
#         # json_file to dict
#         if input is None:
#             log.warning('require, dict or str (file), Nonetype detected')
#         elif isinstance(self.input, str) and os.path.exists(self.input):
#             self.json_file_in = self.input
#             # self.reader() # updated
#         elif isinstance(self.input, dict):
#             self.d = self.input
#             # self.writer()
#         else:
#             log.warning('expect dict or str, get {}'.format(type(self.input).__name__))
# 
# 
#     def _tmp(self):
#         """
#         Create a tmp file to save json object
#         """
#         tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json',
#             delete=False)
#         return tmp.name
# 
# 
#     def to_dict(self):
#         # self.json_file_in = self.input
#         return self.reader()
# 
# 
#     def to_json(self, file=None):
#         return self.writer(file)
# 
# 
#     def reader(self):
#         """
#         Read json file as dict
#         """
#         if isinstance(self.input, dict):
#             self.d = self.input
#         else:
#             try:
#                 with open(self.json_file_in) as r:
#                     if os.path.getsize(self.json_file_in) > 0:
#                         self.d = json.load(r)
#                     else:
#                         self.d = {}
#             except:
#                 log.warning('failed reading json file: {}'.format(self.json_file_in))
#                 self.d = {}
# 
#         return self.d
# 
# 
#     def writer(self, json_file_out=None):
#         """
#         Write d (dict) to file x, in json format
#         """
#         # save to file
#         if json_file_out is None:
#             json_file_out = self._tmp()
# 
#         try:
#             with open(json_file_out, 'wt') as w:
#                 json.dump(self.d, w, indent=4, sort_keys=True)
#         except:
#             log.warning('failed saving file: {}'.format(json_file_out))
# 
#         return json_file_out
# 
# 
# class Genome(object):
#     """
#     List related information of specific genome
#     1. get_fa(), genome fasta
#     2. get_fasize(), genome fasta size
#     3. bowtie_index(), bowtie index, optional, rRNA=True
#     4. bowtie2_index(), bowtie2 index, optional, rRNA=True
#     5. star_index(), STAR index, optional, rRNA=True
#     6. gene_bed(),
#     7. gene_rmsk(),
#     8. gene_gtf(), optional, version='ucsc|ensembl|ncbi'
#     9. te_gtf(), optional, version='ucsc'
#     10. te_consensus(), optional, fruitfly()
#     ...
# 
#     directory structure of genome should be like this:
#     /path-to-data/{genome}/
#         |- bigZips  # genome fasta, fasize, chromosome
#         |- annotation_and_repeats  # gtf, bed, rRNA, tRNA, annotation
#         |- bowtie_index
#         |- bowtie2_index
#         |- STAR_index
#         |- hisat2_index
#         |- phylop100
#         |- ...
# 
#     default: $HOME/data/genome/{genome}
# 
#     """
#     def __init__(self, genome, genome_path=None,
#         repeat_masked_genome=False, **kwargs):
#         assert isinstance(genome, str)
#         self.genome = genome
#         self.repeat_masked_genome = repeat_masked_genome
#         self.kwargs = kwargs
# 
#         # path
#         if genome_path is None:
#             genome_path = os.path.join(str(pathlib.Path.home()),
#                 'data', 'genome')
#         self.genome_path = genome_path
# 
#         # # deprecated
#         # if not supportedGenome(genome):
#         #     log.error('genome not supported: {}'.foramt(genome))
# 
# 
#     def get_fa(self):
#         """
#         Get the fasta file of specific genome
#         {genome}/bigZips/{genome}.fa
#         also check ".gz" file
#         """
#         fa = os.path.join(self.genome_path, self.genome, 'bigZips',
#             self.genome + '.fa')
#         if not os.path.exists(fa):
#             # gencode version
#             fa = os.path.join(self.genome_path, self.genome, 'fasta',
#                 self.genome + '.fa')
# 
#         fa_gz = fa + '.gz'
#         if not os.path.exists(fa):
#             if os.path.exists(fa_gz):
#                 log.error('require to unzip the fasta file: %s' % fa_gz)
#             else:
#                 log.error('fasta file not detected: %s' % fa)
#             return None
#         else:
#             return fa
# 
# 
#     def get_fasize(self):
#         """Get the fasta size file, chromosome size
#         optional, fetch chrom size from ucsc
#         http://hgdownload.cse.ucsc.edu/goldenPath/<db>/bigZips/<db>.chrom.sizes
# 
#         or using UCSC tool: fetchChromSizes
#         fetchChromSizes hg39 > hg38.chrom.sizes
#         """
#         fa = self.get_fa()
#         fa_size = fa + '.chrom.sizes'
# 
#         if not os.path.exists(fa_size):
#             # log.info('Downloading chrom.sizes from UCSC')
#             # url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes' % self.genome
#             # download(url, fa_size)
#             log.warning('file not exists, run samtools faidx to generate it')
#             pysam.faidx(fa) # create *.fa.fai
#             os.rename(fa + '.fai', fa_size)
# 
#         return fa_size
# 
# 
#     def phylop100(self):
#         """Return the phylop100 bigWig file of hg19, only
#         for conservation analysis
#         """
#         p = os.path.join(self.genome_path, self.genome, 'phyloP100way',
#             self.genome + '.100way.phyloP100way.bw')
#         if not os.path.exists(p):
#             p = None
#         return p
# 
# 
#     def gene_bed(self, version='refseq', rmsk=False):
#         """Return the gene annotation in BED format
#         support UCSC, ensembl, gencode
#         """
#         if rmsk:
#             suffix = '.rmsk.bed'
#         else:
#             suffix = '.refseq.bed'
#         g = os.path.join(self.genome_path, self.genome,
#             'annotation_and_repeats', self.genome + suffix)
#         if not os.path.exists(g):
#             g = None
#         return g
# 
# 
#     def gene_gtf(self, version='refseq'):
#         """Return the gene annotation in GTF format
#         support refseq, ensembl, gencode
#         """
#         version = version.lower() #
# 
#         gtf = os.path.join(
#             self.genome_path,
#             self.genome,
#             'annotation_and_repeats',
#             self.genome + '.' + version + '.gtf')
# 
#         if not os.path.exists(gtf):
#             gtf = os.path.join(
#             self.genome_path,
#             self.genome,
#             'gtf',
#             self.genome + '.' + version + '.gtf')
# 
#         if not os.path.exists(gtf):
#             gtf = None
# 
#         return gtf
# 
# 
#     def te(self, format='gtf'):
#         """Return TE annotation of the genome
#         or return TE consensus sequence for the genome (dm3)
#         """
#         # only dm3 supported
#         te_gtf = os.path.join(self.genome_path, self.genome,
#             self.genome + '_transposon',
#             self.genome + '_transposon.' + format)
#         if not os.path.exists(te_gtf):
#             te_gtf = None
# 
#         return te_gtf
# 
# 
#     def piRNA_cluster(self, format='gtf'):
#         """Return TE annotation of the genome
#         or return TE consensus sequence for the genome (dm3)
#         """
#         # only dm3 supported
#         te_gtf = os.path.join(self.genome_path, self.genome,
#             self.genome + '_piRNA_clusters',
#             self.genome + '_piRNA_clusters.' + format)
#         if not os.path.exists(te_gtf):
#             te_gtf = None
# 
#         return te_gtf
# 
# 
# class Bam(object):
#     """
#     Manipulate BAM files
#     - sort
#     - index
#     - merge
#     - count
#     - to_bed
#     - rmdup
#     - ...
# 
#     Using Pysam, Pybedtools, ...
# 
#     code from cgat:
#     """
#     def __init__(self, infile, threads=4):
#         self.bam = infile
#         self.threads = threads
#         # self.bed = self.to_bed()
# 
# 
#     def index(self):
#         """
#         Create index for bam
#         """
#         bai = self.bam + '.bai'
#         if not os.path.exists(bai):
#             pysam.index(self.bam)
# 
#         return os.path.exists(bai)
# 
# 
#     def sort(self, outfile=None, by_name=False, overwrite=False):
#         """
#         Sort bam file by position (default)
#         save to *.sorted.bam (or specify the name)
#         """
#         if outfile is None:
#             outfile = os.path.splitext(self.bam)[0] + '.sorted.bam'
# 
#         if os.path.exists(outfile) and overwrite is False:
#             logging.info('file exists: {}'.format(outfile))
#         else:
#             tmp = pysam.sort('-@', str(self.threads), '-o', outfile, self.bam)
# 
#         return outfile
# 
# 
#     def merge(self):
#         """
#         Merge multiple BAM files using samtools
#         """
# 
#         # pysam.merge('')
#         pass
# 
# 
#     def count(self):
#         """Using samtools view -c"""
#         x = pysam.view('-c', self.bam)
#         return int(x.strip())
# 
# 
#     def to_bed(self, outfile=None):
#         """Convert BAM to BED
#         pybetools
#         """
#         if outfile is None:
#             outfile = os.path.splitext(self.bam)[0] + '.bed'
# 
#         if not os.path.exists(outfile):
#             pybedtools.BedTool(self.bam).bam_to_bed().saveas(outfile)
# 
#         return outfile
# 
# 
#     def rmdup(self, outfile=None, overwrite=False):
#         """
#         Remove duplicates using picard/sambamba
#         sambamba markdup -r --overflow-list-size 800000 raw.bam rmdup.bam
#         picard MarkDuplicates  REMOVE_SEQUENCING_DUPLICATES=True I=in.bam O=outfile.bam
#         """
#         if outfile is None:
#             outfile = os.path.splitext(self.bam)[0] + '.rmdup.bam'
# 
#         sambamba = shutil.which('sambamba')
#         log = outfile + '.sambamba.log'
#         cmd = '{} markdup -r -t {} --overflow-list-size 800000 \
#             --tmpdir="./" {} {} 2> {}'.format(
#             sambamba,
#             str(self.threads),
#             self.bam,
#             outfile,
#             log)
# 
#         if os.path.exists(outfile) and overwrite is False:
#             logging.info('file exists: {}'.format(outfile))
#         else:
#             run_shell_cmd(cmd)
# 
#         return outfile
# 
# 
#     def proper_pair(self, outfile=None, overwrite=False):
#         """
#         Extract proper pair
#         samtools view -f 2
#         """
#         if outfile is None:
#             outfile = os.path.splitext(self.bam)[0] + '.proper_pair.bam'
# 
#         if os.path.exists(outfile) and overwrite is False:
#             logging.info('file exists: {}'.format(outfile))
#         else:
#             pysam.view('-f', '2', '-h', '-b', '-@', str(self.threads),
#                 '-o', outfile, self.bam, catch_stdout=False)
# 
#         return outfile
# 
# 
#     def subset(self, size=20000, subdir=None):
#         """
#         Extract N reads from bam list
#         """
#         if subdir is None:
#             subdir = self._tmp(delete=False)
#         check_path(subdir)
# 
#         # src, dest
#         dest = os.path.join(subdir, os.path.basename(self.bam))
#         
#         # run
#         self.index()        
#         srcfile = pysam.AlignmentFile(self.bam, 'rb')
#         destfile = pysam.AlignmentFile(dest, 'wb', template=srcfile)
#         # counter
#         i = 0
#         for read in srcfile.fetch():
#             i +=1
#             if i > size:
#                 break
#             destfile.write(read)
# 
#         return dest
# 
# 
#     def _tmp(self, delete=True):
#         """
#         Create a tmp filename
#         """
#         tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=delete)
#         return tmp.name
# 
# 
#     ##########################################
#     ## code from cgat: BEGIN
#     ##########################################
#     def isPaired(self, topn=1000):
#         """
#         Check if infile contains paired end reads
# 
#         go through the topn alignments in file,
# 
#         return: True, any of the alignments are paired
#         """
#         samfile = pysam.AlignmentFile(self.bam)
#         n = 0
#         for read in samfile:
#             if read.is_paired:
#                 break
#             n += 1
#             if n == topn:
#                 break
# 
#         samfile.close()
# 
#         return n != topn
# 
# 
#     def getNumReads(self):
#         """
#         Count number of reads in bam file.
# 
#         This methods works through pysam.idxstats.
# 
#         Arguments
#         ---------
#         bamfile : string
#             Filename of :term:`bam` formatted file. The file needs
#             to be indexed.
#         Returns
#         -------
#         nreads : int
#             Number of reads
#         """
# 
#         lines = pysam.idxstats(self.bam).splitlines()
# 
#         try:
#             nreads = sum(
#                 map(int, [x.split("\t")[2]
#                           for x in lines if not x.startswith("#")]))
# 
#         except IndexError as msg:
#             raise IndexError(
#                 "can't get number of reads from bamfile, msg=%s, data=%s" %
#                 (msg, lines))
#         return nreads
# 
# 
#     def estimateInsertSizeDistribution(self, topn=10000, n=10,
#         method="picard", similarity_threshold=1.0, max_chunks=1000):
#         """
#         from pysam
#         Estimate insert size from a subset of alignments in a bam file.
# 
#         Several methods are implemented.
# 
#         picard
#             The method works analogous to picard by restricting the estimates
#             to a core distribution. The core distribution is defined as all
#             values that lie within n-times the median absolute deviation of
#             the full data set.
#         convergence
#             The method works similar to ``picard``, but continues reading
#             `alignments` until the mean and standard deviation stabilize.
#             The values returned are the median mean and median standard
#             deviation encountered.
# 
#         The method `convergence` is suited to RNA-seq data, as insert sizes
#         fluctuate siginificantly depending on the current region
#         being looked at.
# 
#         Only mapped and proper pairs are considered in the computation.
# 
#         Returns
#         -------
#         mean : float
#            Mean of insert sizes.
#         stddev : float
#            Standard deviation of insert sizes.
#         npairs : int
#            Number of read pairs used for the estimation
#         method : string
#            Estimation method
#         similarity_threshold : float
#            Similarity threshold to apply.
#         max_chunks : int
#            Maximum number of chunks of size `alignments` to be used
#            in the convergence method.
# 
#         """
# 
#         assert self.isPaired(self.bam), \
#             'can only estimate insert size from' \
#             'paired bam files'
# 
#         samfile = pysam.AlignmentFile(self.bam)
# 
#         def get_core_distribution(inserts, n):
#             # compute median absolute deviation
#             raw_median = numpy.median(inserts)
#             raw_median_dev = numpy.median(numpy.absolute(inserts - raw_median))
# 
#             # set thresholds
#             threshold_min = max(0, raw_median - n * raw_median_dev)
#             threshold_max = raw_median + n * raw_median_dev
# 
#             # define core distribution
#             return inserts[numpy.logical_and(inserts >= threshold_min,
#                                              inserts <= threshold_max)]
# 
#         if method == "picard":
# 
#             # only get first read in pair to avoid double counting
#             inserts = numpy.array(
#                 [read.template_length for read in samfile.head(n=topn)
#                  if read.is_proper_pair
#                  and not read.is_unmapped
#                  and not read.mate_is_unmapped
#                  and not read.is_read1
#                  and not read.is_duplicate
#                  and read.template_length > 0])
#             core = get_core_distribution(inserts, n)
# 
#             return numpy.mean(core), numpy.std(core), len(inserts)
# 
#         elif method == "convergence":
# 
#             means, stds, counts = [], [], []
#             last_mean = 0
#             iteration = 0
#             while iteration < max_chunks:
# 
#                 inserts = numpy.array(
#                     [read.template_length for read in samfile.head(
#                         n=topn,
#                         multiple_iterators=False)
#                      if read.is_proper_pair
#                      and not read.is_unmapped
#                      and not read.mate_is_unmapped
#                      and not read.is_read1
#                      and not read.is_duplicate
#                      and read.template_length > 0])
#                 core = get_core_distribution(inserts, n)
#                 means.append(numpy.mean(core))
#                 stds.append(numpy.std(core))
#                 counts.append(len(inserts))
#                 mean_core = get_core_distribution(numpy.array(means), 2)
#                 mm = numpy.mean(mean_core)
#                 if abs(mm - last_mean) < similarity_threshold:
#                     break
#                 last_mean = mm
# 
#             return numpy.median(means), numpy.median(stds), sum(counts)
#         else:
#             raise ValueError("unknown method '%s'" % method)
# 
# 
#     def estimateTagSize(self, topn=10, multiple="error"):
#         """
#         Estimate tag/read size from first alignments in file.
# 
#         Arguments
#         ---------
#         bamfile : string
#            Filename of :term:`bam` formatted file
#         alignments : int
#            Number of alignments to inspect
#         multiple : string
#            How to deal if there are multiple tag sizes present.
#            ``error`` will raise a warning, ``mean`` will return the
#            mean of the read lengths found. ``uniq`` will return a
#            unique list of read sizes found. ``all`` will return all
#            read sizes encountered.
# 
#         Returns
#         -------
#         size : int
#            The read size (actual, mean or list of read sizes)
# 
#         Raises
#         ------
#         ValueError
#            If there are multiple tag sizes present and `multiple` is set to
#            `error`.
# 
#         """
#         samfile = pysam.AlignmentFile(self.bam)
#         sizes = [read.rlen for read in samfile.head(topn)]
#         mi, ma = min(sizes), max(sizes)
# 
#         if mi == 0 and ma == 0:
#             sizes = [read.inferred_length for read in samfile.head(alignments)]
#             # remove 0 sizes (unaligned reads?)
#             sizes = [x for x in sizes if x > 0]
#             mi, ma = min(sizes), max(sizes)
# 
#         if mi != ma:
#             if multiple == "error":
#                 raise ValueError('multiple tag sizes in %s: %s' % (bamfile, sizes))
#             elif multiple == "mean":
#                 mi = int(sum(sizes) / len(sizes))
#             elif multiple == "uniq":
#                 mi = list(sorted(set(sizes)))
#             elif multiple == "all":
#                 return sizes
# 
#         return mi
# 
# 
#     def getNumberOfAlignments(self):
#         """
#         Return number of alignments in bamfile.
#         """        
#         if not os.path.exists(self.bam + '.bai'):
#             pysam.index(self.bam)
#         samfile = pysam.AlignmentFile(self.bam)
#         return samfile.mapped
#     ##########################################
#     ## code from cgat: END
#     ##########################################
# 
# ################################################################################
# ## TEMP
# def featureCounts_reader(x, bam_names=False):
#     """
#     Read fc .txt as data.frame
#     """
#     df = pd.read_csv(x, '\t', comment='#')
#     if bam_names:
#         bam_list = df.columns.to_list()[6:] # from 7-column to end
#         bam_list = [os.path.splitext(i)[0] for i in bam_list]
#         return [os.path.basename(i) for i in bam_list]
#     else:
#         return df
# 
# 
# ## design, create combinations
# def design_combinations(seq, n=2, return_index=True):
#     """
#     seq: the group name for each sample
# 
#     example:
#     seq = ['ctl', 'ctl', 'exp', 'exp'] # keep order
# 
#     : return_index=True
#     [[0, 1], [2, 3]]
# 
#     : return_index=False
#     ['ctl', 'exp']
#     """
#     seq_unique = list_uniquer(seq, sorted=False) # keep order
# 
#     if len(seq_unique) >= n:
#         item_pairs = list(combinations(seq_unique, n))
#         # for index
#         index_pairs = []
#         for (a, b) in item_pairs:
#             index_a = [i for i, x in enumerate(seq) if x == a]
#             index_b = [i for i, x in enumerate(seq) if x == b]
#             index_pairs.append([index_a, index_b])
#         return index_pairs if return_index else item_pairs
#     else:
#         return []
# 
# 
# class FeatureCounts(object):
#     """
#     Run featureCounts for GTF + BAM(s)
#     mission()
#     count.txt
#     determin: library type: +-/-+; -+/+-;
#     map reads, scale
#     FPKM/RPKM;
#     
#     
#     example:
#     FeatureCounts(gtf=a, bam_list=b, outdir=c).run()
#     gene file format: bed, gtf, saf
#     
#     bed_to_saf
#     """
#     def __init__(self, **kwargs):
#         """
#         arguments:
#         gtf:
#         bam_list: (index)
#         outdir:
#         strandness: 0=no, 1=sens, 2=anti, 3=both
#         smpname: count.txt default
#         threads: 4
#         overwrite: False
#         """
#         self = update_obj(self, kwargs, force=True)
#         self.init_args()
#         self.init_files()
#         Toml(self.__dict__).to_toml(self.config_toml)
# 
# 
#     def init_args(self):
#         args_init = {
#             'gtf': None,
#             'bam_list': None,
#             'outdir': None,
#             'prefix': None,
#             'strandness': 0,
#             'threads': 4,
#             'overwrite': False
#         }
#         self = update_obj(self, args_init, force=False)
# 
#         # gtf/bed/saf
#         if not isinstance(self.gtf, str):
#             raise ValueError('--gtf, expect str, got {}'.format(
#                 type(self.gtf).__name__))
# 
#         if not file_exists(self.gtf):
#             raise ValueError('--gtf, file not exists {}'.format(self.gtf))
# 
#         # bam files
#         if isinstance(self.bam_list, str):
#             self.bam_list = [self.bam_list]
#         elif isinstance(self.bam_list, list):
#             pass
#         else:
#             raise ValueError('bam_list, str or list, got {}'.format(
#                 type(self.bam_list).__name__))
# 
#         if not all(file_exists(self.bam_list)):
#             raise ValueError('--bam-list, file not exists:')
#         # index
#         [Bam(i).index() for i in self.bam_list]
# 
#         # output
#         if not isinstance(self.outdir, str):
#             # self.outdir = str(pathlib.Path.cwd())
#             self.outdir = self._tmp(dir=True)
#         check_path(self.outdir)
# 
#         # prefix
#         if not isinstance(self.prefix, str):
#             self.prefix = 'count.txt'
# 
#         # abs path
#         self.gtf = file_abspath(self.gtf)
#         self.bam_list = file_abspath(self.bam_list)
#         self.outdir = file_abspath(self.outdir)
# 
# 
#     def init_files(self):
#         self.config_dir = os.path.join(self.outdir, 'config')
# 
#         default_files = {
#             # 'config_txt': self.config_dir + '/config.txt',
#             # 'config_pickle': self.config_dir + '/config.pickle',
#             # 'config_json': self.config_dir + '/config.json',
#             'config_toml': self.config_dir + '/config.toml',
#             'count_txt': self.outdir + '/' + self.prefix,
#             'summary': self.outdir + '/' + self.prefix + '.summary',
#             'log': self.outdir + '/' + self.prefix + '.featureCounts.log',
#             'stat': self.outdir + '/' + self.prefix + '.featureCounts.stat'
#         }
#         self = update_obj(self, default_files, force=True) # key
#         check_path(self.config_dir)
# 
#         # check GTF/SAF file
#         # Convert BED to SAF
#         self.gtf_ext = os.path.splitext(self.gtf)[1].lower()
#         if self.gtf_ext in ['.bed', '.narrowpeak', '.broadpeak']:
#             self.saf = os.path.join(self.outdir, file_prefix(self.gtf)[0] + '.saf')
#             self.bed_to_saf(self.gtf, self.saf)
# 
# 
#     def _tmp(self, dir=False):
#         """
#         Create a tmp file to save json object
#         """
#         if dir:
#             tmp = tempfile.TemporaryDirectory()
#         else:
#             tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=False)
#         return tmp.name
# 
# 
#     def is_pe(self):
#         """
#         Check whether the input bam files are Paired or Single file
#         Bam().isPaired()
#         """
#         return all([Bam(i).isPaired() for i in self.bam_list])
# 
# 
#     def bed_to_saf(self, bed_in, saf_out):
#         """
#         SAF format: 
#         GeneID Chr Start End Strand
# 
#         see: https://www.biostars.org/p/228636/#319624
#         """
#         if file_exists(saf_out) and not self.overwrite:
#             log.info('bed_to_saf() skipped, file exists: {}'.format(saf_out))
#         else:
#             try:
#                 with open(bed_in, 'rt') as r, open(saf_out, 'wt') as w:
#                     for line in r:
#                         tabs = line.strip().split('\t')
#                         chr, start, end = tabs[:3]
#                         # strand
#                         strand = '.'
#                         if len(tabs) > 5:
#                             s = tabs[5]
#                             if s in ['+', '-', '.']:
#                                 strand = s
#                         # id
#                         if len(tabs) > 3:
#                             name = os.path.basename(tabs[3])
#                         else:
#                             name = '_'.join([chr, start, end, strand])
#                         # output
#                         w.write('\t'.join([name, chr, start, end, strand])+'\n')
#             except:
#                 log.error('bed_to_saf() failed, see: {}'.format(saf_out))
# 
# 
#     def get_arg_gtf(self):
#         """
#         gtf
#         saf
#         bed -> saf
#         """
#         # check gtf or saf or bed
#         if self.gtf_ext in ['.gtf']:
#             arg = '-a {} -F GTF -t exon -g gene_id '.format(self.gtf)
#         elif self.gtf_ext in ['.bed', '.narrowpeak', '.broadpeak']:
#             arg = '-a {} -F SAF'.format(self.saf) # SAF
#         elif self.gtf_ext in ['.saf']:
#             arg = '-a {} -F SAF'.format(self.gtf)
#         else:
#             arg = ''
# 
#         return arg
# 
# 
#     def get_arg_pe(self):
#         """
#         arguments for Paired-End fq
#         """
#         return '-p -C -B' if self.is_pe() else ''
# 
# 
#     def get_cmd(self):
#         """
#         prepare args for featureCounts
#         """
#         return ' '.join([
#             '{}'.format(shutil.which('featureCounts')),
#             '-s {}'.format(self.strandness),
#             self.get_arg_gtf(),
#             self.get_arg_pe(),
#             '-o {}'.format(self.count_txt),
#             '-T {}'.format(self.threads),
#             '-M -O --fraction',
#             ' '.join(self.bam_list),
#             '2> {}'.format(self.log)
#             ])
# 
# 
#     def wrap_log(self):
#         """
#         save output file to log,
#         """
#         try:
#             df = pd.read_csv(self.summary, '\t', index_col=0)
#             df.columns = list(map(os.path.basename, df.columns.to_list()))
#             total = df.sum(axis=0, skipna=True)
#             assign = df.loc['Assigned', ]
#             assign_pct = assign / total
#             assign_pct = assign_pct.round(decimals=4)
#             assign_df = assign_pct.to_frame('assigned')
#             assign_df['strandness'] = self.strandness
#             # save to json
#             assign_df.to_csv(self.stat, sep='\t', index=True, header=False)
#             assign_df.to_json(self.stat)
# 
#             # mimimal value
#             assign_min = assign_pct.min()
#             if assign_min < 0.50:
#                 log.warning('Caution: -s {}, {:.2f}% assigned, see {}'.format(
#                     self.strandness, assign_min, self.summary))
#             print(assign_df)
#         except:
#             log.warning('reading file failed: {}'.format(self.summary))
#             total, assign, assign_pct = [1, 0, 0]
# 
#         df = pd.DataFrame([total, assign, assign_pct]).T
#         df.columns = ['total', 'map', 'pct']
#         return df
# 
# 
#     def run(self):
#         """
#         run featureCounts
#         """
#         cmd = self.get_cmd()
# 
#         # save cmd
#         cmd_shell = os.path.join(self.outdir, self.prefix +'.cmd.sh')
#         with open(cmd_shell, 'wt') as w:
#             w.write(cmd + '\n')
# 
#         if file_exists(self.count_txt) and not self.overwrite:
#             log.info("FeatureCounts() skippped, file exists: {}".format(
#                 self.count_txt))
#         else:
#             try:
#                 run_shell_cmd(cmd)
#             except:
#                 log.error('FeatureCounts() failed, see: {}'.format(self.log))
# 
#         return self.wrap_log()
# 
# 
# 
# def symlink(src, dest, absolute_path=True):
#     """
#     Create symlinks within output dir
#     ../src
#     """
#     if src is None or dest is None:
#         log.warning('symlink skipped: {}, to: {}'.format(src, dest))
#     elif file_exists(dest):
#         log.warning('symlink skipped, target exists...'.format(dest))
#     else:
#         if absolute_path:
#             # support: ~, $HOME,
#             srcname = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
#         else:
#             # only for directories within the same folder
#             srcname = os.path.join('..', os.path.basename(src))
# 
#         if not os.path.exists(dest):
#             os.symlink(srcname, dest)
# 
