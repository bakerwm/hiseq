
# -*- coding: utf-8 -*-

"""
Common functions for pipeline construction
file modification,
...
"""

import os
import sys
import re
import gzip
import shutil
import json
import pickle
import fnmatch
import tempfile
import logging
import functools
import subprocess
import pysam
import pybedtools
import pathlib
import binascii
import datetime
import pandas as pd
from itertools import combinations
# from .args import args_init, ArgumentsInit
### local test ###
# from args import args_init, ArgumentsInit # for local test


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)


# Decorator:
class Logger(object):
    def __init__(self, level='INFO'):
        logging.basicConfig(
            format='[%(asctime)s %(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            stream=sys.stdout)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(level)

    def __call__(self, fn):
        @functools.wraps(fn)
        def decorated(*args, **kwargs):
            try:
                self.logger.info('{0} - {1} - {2}'.format(
                    fn.__name__,
                    args,
                    kwargs))
                result = fn(*args, **kwargs)
                self.logger.info(result)
                # return result
            except Exception as ex:
                self.logger.info('Exception {0}'.format(ex))
                raise ex
            return result
        return decorated


## 1. files and path ##
def listdir(path, full_name=True, recursive=False, include_dir=False):
    """
    List all the files within the path
    """
    out = []
    for root, dirs, files in os.walk(path):
        if full_name:
            dirs = [os.path.join(root, d) for d in dirs]
            files = [os.path.join(root, f) for f in files]
        out += files

        if include_dir:
            out += dirs

        if recursive is False:
            break

    return sorted(out)


def listfile(path='.', pattern='*', full_name=True, recursive=False):
    """
    Search files by the pattern, within directory
    fnmatch.fnmatch()

    pattern:

    *       matches everything
    ?       matches any single character
    [seq]   matches any character in seq
    [!seq]  matches any char not in seq

    An initial period in FILENAME is not special.
    Both FILENAME and PATTERN are first case-normalized
    if the operating system requires it.
    If you don't want this, use fnmatchcase(FILENAME, PATTERN).

    example:
    listfile('./', '*.fq')
    """
    fn_list = listdir(path, full_name, recursive, include_dir=False)
    fn_list = [f for f in fn_list if fnmatch.fnmatch(f, pattern)]
    return sorted(fn_list)


def list_fq_files(path, pattern='*'):
    """
    Parse fastq files within path, using the prefix
    PE reads
    SE reads

    *.fastq
    *.fq
    *.fastq.gz
    *.fq.gz

    _1.
    _2.
    """
    # all fastq files: *f[astq]+(.gz)?
    fq_list = listfile('*q.gz', path) # *fastq.gz, *fq.gz
    fq_list.extend(listfile('*q', path)) # *fastq, *fq

    # filter
    if pattern == '*':
        hit_list = fq_list
    else:
        p = re.compile(r'(_[12])?.f(ast)?q(.gz)?$')
        hit_list = [f for f in all_files if p.search(f) and x in f]

    # chk1
    if len(hit_list) == 0:
        log.error('no fastq files found: {}'.format(path))

    # determine SE or PE
    r0 = r1 = r2 = []
    for i in hit_list:
        p1 = re.compile(r'_[rR]?1.f(ast)?q(.gz)?$') # read1
        p2 = re.compile(r'_[rR]?2.f(ast)?q(.gz)?$') # read2
        if p1.search(i):
            r1.append(i)
        elif p2.search(i):
            r2.append(i)
        else:
            r0.append(i)

    # chk2
    if len(r2) > 0 and not len(r1) == len(r2):
        log.error('read1 and read2 not equal: \nread1: {}\nread2: {}'.format(r1, r2))

    # organize [[r1, r2], [r1, r2], ...]
    out_list = []
    for i, j in zip(r1, r2):
        out_list.append([i, j])

    for i in r0:
        out_list.append([i, None])

    # if len(r1) > 6:
    #     log.warning('too many records matched : {} \n{}'.format(pattern, r1))

    return out_list


def is_gz(filepath):
    if os.path.exists(filepath):
        with open(filepath, 'rb') as test_f:
            return binascii.hexlify(test_f.read(2)) == b'1f8b'
    else:
        if filepath.endswith('.gz'):
            return True
        else:
            return False


def file_prefix(fn, with_path=False):
    """
    extract the prefix of a file
    remove extensions
    .gz, .fq.gz
    """
    assert isinstance(fn, str)
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.endswith('gz') or px.endswith('.bz2'):
        px = os.path.splitext(p1)[1] + px
        p1 = os.path.splitext(p1)[0]
    if not with_path:
        p1 = os.path.basename(p1)
    return [p1, px]


def symlink(src, dest, absolute_path=True):
    """
    Create symlinks within output dir
    ../src
    """
    if src is None or dest is None:
        log.warning('symlink skipped: {}, to: {}'.format(src, dest))
    else:
        if absolute_path:
            # support: ~, $HOME,
            srcname = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        else:
            # only for directories within the same folder
            srcname = os.path.join('..', os.path.basename(src))

        if not os.path.exists(dest):
            os.symlink(srcname, dest)


def check_file(x, show_log=False):
    """
    if x (file, list) exists or not
    """
    if isinstance(x, str):
        flag = 'ok' if os.path.exists(x) else 'failed'
        if show_log is True:
            log.info('{:<6s} : {}'.format(flag, x))
        return os.path.exists(x)
    elif isinstance(x, list):
        return all([check_file(i, show_log=show_log) for i in x])
    else:
        log.warning('expect str and list, not {}'.format(type(x)))
        return None


def check_path(x, show_log=False, create_dirs=True):
    """
    Check if x is path, Create path
    """
    if isinstance(x, str):
        if os.path.isdir(x):
            tag = True
        else:
            if create_dirs is True:
                try:
                    os.makedirs(x)
                    tag = True
                except:
                    tag = False
            else:
                tag = False
        # show log
        flag = 'ok' if tag is True else 'failed'
        if show_log is True:
            log.info('{:<6s} : {}'.format(flag, x))
        return tag
    elif isinstance(x, list):
        return all([check_path(i, show_log, create_dirs) for i in x])
    else:
        log.warning('expect str and list, not {}'.format(type(x)))
        return None


def fq_name(fq, include_path=False, pe_fix=False):
    """
    parse the name of fastq file:
    .fq.gz
    .fastq.gz
    .fq
    .fastq
    (also for fasta, fa)
    """
    # if isinstance(x, str):
    #     fname = file_prefix(x)[0]
    #     fname = re.sub('[._][rR]?1$', '', fname)
    #     fname = re.sub('_\d$', '', fname)
    # elif isinstance(x, list):
    #     fname = [fq_name(f) for f in x]
    # elif x is None:
    #     log.warning('None type detected')
    #     fname = None
    # else:
    #     log.warning('unknown input detected')
    #     fname = None
    # return fname
    ###############
    # fastq.gz, fastq, fasta.gz, fa.gz
    # p1 = re.compile('(_[12])?[.](fast|f)[aq](.gz)?$', re.IGNORECASE)
    p1 = re.compile('[.](fast|f)[aq](.gz)?$', re.IGNORECASE)
    p2 = re.compile('[._][12]$', re.IGNORECASE)
    if isinstance(fq, str):
        fq = fq if include_path is True else os.path.basename(fq)
        fq_tmp = re.sub(p1, '', fq) # r1_1.fq.gz : r1_1
        if pe_fix is True:
            fq_tmp = re.sub(p2, '', fq_tmp) # r1_1.fq.gz: r1
        return fq_tmp
    elif isinstance(fq, list):
        return [fq_name(x, include_path=include_path, pe_fix=pe_fix) for x in fq]
    else:
        log.warning('unknown type found: {}'.format(type(fq)))
        return fq



def fq_name_rmrep(fq, include_path=False, pe_fix=True):
    """
    x, filename, or list

    Remove the *.rep[123], *.REP[123] from tail
    """
    # if isinstance(x, str):
    #     return fq_name(x).rstrip('rep|REP|r|R||_|.|1|2')
    # elif isinstance(x, list):
    #     return [fq_name_rmrep(i) for i in x]
    # else:
    #     pass
    ##################
    p1 = re.compile('[._](rep|r)[0-9]+$', re.IGNORECASE) # _rep1, _r1
    if isinstance(fq, str):
        return re.sub(p1, '', fq_name(fq, include_path, pe_fix))
    elif isinstance(fq, list):
        return [fq_name_rmrep(x, include_path, pe_fix) for x in fq]
    else:
        log.warning('unknown type found: {}'.format(type(fq)))
        return fq


def file_abspath(file):
    """
    Create os.path.abspath() for files, directories
    """
    if file is None:
        return None
    elif isinstance(file, str):
        return os.path.abspath(file)
    elif isinstance(file, list):
        return [file_abspath(x) for x in file]
    else:
        log.warning('unknown type found: {}'.format(type(file)))
        return file


def file_exists(file, isfile=True):
    """
    os.path.exists() for files, os.path.isfile()
    """
    if file is None:
        return None
    elif isinstance(file, str):
        chk1 = os.path.exists(file) 
        chk2 = os.path.isfile(file)
        chk3 = os.path.isdir(file) # candidate
        return chk1 and chk2 if isfile else chk1
    elif isinstance(file, list):
        return [file_exists(x, isfile) for x in file]
    else:
        log.warning('unknown type found: {}'.format(type(file)))
        return file


def merge_names(x):
    """
    Get the name of replictes
    common in left-most
    """
    assert isinstance(x, list)
    name_list = [os.path.basename(i) for i in x]
    name_list = [re.sub('.rep[0-9].*$', '', i) for i in name_list]
    return list(set(name_list))[0]


## 2. commandline ##
def run_shell_cmd(cmd):
    """This command is from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py
    """
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid,
        pgid,
        rc,
        stderr.strip(),
        stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')


def gzip_cmd(src, dest, decompress=True, rm=True):
    """
    Gzip Compress or Decompress files using gzip module in python
    rm, True/False, whether remove old file

    # check the src file by extension: .gz
    """
    if os.path.exists(dest):
        log.warning('file exists, skipped - {}'.format(dest))
    else:
        if decompress:
            if is_gz(src):
                with gzip.open(src, 'rb') as r, open(dest, 'wb') as w:
                    shutil.copyfileobj(r, w)
            else:
                log.warning('not a gzipped file: {}'.format(src))
                shutil.copy(src, dest)
        else:
            if is_gz(src):
                log.warning('input is gzipped file, no need gzip')
                shutil.copy(src, dest)
            else:
                with open(src, 'rb') as r, gzip.open(dest, 'wb') as w:
                    shutil.copyfileobj(r, w)

    # output
    if rm is True:
        os.remove(src)

    return dest


def list_uniquer(seq, sorted=True, idfun=None):
    """
    seq: a list with items
    sorted: whether sort the output(unique)

    get the unique of inlist

    see1: Markus
    remove duplicates from a list while perserving order
    https://stackoverflow.com/a/480227/2530783

    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    see2: ctcherry
    https://stackoverflow.com/a/89202/2530783

    def f5(seq, idfun=None):
        # order preserving
        if idfun is None:
            def idfun(x): return x
        seen = {}
        result = []
        for item in seq:
            marker = idfun(item)
            # in old Python versions:
            # if seen.has_key(marker)
            # but in new ones:
            if marker in seen: continue
            seen[marker] = 1
            result.append(item)
        return result
    """
    if idfun is None:
        def idfun(x): return x # for None

    if not isinstance(seq, list):
        log.error('list required, but get {}'.format(type(seq)))
        return [] # blank
    elif sorted is True:
        return list(set(seq))
    else:
        seen = set()
        return [x for x in seq if x not in seen and not seen.add(x)]


def args_checker(d, x, update=False):
    """Check if dict and x are consitent
    d is dict
    x is pickle file
    """
    assert isinstance(d, dict)
    flag = None
    if os.path.exists(x):
        # read file to dict
        with open(x, 'rb') as fh:
            d_checker = pickle.load(fh)
        if d == d_checker:
            flag = True
        else:
            if update:
                with open(x, 'wb') as fo:
                    pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    elif isinstance(x, str):
        # save dict to new file
        with open(x, 'wb') as fo:
            pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        log.error('illegal x= argument: %s' % x)

    return flag


def args_logger(d, x, overwrite=False):
    """Format dict, save to file
        key: value
    """
    assert isinstance(d, dict)
    n = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
    if os.path.exists(x) and overwrite is False:
        return True
    else:
        with open(x, 'wt') as fo:
            fo.write('\n'.join(n) + '\n')
        return '\n'.join(n)


## 3. utils ##
def sam_flag_check(query, subject):
    """
    Check two numbers, query (for filtering) in subject or not
    convert to binary mode
    q: 0000010  (2)
    s: 1011011  (91)
    q in s
    range: 0 - 2048 (SAM flag)
    """
    def to_bin(n):
        return '{0:012b}'.format(n)

    # convert to binary mode
    q = to_bin(eval(query))
    s = to_bin(eval(subject))

    # check q, s
    flag = True
    for j, k in zip(q[::-1], s[::-1]):
        if not j == '1':
            continue
        if eval(j) - eval(k) > 0:
            flag = False
            break

    return flag



################################################################################
## functions for pipeline

def in_dict(d, k):
    """
    Check the keys in dict or not
    """
    assert isinstance(d, dict)
    if isinstance(k, str):
        k_list = [k]
    elif isinstance(k, list):
        k_list = list(map(str, k))
    else:
        log.warning('expect str and list, not {}'.format(type(k)))
        return False

    return all([i in d for i in k_list])


def in_attr(x, a, return_values=True):
    """
    Check a (attributes) in object a or not
    return the values or not
    """
    if isinstance(a, str):
        a_list = [a]
    elif isinstance(a, list):
        a_list = list(map(str, a))
    else:
        log.warning('expect str and list, not {}'.format(type(a)))
        return False

    # status
    status = all([hasattr(x, i) for i in a_list])

    if status and return_values:
        # values
        return [getattr(x, i) for i in a_list]
    else:
        return status


def dict_to_log(d, x, overwrite=False):
    """
    Convert dict to log style
        key | value
    """
    assert isinstance(d, dict)
    logout = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
    if overwrite is True or not os.path.exists(x):
        with open(x, 'wt') as w:
            w.write('\n'.join(logout) + '\n')

    return '\n'.join(logout)


def dict_to_pickle(d, x):
    """
    Convert dict to pickle
    """
    assert isinstance(d, dict)
    with open(x, 'wb') as w:
        pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_to_dict(x):
    """
    Convert pickle file to dict
    """
    with open(x, 'rb') as r:
        return pickle.load(r)


class Dict2Obj(object):
    """
    >>> d = {'a': 1, 'b': 2}
    >>> b = Dict2obj(**d)
    >>> print(b.a)
    1

    class AttributeDict(object):
    A class to convert a nested Dictionary into an object with key-values
    accessibly using attribute notation (AttributeDict.attribute) instead of
    key notation (Dict["key"]). This class recursively sets Dicts to objects,
    allowing you to recurse down nested dicts (like: AttributeDict.attr.attr)
    from: http://databio.org/posts/python_AttributeDict.html
    by Nathan Sheffield
    """
    def __init__(self, **entries):
        self.add_entries(**entries)

    def add_entries(self, **entries):
        for key, value in entries.items():
            if type(value) is dict:
                self.__dict__[key] = AttributeDict(**value)
            else:
                self.__dict__[key] = value

    def __getitem__(self, key):
        """
        Provides dict-style access to attributes
        """
        return getattr(self, key)


class Dict2Class(object):
    """
    Turns a dictionary into a class
    from: https://www.blog.pythonlibrary.org/2014/02/14/python-101-how-to-change-a-dict-into-a-class/
    by Mike
    """
    #----------------------------------------------------------------------
    def __init__(self, d):
        """Constructor"""
        for k, v in d.items():
            setattr(self, k, v)

        # for key in d:
        #     setattr(self, key, d[key])


def bed2gtf(infile, outfile):
    """Convert BED to GTF
    chrom chromStart chromEnd name score strand
    """
    with open(infile) as r, open(outfile, 'wt') as w:
        for line in r:
            fields = line.strip().split('\t')
            start = int(fields[1]) + 1
            w.write('\t'.join([
                fields[0],
                'BED_file',
                'gene',
                str(start),
                fields[2],
                '.',
                fields[5],
                '.',
                'gene_id "{}"; gene_name "{}"'.format(fields[3], fields[3])
                ]) + '\n')
    return outfile


### deprecated - BEGIN ###
def listfiles(path, full_name=True, recursive=False, include_dir=False):
    """
    List all the files within the path
    """
    out = []
    for root, dirs, files in os.walk(path):
        if full_name:
            dirs = [os.path.join(root, d) for d in dirs]
            files = [os.path.join(root, f) for f in files]
        out += files

        if include_dir:
            out += dirs

        if recursive is False:
            break
    return out


def listfiles2(pattern, path='.', full_name=True, recursive=False):
    """
    List all the files in specific directory
    fnmatch.fnmatch()

    pattern:

    *       matches everything
    ?       matches any single character
    [seq]   matches any character in seq
    [!seq]  matches any char not in seq

    An initial period in FILENAME is not special.
    Both FILENAME and PATTERN are first case-normalized
    if the operating system requires it.
    If you don't want this, use fnmatchcase(FILENAME, PATTERN).

    example:
    listfiles('*.fq', './')
    """
    fn_list = listfiles(path, full_name, recursive, include_dir=False)
    fn_list = [f for f in fn_list if fnmatch.fnmatch(f, pattern)]
    return fn_list


def is_path(path, create = True):
    """
    Check path, whether a directory or not
    if not, create it
    """
    assert isinstance(path, str)
    if os.path.exists(path):
        return True
    else:
        if create:
            try:
                os.makedirs(path)
                return True
            except IOError:
                log.error('failed to create directories: %s' % path)
        else:
            return False
### deprecated - END ###


class Json(object):
    """
    Wrapper for *.json file
    1.json to dict
    2.dict to json -> file
    """
    def __init__(self, input, **kwargs):
        self.input = input
        self.mission()


    def mission(self):
        """
        Check what to do, based on the input args
        """
        # json_file to dict
        if input is None:
            log.warning('require, dict or str (file), Nonetype detected')
        elif isinstance(self.input, str) and os.path.exists(self.input):
            self.json_file_in = self.input
            # self.reader() # updated
        elif isinstance(self.input, dict):
            self.d = self.input
            # self.writer()
        else:
            log.warning('expect dict or str, get {}'.format(type(self.input).__name__))


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json',
            delete=False)
        return tmp.name


    def to_dict(self):
        # self.json_file_in = self.input
        return self.reader()


    def to_json(self, file=None):
        return self.writer(file)


    def reader(self):
        """
        Read json file as dict
        """
        if isinstance(self.input, dict):
            self.d = self.input
        else:
            try:
                with open(self.json_file_in) as r:
                    if os.path.getsize(self.json_file_in) > 0:
                        self.d = json.load(r)
                    else:
                        self.d = {}
            except:
                log.warning('failed reading json file: {}'.format(self.json_file_in))
                self.d = {}

        return self.d


    def writer(self, json_file_out=None):
        """
        Write d (dict) to file x, in json format
        """
        # save to file
        if json_file_out is None:
            json_file_out = self._tmp()

        try:
            with open(json_file_out, 'wt') as w:
                json.dump(self.d, w, indent=4, sort_keys=True)
        except:
            log.warning('failed saving file: {}'.format(json_file_out))

        return json_file_out



# class Json(object):

#     def __init__(self, x):
#         """
#         x
#           - dict, save to file
#           - json, save to file
#           - file, read as dict
#         Save dict to json file
#         Read from json file as dict
#         ...
#         """
#         self.x = x # input

#         if isinstance(x, Json):
#             self.dict = x.dict
#         elif isinstance(x, dict):
#             # input a dict,
#             # save to file
#             self.dict = x
#         elif os.path.exists(x):
#             # a file saving json content
#             self.dict = self.reader()
#         else:
#             raise Exception('unknown objec: {}'.format(x))

#     def _tmp(self):
#         """
#         Create a tmp file to save json object
#         """
#         tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json',
#             delete=False)
#         return tmp.name


#     def reader(self):
#         """
#         Read json file as dict
#         """
#         if os.path.getsize(self.x) > 0:
#             with open(self.x) as r:
#                 d = json.load(r)
#         else:
#             d = {}

#         return d


#     def writer(self, f=None):
#         """
#         Write d (dict) to file x, in json format
#         """
#         # save to file
#         if f is None:
#             f = self._tmp()

#         assert isinstance(f, str)
#         # assert os.path.isfile(f)

#         if isinstance(self.x, dict):
#             with open(f, 'wt') as w:
#                 json.dump(self.x, w, indent=4, sort_keys=True)

#         return f


class Genome(object):
    """
    List related information of specific genome
    1. get_fa(), genome fasta
    2. get_fasize(), genome fasta size
    3. bowtie_index(), bowtie index, optional, rRNA=True
    4. bowtie2_index(), bowtie2 index, optional, rRNA=True
    5. star_index(), STAR index, optional, rRNA=True
    6. gene_bed(),
    7. gene_rmsk(),
    8. gene_gtf(), optional, version='ucsc|ensembl|ncbi'
    9. te_gtf(), optional, version='ucsc'
    10. te_consensus(), optional, fruitfly()
    ...

    directory structure of genome should be like this:
    /path-to-data/{genome}/
        |- bigZips  # genome fasta, fasize, chromosome
        |- annotation_and_repeats  # gtf, bed, rRNA, tRNA, annotation
        |- bowtie_index
        |- bowtie2_index
        |- STAR_index
        |- hisat2_index
        |- phylop100
        |- ...

    default: $HOME/data/genome/{genome}

    """
    def __init__(self, genome, genome_path=None,
        repeat_masked_genome=False, **kwargs):
        assert isinstance(genome, str)
        self.genome = genome
        self.repeat_masked_genome = repeat_masked_genome
        self.kwargs = kwargs

        # path
        if genome_path is None:
            genome_path = os.path.join(str(pathlib.Path.home()),
                'data', 'genome')
        self.genome_path = genome_path

        # # deprecated
        # if not supportedGenome(genome):
        #     log.error('genome not supported: {}'.foramt(genome))


    def get_fa(self):
        """
        Get the fasta file of specific genome
        {genome}/bigZips/{genome}.fa
        also check ".gz" file
        """
        fa = os.path.join(self.genome_path, self.genome, 'bigZips',
            self.genome + '.fa')
        if not os.path.exists(fa):
            # gencode version
            fa = os.path.join(self.genome_path, self.genome, 'fasta',
                self.genome + '.fa')

        fa_gz = fa + '.gz'
        if not os.path.exists(fa):
            if os.path.exists(fa_gz):
                log.error('require to unzip the fasta file: %s' % fa_gz)
            else:
                log.error('fasta file not detected: %s' % fa)
            return None
        else:
            return fa


    def get_fasize(self):
        """Get the fasta size file, chromosome size
        optional, fetch chrom size from ucsc
        http://hgdownload.cse.ucsc.edu/goldenPath/<db>/bigZips/<db>.chrom.sizes

        or using UCSC tool: fetchChromSizes
        fetchChromSizes hg39 > hg38.chrom.sizes
        """
        fa = self.get_fa()
        fa_size = fa + '.chrom.sizes'

        if not os.path.exists(fa_size):
            # log.info('Downloading chrom.sizes from UCSC')
            # url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes' % self.genome
            # download(url, fa_size)
            log.warning('file not exists, run samtools faidx to generate it')
            pysam.faidx(fa) # create *.fa.fai
            os.rename(fa + '.fai', fa_size)

        return fa_size


    def phylop100(self):
        """Return the phylop100 bigWig file of hg19, only
        for conservation analysis
        """
        p = os.path.join(self.genome_path, self.genome, 'phyloP100way',
            self.genome + '.100way.phyloP100way.bw')
        if not os.path.exists(p):
            p = None
        return p


    def gene_bed(self, version='refseq', rmsk=False):
        """Return the gene annotation in BED format
        support UCSC, ensembl, gencode
        """
        if rmsk:
            suffix = '.rmsk.bed'
        else:
            suffix = '.refseq.bed'
        g = os.path.join(self.genome_path, self.genome,
            'annotation_and_repeats', self.genome + suffix)
        if not os.path.exists(g):
            g = None
        return g


    def gene_gtf(self, version='refseq'):
        """Return the gene annotation in GTF format
        support refseq, ensembl, gencode
        """
        version = version.lower() #

        gtf = os.path.join(
            self.genome_path,
            self.genome,
            'annotation_and_repeats',
            self.genome + '.' + version + '.gtf')

        if not os.path.exists(gtf):
            gtf = os.path.join(
            self.genome_path,
            self.genome,
            'gtf',
            self.genome + '.' + version + '.gtf')

        if not os.path.exists(gtf):
            gtf = None

        return gtf


    def te_gtf(self, format='gtf'):
        """Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        # only dm3 supported
        te_gtf = os.path.join(self.genome_path, self.genome,
            self.genome + '_transposon',
            self.genome + '_transposon.gtf')
        if not os.path.exists(te_gtf):
            te_gtf = None

        return te_gtf


class Bam(object):
    """
    Manipulate BAM files
    - sort
    - index
    - merge
    - count
    - to_bed
    - rmdup
    - ...

    Using Pysam, Pybedtools, ...

    code from cgat:
    """
    def __init__(self, infile, threads=4):
        self.bam = infile
        self.threads = threads
        # self.bed = self.to_bed()


    def index(self):
        """
        Create index for bam
        """
        bai = self.bam + '.bai'
        if not os.path.exists(bai):
            pysam.index(self.bam)

        return os.path.exists(bai)


    def sort(self, outfile=None, by_name=False, overwrite=False):
        """
        Sort bam file by position (default)
        save to *.sorted.bam (or specify the name)
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.sorted.bam'

        if os.path.exists(outfile) and overwrite is False:
            logging.info('file exists: {}'.format(outfile))
        else:
            tmp = pysam.sort('-@', str(self.threads), '-o', outfile, self.bam)

        return outfile


    def merge(self):
        """
        Merge multiple BAM files using samtools
        """

        # pysam.merge('')
        pass


    def count(self):
        """Using samtools view -c"""
        return pysam.view('-c', self.bam)


    def to_bed(self, outfile=None):
        """Convert BAM to BED
        pybetools
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.bed'

        if not os.path.exists(outfile):
            pybedtools.BedTool(self.bam).bam_to_bed().saveas(outfile)

        return outfile


    def rmdup(self, outfile=None, overwrite=False):
        """
        Remove duplicates using picard/sambamba
        sambamba markdup -r --overflow-list-size 800000 raw.bam rmdup.bam
        picard MarkDuplicates  REMOVE_SEQUENCING_DUPLICATES=True I=in.bam O=outfile.bam
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.rmdup.bam'

        sambamba = shutil.which('sambamba')
        log = outfile + '.sambamba.log'
        cmd = '{} markdup -r -t {} --overflow-list-size 800000 \
            --tmpdir="./" {} {} 2> {}'.format(
            sambamba,
            str(self.threads),
            self.bam,
            outfile,
            log)

        if os.path.exists(outfile) and overwrite is False:
            logging.info('file exists: {}'.format(outfile))
        else:
            run_shell_cmd(cmd)

        return outfile


    def proper_pair(self, outfile=None, overwrite=False):
        """
        Extract proper pair
        samtools view -f 2
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.proper_pair.bam'

        if os.path.exists(outfile) and overwrite is False:
            logging.info('file exists: {}'.format(outfile))
        else:
            pysam.view('-f', '2', '-h', '-b', '-@', str(self.threads),
                '-o', outfile, self.bam, catch_stdout=False)

        return outfile


    ##########################################
    ## code from cgat: BEGIN
    ##########################################
    def isPaired(self, topn=1000):
        """
        Check if infile contains paired end reads

        go through the topn alignments in file,

        return: True, any of the alignments are paired
        """
        samfile = pysam.AlignmentFile(self.bam)
        n = 0
        for read in samfile:
            if read.is_paired:
                break
            n += 1
            if n == topn:
                break

        samfile.close()

        return n != topn


    def getNumReads(self):
        """
        Count number of reads in bam file.

        This methods works through pysam.idxstats.

        Arguments
        ---------
        bamfile : string
            Filename of :term:`bam` formatted file. The file needs
            to be indexed.
        Returns
        -------
        nreads : int
            Number of reads
        """

        lines = pysam.idxstats(self.bam).splitlines()

        try:
            nreads = sum(
                map(int, [x.split("\t")[2]
                          for x in lines if not x.startswith("#")]))

        except IndexError as msg:
            raise IndexError(
                "can't get number of reads from bamfile, msg=%s, data=%s" %
                (msg, lines))
        return nreads


    def estimateInsertSizeDistribution(self, topn=10000, n=10,
        method="picard", similarity_threshold=1.0, max_chunks=1000):
        """
        Estimate insert size from a subset of alignments in a bam file.

        Several methods are implemented.

        picard
            The method works analogous to picard by restricting the estimates
            to a core distribution. The core distribution is defined as all
            values that lie within n-times the median absolute deviation of
            the full data set.
        convergence
            The method works similar to ``picard``, but continues reading
            `alignments` until the mean and standard deviation stabilize.
            The values returned are the median mean and median standard
            deviation encountered.

        The method `convergence` is suited to RNA-seq data, as insert sizes
        fluctuate siginificantly depending on the current region
        being looked at.

        Only mapped and proper pairs are considered in the computation.

        Returns
        -------
        mean : float
           Mean of insert sizes.
        stddev : float
           Standard deviation of insert sizes.
        npairs : int
           Number of read pairs used for the estimation
        method : string
           Estimation method
        similarity_threshold : float
           Similarity threshold to apply.
        max_chunks : int
           Maximum number of chunks of size `alignments` to be used
           in the convergence method.

        """

        assert self.isPaired(self.bam), \
            'can only estimate insert size from' \
            'paired bam files'

        samfile = pysam.AlignmentFile(self.bam)

        def get_core_distribution(inserts, n):
            # compute median absolute deviation
            raw_median = numpy.median(inserts)
            raw_median_dev = numpy.median(numpy.absolute(inserts - raw_median))

            # set thresholds
            threshold_min = max(0, raw_median - n * raw_median_dev)
            threshold_max = raw_median + n * raw_median_dev

            # define core distribution
            return inserts[numpy.logical_and(inserts >= threshold_min,
                                             inserts <= threshold_max)]

        if method == "picard":

            # only get first read in pair to avoid double counting
            inserts = numpy.array(
                [read.template_length for read in samfile.head(n=topn)
                 if read.is_proper_pair
                 and not read.is_unmapped
                 and not read.mate_is_unmapped
                 and not read.is_read1
                 and not read.is_duplicate
                 and read.template_length > 0])
            core = get_core_distribution(inserts, n)

            return numpy.mean(core), numpy.std(core), len(inserts)

        elif method == "convergence":

            means, stds, counts = [], [], []
            last_mean = 0
            iteration = 0
            while iteration < max_chunks:

                inserts = numpy.array(
                    [read.template_length for read in samfile.head(
                        n=topn,
                        multiple_iterators=False)
                     if read.is_proper_pair
                     and not read.is_unmapped
                     and not read.mate_is_unmapped
                     and not read.is_read1
                     and not read.is_duplicate
                     and read.template_length > 0])
                core = get_core_distribution(inserts, n)
                means.append(numpy.mean(core))
                stds.append(numpy.std(core))
                counts.append(len(inserts))
                mean_core = get_core_distribution(numpy.array(means), 2)
                mm = numpy.mean(mean_core)
                if abs(mm - last_mean) < similarity_threshold:
                    break
                last_mean = mm

            return numpy.median(means), numpy.median(stds), sum(counts)
        else:
            raise ValueError("unknown method '%s'" % method)


    def estimateTagSize(self, topn=10, multiple="error"):
        """
        Estimate tag/read size from first alignments in file.

        Arguments
        ---------
        bamfile : string
           Filename of :term:`bam` formatted file
        alignments : int
           Number of alignments to inspect
        multiple : string
           How to deal if there are multiple tag sizes present.
           ``error`` will raise a warning, ``mean`` will return the
           mean of the read lengths found. ``uniq`` will return a
           unique list of read sizes found. ``all`` will return all
           read sizes encountered.

        Returns
        -------
        size : int
           The read size (actual, mean or list of read sizes)

        Raises
        ------
        ValueError
           If there are multiple tag sizes present and `multiple` is set to
           `error`.

        """
        samfile = pysam.AlignmentFile(self.bam)
        sizes = [read.rlen for read in samfile.head(topn)]
        mi, ma = min(sizes), max(sizes)

        if mi == 0 and ma == 0:
            sizes = [read.inferred_length for read in samfile.head(alignments)]
            # remove 0 sizes (unaligned reads?)
            sizes = [x for x in sizes if x > 0]
            mi, ma = min(sizes), max(sizes)

        if mi != ma:
            if multiple == "error":
                raise ValueError('multiple tag sizes in %s: %s' % (bamfile, sizes))
            elif multiple == "mean":
                mi = int(sum(sizes) / len(sizes))
            elif multiple == "uniq":
                mi = list(sorted(set(sizes)))
            elif multiple == "all":
                return sizes

        return mi


    def getNumberOfAlignments(self):
        """return number of alignments in bamfile.
        """
        samfile = pysam.AlignmentFile(self.bam)
        return samfile.mapped
    ##########################################
    ## code from cgat: END
    ##########################################


# ## for index
# ## to-do
# ##   - build index (not recommended)
# ##
# class AlignIndex(object):

#     def __init__(self, aligner='bowtie', index=None, **kwargs):
#         """
#         Required args:
#           - aligner
#           - index (optional)
#           - genome
#           - group : genome, rRNA, transposon, piRNA_cluster, ...
#           - genome_path
#         """
#         ## init
#         # args = args_init(kwargs, align=True) # init
#         # args = ArgumentsInit(kwargs, align=True)
#         self.args = ArgumentsInit(kwargs, trim=True).dict.__dict__
#         self.args.pop('args_input', None)
#         self.args.pop('cmd_input', None)
#         self.args.pop('dict', None)

#         self.aligner = aligner

#         if isinstance(index, str):
#             # index given
#             self.index = index.rstrip('/') #
#             self.name = self.get_name()
#             self.aligner_supported = self.get_aligner() # all
#             self.check = index if self.is_index() else None
#         else:
#             # index not defined, search required
#             # log.warning('index=, not defined; .search() required')
#             pass

#         # self.kwargs = args


#     def get_aligner(self, index=None):
#         """
#         Search the available index for aligner:
#         bowtie, [*.[1234].ebwt,  *.rev.[12].ebwt]
#         bowtie2, [*.[1234].bt2, *.rev.[12].bt2]
#         STAR,
#         bwa,
#         hisat2,
#         """
#         if index is None:
#             index = self.index

#         bowtie_files = [index + i for i in [
#             '.1.ebwt',
#             '.2.ebwt',
#             '.3.ebwt',
#             '.4.ebwt',
#             '.rev.1.ebwt',
#             '.rev.2.ebwt']]

#         bowtie2_files = [index + i for i in [
#             '.1.bt2',
#             '.2.bt2',
#             '.3.bt2',
#             '.4.bt2',
#             '.rev.1.bt2',
#             '.rev.2.bt2']]

#         hisat2_files = ['{}.{}.ht2'.format(index, i) for i in range(1, 9)]


#         bwa_files = [index + i for i in [
#             '.sa',
#             '.amb',
#             '.ann',
#             '.pac',
#             '.bwt']]

#         STAR_files = [os.path.join(index, i) for i in [
#             'SAindex',
#             'Genome',
#             'SA',
#             'chrLength.txt',
#             'chrNameLength.txt',
#             'chrName.txt',
#             'chrStart.txt',
#             'genomeParameters.txt']]

#         ## check exists
#         bowtie_chk = [os.path.exists(i) for i in bowtie_files]
#         bowtie2_chk = [os.path.exists(i) for i in bowtie2_files]
#         hisat2_chk = [os.path.exists(i) for i in hisat2_files]
#         bwa_chk = [os.path.exists(i) for i in bwa_files]
#         STAR_chk = [os.path.exists(i) for i in STAR_files]

#         ## check file exists
#         aligner = []

#         if all(bowtie_chk):
#             aligner.append('bowtie')
#         elif all(bowtie2_chk):
#             aligner.append('bowtie2')
#         elif all(hisat2_chk):
#             aligner.append('hisat2')
#         elif all(bwa_chk):
#             aligner.append('bwa')
#         elif all(STAR_chk):
#             aligner.append('STAR')
#         else:
#             pass

#         return aligner


#     def is_index(self, index=None):
#         """
#         Check if index support for aligner
#         """
#         if index is None:
#             index = self.index
#         return self.aligner in self.get_aligner(index)

#         # return self.aligner in self.aligner_supported if
#         #    index is None else self.aligner in self.get_aligner(index)


#     def get_name(self, index=None):
#         """
#         Get the name of index
#         basename: bowtie, bowtie2, hisqt2, bwa
#         folder: STAR
#         """
#         if index is None:
#             index = self.index
#         if os.path.isdir(index):
#             # STAR
#             iname = os.path.basename(index)
#             # iname = os.path.basename(os.path.dirname(index))
#         elif os.path.basename(index) == 'genome':
#             # ~/data/genome/dm3/bowtie2_index/genome
#             # bowtie, bowtie2, bwa, hisat2
#             # iname = os.path.basename(index)
#             iname = os.path.basename(os.path.dirname(os.path.dirname(index)))
#         else:
#             iname = os.path.basename(index)

#         return iname


#     def search(self, genome=None, group='genome'):
#         """
#         Search the index for aligner: STAR, bowtie, bowtie2, bwa, hisat2
#         para:

#         *genome*    The ucsc name of the genome, dm3, dm6, mm9, mm10, hg19, hg38, ...
#         *group*      Choose from: genome, rRNA, transposon, piRNA_cluster, ...

#         structure of genome_path:
#         default: {HOME}/data/genome/{genome_version}/{aligner}/

#         path-to-genome/
#             |- Bowtie_index /
#                 |- genome
#                 |- rRNA
#                 |- MT_trRNA
#             |- transposon
#             |- piRNA cluster

#         """
#         args = self.args.copy()

#         g_path = args.get('genome_path', './')
#         i_prefix = os.path.join(g_path, genome, self.aligner + '_index')

#         # all index
#         d = {
#             'genome': [
#                 os.path.join(i_prefix, 'genome'),
#                 os.path.join(i_prefix, genome)],
#             'genome_rm': [
#                 os.path.join(i_prefix, 'genome_rm'),
#                 os.path.join(i_prefix, genome + '_rm')],
#             'MT_trRNA': [
#                 os.path.join(i_prefix, 'MT_trRNA')],
#             'rRNA': [
#                 os.path.join(i_prefix, 'rRNA')],
#             'chrM': [
#                 os.path.join(i_prefix, 'chrM')],
#             'structural_RNA': [
#                 os.path.join(i_prefix, 'structural_RNA')],
#             'te': [
#                 os.path.join(i_prefix, 'transposon')],
#             'piRNA_cluster': [
#                 os.path.join(i_prefix, 'piRNA_cluster')],
#             'miRNA': [
#                 os.path.join(i_prefix, 'miRNA')],
#             'miRNA_hairpin': [
#                 os.path.join(i_prefix, 'miRNA_hairpin')]}

#         # hit
#         i_list = d.get(group, ['temp_temp'])
#         i_list = [i for i in i_list if self.is_index(i)]

#         return i_list[0] if len(i_list) > 0 else None


################################################################################
## functions


################################################################################
## TEMP

def featureCounts_reader(x, bam_names=False):
    """
    Read fc .txt as data.frame
    """
    df = pd.read_csv(x, '\t', comment='#')
    if bam_names:
        bam_list = df.columns.to_list()[6:] # from 7-column to end
        bam_list = [os.path.splitext(i)[0] for i in bam_list]
        return [os.path.basename(i) for i in bam_list]
    else:
        return df


## design, create combinations
def design_combinations(seq, n=2, return_index=True):
    """
    seq: the group name for each sample

    example:
    seq = ['ctl', 'ctl', 'exp', 'exp'] # keep order

    : return_index=True
    [[0, 1], [2, 3]]

    : return_index=False
    ['ctl', 'exp']
    """
    seq_unique = list_uniquer(seq, sorted=False) # keep order

    if len(seq_unique) >= n:
        item_pairs = list(combinations(seq_unique, n))
        # for index
        index_pairs = []
        for (a, b) in item_pairs:
            index_a = [i for i, x in enumerate(seq) if x == a]
            index_b = [i for i, x in enumerate(seq) if x == b]
            index_pairs.append([index_a, index_b])
        return index_pairs if return_index else item_pairs
    else:
        return []
