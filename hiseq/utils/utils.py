#!/usr/bin/env python3

"""Functions for string, list, dict, ...
transform to another version 
attribution of the vectors
"""



import os
import sys
import json
import yaml
import toml
import pickle
import glob
import logging
import collections
import subprocess
import random
import string
from PIL import Image # pillow
from datetime import datetime
from dateutil import tz
from itertools import combinations
from hiseq.utils.file import file_exists, list_dir
# from .helper import log
# from .file import file_exists, file_abspath
from difflib import SequenceMatcher



logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


## Temp
def find_longest_common_str(s1, s2):
    if isinstance(s1, str) and isinstance(s2, str):
        m = SequenceMatcher(None, s1, s2) # match
        l = m.find_longest_match(0, len(s1), 0, len(s2))
        out = s1[l.a:(l.a+l.size)]
    else:
        log.error('only support str, got s1={} s2={}'.type(
            type(s1).__name__, type(s2).__name__))
        out = None
    return out




################################################################################
## functions for hiseq-global

def list_hiseq_file(x, keys='bam', hiseq_type='r1'):
    """
    Return the specific file from project_dir of hiseq

    Parameters
    ----------
    x:  str
        Directory of the hiseq project

    keys:  str
        keys of the file in hiseq project, eg: 'bam', 'align_dir', ...
        default: ['bam']

    hiseq_type:  str
        Type of the hiseq, ['r1', 'rn', 'rx', 'rt', 'rp']
        see: HiseqReader() for details; default: ['r1']
    """
    out = []
    prj_dirs = list_hiseq_dir(x, hiseq_type) # all dirs
    for d in prj_dirs:
        a = read_hiseq(d)
        out.append(getattr(a, keys, None))
    return out



def list_hiseq_dir(x, hiseq_type='r1'):
    """
    Return the project_dir of hiseq

    Parameters
    ----------
    x:  str
        Directory of the hiseq project

    hiseq_type:  str
        Type of the hiseq, ['r1', 'rn', 'rx', 'rt', 'rp']
        see: HiseqReader() for details; default: ['r1']
    """
    a = read_hiseq(x)
    out = []
    if a.is_hiseq:
        if hiseq_type == 'r1':
            if a.is_hiseq_r1:
                out.append(x)
            elif a.is_hiseq_rn:
                # ATACseq; CnR, ChIPseq, ... !!!!to-do
                for i in a.rep_list:
                    out.extend(list_hiseq_dir(i, 'r1'))
            elif a.is_hiseq_rx:
                # list all directories in project
                dir_list = list_dir(a.project_dir, include_dir=True)
                for i in dir_list:
                    out.extend(list_hiseq_dir(i, 'r1'))
            else:
                pass
        elif hiseq_type == 'rn':
            if a.is_hiseq_rn:
                out.append(x)
            elif a.is_hiseq_rx:
                # list all directories in project
                dir_list = list_dir(a.project_dir, include_dir=True)
                for i in dir_list:
                    out.extend(list_hiseq_dir(i, 'rn'))
            else:
                pass
        elif hiseq_type == 'rx':
            if a.is_hiseq_rx:
                out.append(x)
        else:
            pass
    out = list(set(out)) # remove duplicates
    return sorted(out)


def read_hiseq(x, hiseq_type=True):
    """
    Parameters
    ---------
    x:  str
        Path to Hiseq directory

    hiseq_type: str
        Check the x is one-of-hiseq_types, could be head/tail;
        atacseq_r1, r1, rn, rnaseq_rx, ...
        default: [True], do

    Read config from hiseq directory
    """
    a = HiseqReader(x)
    if hasattr(a, 'hiseq_type'):
        a_hiseq_type = a.hiseq_type
        k = False
        if hiseq_type is True:
            k = True #pass
        elif isinstance(hiseq_type, str):
            k = any([
                a_hiseq_type == hiseq_type,
                a_hiseq_type.startswith(hiseq_type),
                a_hiseq_type.endswith(hiseq_type),
            ])
        else:
            pass
        # check
        if not k:
            log.warning('hiseq_dir not match, expect {}, got {}'.format(
                hiseq_type, a_hiseq_type))
    else:
        # raise ValueError('not a hiseq dir: {}'.format(x))
        # log.warning('not a hiseq dir: {}'.format(x))
        pass # reduce messagage #
    return a


class HiseqReader(object):
    """
    Parameters
    ---------
    x: str
        Path to Hiseq directory
    """
    def __init__(self, x):
        self.x = x
        self.init_args()


    def init_args(self):
        self.is_hiseq = False # init
        d = self.load()
        if isinstance(d, dict):
            hiseq_types = [
                'atacseq_type', 'rnaseq_type', 'hiseq_type',
                'align_type'
            ] # !!!! to-be-update
            # which hiseq
            hiseq_type = None
            for h in hiseq_types:
                if h in d:
                    hiseq_type = h
                    break
            # which type
            if isinstance(hiseq_type, str):
                a = d.get(hiseq_type, None)
            else:
                a = None
            # check output
            if isinstance(a, str):
                self.is_hiseq = True
                self.hiseq_type = a
                self.is_hiseq_r1 = a.endswith('_r1') # fq
                self.is_hiseq_rn = a.endswith('_rn') # group
                self.is_hiseq_rx = a.endswith('_rx') # group vs group
                self.is_hiseq_rt = a.endswith('_rt') # !!
                self.is_hiseq_rp = a.endswith('_rp') # report
            else:
                log.warning('unknown hiseq dir: {}'.format(self.x))
        else:
            # log.warning('not a hiseq dir: {}'.format(self.x))
            pass # reduce message log
        # update args
        self.args = update_obj(self, d, force=True)


    def list_config(self):
        """
        List the config files
        Support: toml, pickle, yaml, json, ... [priority]
        return file list
        # hiseq
        hiseq
          |-config
          |   |-config.toml

        # alignment
        align_dir
          |- smp_nmae
          |    |- index
          |    |    |- config.pickle
        """
        c_files = ['config.' + i for i in ['toml', 'pickle', 'yaml', 'json']]
        # search config files
        out = None
        for f in c_files:
            c1 = os.path.join(self.x, 'config', f)
            c2 = os.path.join(self.x, '*', '*', f)
            c1x = glob.glob(c1)
            c2x = glob.glob(c2)
            if len(c1x) > 0:
                out = c1x[0]
                break
            elif len(c2x) > 0:
                out = c2x[0]
                break
            else:
                continue
        return out


    def load(self):
        if isinstance(self.x, str) and file_exists(self.x):
            config_file = self.list_config()
            out = Config().load(config_file)
        else:
            out = None
        return out


################################################################################
## tmp functions
def gen_random_string(slen=10):
    return ''.join(random.sample(string.ascii_letters + string.digits, slen))


def print_dict(d):
    d = dict(sorted(d.items(), key=lambda x:x[0]))
    # d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


def init_cpu(threads=1, parallel_jobs=1):
    """
    The number of threads, parallel_jobs
    """
    n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()
    max_jobs = int(n_cpu / 4.0)
    ## check parallel_jobs (max: 1/4 of n_cpus)
    if parallel_jobs > max_jobs: 
        log.warning('Too large, change parallel_jobs from {} to {}'.format(
            parallel_jobs, max_jobs))
        parallel_jobs = max_jobs
    ## check threads
    max_threads = int(0.8 * n_cpu / parallel_jobs)
    if threads * parallel_jobs > 0.8 * n_cpu:
        log.warning('Too large, change threads from {} to {}'.format(
            threads, max_threads))
        threads = max_threads
    return (threads, parallel_jobs)


################################################################################












def get_date(timestamp=False):
    """
    Return the current date in UTC.timestamp or local-formated-string
    calculation in UTC
    return in local
    
    switch between timezone by datautil.tz
    
    Example:
    >>> from datetime import datetime
    >>> from dateutil import tz
    >>> get_date()
    '2021-05-18 17:08:53'
    
    >>> get_date(True)
    1621328957.280303
    
    # convert timestamp to local time
    >>> ts = get_date(True)
    >>> datetime.fromtimestamp(ts, tz.tzlocal())
    
    Arguments
    ---------
    ts:  float
        The timestamp (microseconds), if not, get the current time
    
    in_seconds:  bool
        Return the date in seconds (microseconds?!)
    """
    now = datetime.now(tz.tzlocal())
    if isinstance(timestamp, bool) and timestamp:
        out = now.timestamp()
    else:
        now = now.astimezone(tz.tzlocal()) # to local
        out = now.strftime('%Y-%m-%d %H:%M:%S') # YY-mm-dd H:M:S
    return out


def run_shell_cmd(cmd):
    """This command is from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py

    save log to file
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
            log.error(err_str)
    return (rc, stdout.strip('\n'), stderr.strip('\n'))


def is_cmd(x):
    """Check if the executable command"""
    return lambda i: shutil.which(i) is not None


##########################
## check files          ##
##########################
def update_obj(obj, d, force=True, remove=False):
    """Update the object, by dict
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)
    return obj



def unique_list(seq, sorted=True, idfun=None):
    """Get the unique items of list
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
        out = [] # blank
    elif sorted is True:
        out = list(set(seq))
    else:
        seen = set()
        out = [x for x in seq if x not in seen and not seen.add(x)]
    return out


def combination(x, n=2, return_index=True):
    """Generate the combination for a list of items
    
    Parameters
    ----------
    x : list 
        A list of items

    example
    >>> x = ['ctl', 'ctl', 'exp', 'exp'] # keep order
    >>> combination(x, return_index=True)
    [[0, 1], [2, 3]]

    >>> combination(x, return_index=False)
    ['ctl', 'exp']
    """
    xu = unique_list(x, sorted=False) # keep order
    out = []
    if len(xu) >= n:
        item_pairs = list(combinations(xu, n))
        # for index
        index_pairs = []
        for (a, b) in item_pairs:
            index_a = [i for i, j in enumerate(x) if j == a]
            index_b = [i for i, j in enumerate(x) if j == b]
            index_pairs.append([index_a, index_b])
        out = index_pairs if return_index else item_pairs
    return out


class Config(object):
    """Working with config, in dict/yaml/toml/pickle formats
    load/dump
    
    Example:
    1. write to file
    >>> Config(d).dump('out.json')
    >>> Config(d).dump('out.toml')
    >>> Config(d).dump('out.pickle')
    
    2. load from file
    >>> d = Config().load('in.yaml')

    read/write data
    """
    def __init__(self, x=None, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.x = x


    def load(self, x=None):
        """Read data from x, auto-recognize the file-type
        toml
        json
        pickle
        txt
        ...
        """
        if x == None:
            x = self.x # dict or str

        if x is None:
            x_dict = None # {} ?
        elif isinstance(x, dict):
            x_dict = dict(sorted(x.items(), key=lambda i:i[0]))
        elif isinstance(x, str):
            reader = self.get_reader(x)
            if reader is None:
                x_dict = None
                log.error('unknown x, {}'.format(x))
            else:
                x_dict = reader(x)
        else:
            x_dict = None
            log.warning('dump(x=) dict,str expect, got {}'.format(
                type(x).__name__))

        return x_dict


    def dump(self, d=None, x=None):
        """Write data to file x, auto-recognize the file-type
        d str or dict, data
        x str file to save data(dict)

        toml
        json
        pickle
        txt
        ...
        """
        if d is None:
            d = self.load(self.x)
        # make sure: dict
        if isinstance(x, str):
            writer = self.get_writer(x)
            if writer is None:
                log.error('unknown x, {}'.format(x))
            else:
                writer(d, x)
        else:
            log.warning('dump(x=) expect str, got {}'.format(
                type(x).__name__))


    def guess_format(self, x):
        """Guess the file format, by file extension
    
        file format:
        - toml
        - yaml
        - json
        - pickle

        data format:
        - dict
        """
        formats = {
            'json': 'json',
            'yaml': 'yaml',
            'yml': "yaml",
            'toml': 'toml',
            'pickle': 'pickle'
        }

        if isinstance(x, str):
            x_ext = os.path.splitext(x)[1]
            x_ext = x_ext.lstrip('.').lower()
            x_format = formats.get(x_ext, None)
        elif isinstance(x, dict):
            x_format = 'dict'
        else:
            x_format = None

        return x_format


    def get_reader(self, x):
        """Get the reader for file x, based on the file extension
    
        could be: json/yaml/toml/pickle
        """
        x_format = self.guess_format(x)
        readers = {
            'json': self.from_json,
            'yaml': self.from_yaml,
            'toml': self.from_toml,
            'pickle': self.from_pickle
        }
        return readers.get(x_format, None)


    def get_writer(self, x):
        """Get the reader for file x, based on the file extension

        could be: json/yaml/toml/pickle
        """
        x_format = self.guess_format(x)
        writers = {
            'json': self.to_json,
            'yaml': self.to_yaml,
            'toml': self.to_toml,
            'pickle': self.to_pickle
        }
        return writers.get(x_format, None)


    def from_json(self, x):
        """Loding data from JSON file
        x should be file
        """
        d = None
        if os.path.exists(x):
            try:
                with open(x, 'r') as r:
                    if os.path.getsize(x) > 0:
                        d = json.load(r) # sorted by key
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_json() failed, {}'.format(exc))
            finally:
                return d
        else:
            log.error('from_json() failed, file not exists: {}'.format(x))


    def from_yaml(self, x):
        """Loding data from YAML file
        x should be file
        """
        d = None
        if os.path.exists(x):
            try:
                with open(x, 'r') as r:
                    if os.path.getsize(x) > 0:
                        d = yaml.load(r, Loader=yaml.FullLoader)
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_yaml() failed, {}'.format(exc))
            finally:
                return d
        else:
            log.error('from_yaml() failed, file not exists: {}'.format(x))
        # with open(x, 'r') as r:
        #     try:
        #         d = yaml.safe_load(r)
        #         return collections.OrderedDict(sorted(d.items()))
        #     except yaml.YAMLError as exc:
        #         log.warning(exc)


    def from_toml(self, x):
        """Loding data from TOML file
        x should be file
        """
        d = None
        if os.path.exists(x):
            try:
                with open(x, 'r') as r:
                    if os.path.getsize(x) > 0:
                        d = toml.load(x)
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_toml() failed, {}'.format(exc))
            finally:
                return d
        else:
            log.error('from_toml() failed, file not exists: {}'.format(x))


    def from_pickle(self, x):
        """Loding data from pickle file
        x should be file
        """
        d = None
        if os.path.exists(x):
            try:
                with open(x, 'rb') as r:
                    if os.path.getsize(x) > 0:
                        d = pickle.load(r)
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_pickle() failed, {}'.format(exc))
            finally:
                return d
        else:
            log.error('from_pickle() failed, file not exists: {}'.format(x))


    def to_json(self, d, x):
        """Writing data to JSON file
        d dict, data to file
        x None or str, path to JSON file, or return string
        """
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_json(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_json(d=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_json(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wt') as w:
                    json.dump(d, w, indent=4, sort_keys=True)
                # return x
            except Exception as exc:
                log.error('to_json() failed, {}'.format(exc))


    def to_yaml(self, d, x):
        """Writing data to YAML file
        d dict, data to file
        x str, path to YAML file

        yaml.dump(), does not support OrderedDict
        Solution: OrderedDict -> json -> dict
        """
        # check
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_yaml(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_yaml(d=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_yaml(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wt') as w:
                    yaml.dump(d, w)
            except:
                log.warning('saving as YOML failed, use TOML instead')
                x_toml = os.path.splitext(x)[0] + '.toml'
                with open(x_toml, 'wt') as w:
                    toml.dump(d, w)                
#             except Exception as exc:
#                 log.error('to_yaml() failed, {}'.format(exc))


    def to_toml(self, d, x):
        """Writing data to TOML file
        d dict, data to file
        x str, path to TOML file
        """        
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_toml(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_toml(d=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_toml(d=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wt') as w:
                    toml.dump(d, w)
                # return x
            except Exception as exc:
                log.error('to_toml() failed, {}'.format(exc))


    def to_pickle(self, d, x):
        """Writing data to pickle file
        d dict, data to file
        x str, path to pickle file
        """        
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_pickle(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_pickle(x=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_pickle(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wb') as w:
                    pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)
                # return x
            except Exception as exc:
                log.error('to_pickle() failed, {}'.format(exc))


    def to_log(self, d, x, stdout=False):
        """Writing data to log file: key: value format
        d dict, data to file
        x str, path to pickle file
        """
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_log(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_log(x=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_log(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                # organize msg
                msg = []
                for k, v in d.items():
                    if isinstance(v, str) or isinstance(v, numbers.Number) or isinstance(v, bool):
                        v = str(v)
                    elif isinstance(v, list):
                        v = ', '.join(map(str, v))
                    else:
                        v = '...' # skip
                    msg.append('{:30s} | {:<40s}'.format(k, v))
                # save
                with open(x, 'wt') as w:
                    w.write('\n'.join(msg) + '\n')
                if stdout:
                    print('\n'.join(msg))
                # return x
            except Exception as exc:
                log.error('to_log() failed, {}'.format(exc))


    def _tmp(self, suffix='.txt'):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
            delete=False)
        return tmp.name
    
    
    


def convert_image(x, out_fmt='PNG'):
    if not out_fmt in ['PNG', 'JPEG', "TIFF"]:
        log.error('out_fmt: [PNG|JPEG|TIFF], {} got'.format(out_fmt))
    out_ext = out_fmt.lower()
    out_img = os.path.splitext(x)[0] + '.' + out_ext
    # read/write
    if os.path.exists(out_img):
        log.warning('file exists, skipping ...: {}'.format(out_img))
    else:
        img = Image.open(x)
        img.save(out_img, out_fmt)

        
