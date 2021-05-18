#!/usr/bin/env python3

"""Functions for string, list, dict, ...
transform to another version 
attribution of the vectors
"""



import os
import sys
import json
import json
import yaml
import toml
import pickle
import logging
import collections
import subprocess
from datetime import datetime
from dateutil import tz
from itertools import combinations
# from .helper import log
# from .file import file_exists, file_abspath


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')



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
            x_dict = collections.OrderedDict(sorted(x.items()))
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
                        d = json.load(r)
                        d = collections.OrderedDict(sorted(d.items()))
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
                        d = collections.OrderedDict(sorted(d.items()))
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
                        d = collections.OrderedDict(sorted(d.items()))
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
                        d = collections.OrderedDict(sorted(d.items()))
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
        x_yaml = x
        x = os.path.splitext(x_yaml)[0] + '.toml'
        log.warning('OrderedDict is not supported in YAML, save as TOML instead: {}'.format(x))
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
                    toml.dump(d, w)
                # return x
            except Exception as exc:
                log.error('to_yaml() failed, {}'.format(exc))


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
    