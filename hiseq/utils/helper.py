
# -*- coding: utf-8 -*-

"""
Common functions for pipeline construction
file modification, 
...
"""

import os
import sys
import gzip
import shutil
import json
import pickle
import fnmatch
import logging
import functools
import subprocess


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)


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
                self.logger.info('{0} - {1} - {2}'.format(fn.__name__, args, kwargs))
                result = fn(*args, **kwargs)
                self.logger.info(result)
                # return result
            except Exception as ex:
                self.logger.info('Exception {0}'.format(ex))
                raise ex
            return result
        return decorated


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
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(pid, pgid, rc,
        stderr.strip(), stdout.strip())
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


def file_prefix(fn, with_path = False):
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


def gzip_cmd(f, t, decompress=True, rm=True):
    """
    Gzip Compress or Decompress files using gzip module in python 
    rm, True/False, whether remove old file

    # check the f file by extension: .gz
    """
    if decompress:
        if f.endswith('.gz'):
            raise Exception('not a gzipped file: {}'.format(f))
        with gzip.open(f, 'rb') as r, open(t, 'wb') as w:
            shutil.copyfileobj(r, w)
    else:
        with open(f, 'rb') as r, gzip.open(t, 'wb') as w:
            shutil.copyfileobj(r, w)

    # output
    if rm is True:
        os.remove(f)

    return t


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


class Json(object):

    def __init__(self, x):
        """
        x 
          - dict, save to file
          - json, save to file
          - file, read as dict
        Save dict to json file
        Read from json file as dict
        ...
        """
        self.x = x # input

        if isinstance(x, Json):
            self.dict = x.dict
        elif isinstance(x, dict):
            # input a dict, 
            # save to file
            self.dict = x
        elif os.path.exists(x):
            # a file saving json content
            self.dict = self.reader()
        else:
            raise Exception('unknown objec: {}'.format(x))

    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json', delete=False)
        return tmp.name


    def reader(self):
        """
        Read json file as dict
        """
        if os.path.getsize(self.x) > 0:
            with open(self.x) as r:
                d = json.load(r)
        else:
            d = {}

        return d


    def writer(self, f=None):
        """
        Write d (dict) to file x, in json format
        """
        # save to file
        if f is None:
            f = self._tmp()

        assert isinstance(f, str)
        # assert os.path.isfile(f)

        if isinstance(self.x, dict):
            with open(f, 'wt') as w:
                json.dump(self.x, w, indent=4, sort_keys=True)

        return f


