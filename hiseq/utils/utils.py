#!/usr/bin/env python3

"""Functions for string, list, dict, ...
transform to another version 
attribution of the vectors
"""



import os
import sys
import logging
from itertools import combinations


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


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


