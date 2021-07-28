#!/usr/bin/env python3

"""General functions for file manipulation

1. file_exists: exists 
2. file_copy: copy
3. file_symlink: symlink
4. file_remove: remove 
5. file_type: type
6. ...

##############
## re-write ##
check_fx()
check_fx_paired()


##############
## re-write ##
def check_path( ):   -> check_path(x, **kwargs)
def path_remove( ):  -> remove_path(x, **kwargs)
def check_file( ):   -> check_file(x, **kwargs)
def file_remove( ):  -> remove_file(x, **kwargs)
def file_copy( ):    -> copy_file(src, dest, **kwargs)
def file_symlink( ): -> symlink_file(src, dest, **kwargs)
def file_prefix( ):  -> file_prefix(x)
def file_abspath( ): -> file_abspath(x)
def file_exists( ):  -> file_exists(x)
def listdir( ):      -> list_dir(x)
def listfile( ):     -> list_file(x, pattern='*')
def list_fx( ):      -> list_fx(x, pattern='*')
def is_gz( ):        -> file_is_gzippped()
def file_row_counter( ): -> file_nrows(x)
def fq_name( ):      -> fx_name(x, fix_pe=False)
def fq_name_rmrep(): -> fx_name(x, fix_rep)
def gzip_cmd( ):     -> gzip_file(src, dest, **kwargs)
def list_uniquer( ):
def update_obj( ):
def run_shell_cmd(cmd):
def design_combinations( ): -> string.combination()
def sam_flag_check( ):      -> bam.is_sam_flag()
def bed2gtf( ):             -> bed.bed2gtf()
def featureCounts_reader( ):-> read_featurecounts()


################
## deprecated ##
def index_checker( ):   -> del
def check_fq( ):        -> check_fx()
def fq_paired( ):       -> check_fx_paired()
def dict_to_log( ):     -> Config()
def dict_to_pickle( ):  -> Config()
def pickle_to_dict(x):  -> Config()
def listfiles( ):       -> list_file()
def listfiles2( ):      -> list_file()
def is_path( ):         -> list_dir()
def list_fq_files( ):   -> list_fx()
def in_dict( ):         -> del
def in_attr( ):         -> del
def args_checker( ):    -> Config()
def args_logger( ):     -> Config()
def symlink( ):         -> symlink_file()
def merge_names(x):
"""

import os
import sys
import re
import pathlib
import shutil
import logging
import fnmatch
import binascii
import subprocess
import pyfastx
import pysam
from xopen import xopen
from hiseq.utils.seq import Fastx


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')



##########################
## check files          ##
##########################
def check_fx(fx, **kwargs):
    """
    Check the fastq/a files
    1. file exist
    2. fq1 required
    3. fq type: fasta/q
    
    Keyword Parameters
    ------------------
    check_empty : bool
        Check the file is empty or not
        
    show_error : bool
        Display the error message
    """
    show_error = kwargs.get('show_error', False)
    out = False
    if isinstance(fx, str):
        if check_file(fx, **kwargs):
            try:
                fx_type = Fastx(fx).format # fasta/q
                out = fx_type in ['fasta', 'fastq']
            except ValueError as err:
                log.info('Failed to read file, with error: {}'.format(err))
        else:
            if show_error:
                log.error('fx failed, {}'.format(fx))
    elif isinstance(fx, list):
        out = [check_fx(i) for i in fx]
    else:
        pass
    return out


def check_fx_paired(fq1, fq2, **kwargs):
    """
    Check the fq1 and fq2 are paired or not
    the name of fq1, fq2 are the same

    fq1: @the-name/1
    fq2: @the-name/2

    scan the first read from the two files
    
    Keyword Parameters
    ------------------
    check_empty : bool
        Check the file is empty or not
        
    show_error : bool
        Display the error message
    """
    if isinstance(fq1, str) and isinstance(fq2, str):
        # name
        chk1 = fx_name(fq1, fix_pe=True) == fx_name(fq2, fix_pe=True)
        # exists
        chk2 = all(check_fx([fq1, fq2], **kwargs))
        # paired
        try:
            fx1 = pyfastx.Fastx(fq1)
            fx2 = pyfastx.Fastx(fq2)
            for a,b in zip(fx1, fx2):
                chk3 = a[0][:-1] == b[0][:-1]
                break
        except:
            chk3 = False
        # output
        out = all([chk1, chk2, chk3])
    elif isinstance(fq1, list) and isinstance(fq2, list):
        if len(fq1) == len(fq2) and len(fq1) > 0:
            out = all([check_fx_paired(f1, f2, **kwargs) for f1,f2 in zip(fq1, fq2)])
        else:
            out = False
    else:
        log.error('illegal fq1,fq2; str,list expect, got {}, {}'.format(
            type(fq1).__name__, type(fq2).__name__))
        out = False
    return out


def check_fx_args(fq1, fq2=None, **kwargs):
    """
    Check the fastx, both str or list; fq2 could be None.
    
    Parameters
    ----------
    fq1 : str or list
        read1 of PE reads, or SE
        
    fq2 : None, str or list
        read2 of PE reads
        
    check:
    1. file exists
    2. file type
    3. fq paired
    4. check_empty
    """
    if isinstance(fq1, str):
        fq1 = [fq1]
    if isinstance(fq2, str):
        fq2 = [fq2]
    if not isinstance(fq1, list):
        log.error('fq1 expect str or list, got {}'.format(
            type(fq1).__name__))
        return None
    # check fq1: message
    c1 = isinstance(fq1, list)
    c1e = all(file_exists(fq1))
    c1x = all([c1, c1e])
    # check fq2:
    c2 = isinstance(fq2, list)
    if c2:
        c2e = all(file_exists(fq2))
        c2p = check_fx_paired(fq1, fq2)
        c2x = all([c2, c2e, c2p])
    elif fq2 is None:
        c2e = c2p = False
        c2x = True # skipped
    else:
        c2x = c2e = c2p = False # force
    # final
    out = all([c1x, c2x])
    if not out:
        msg = '\n'.join([
            '='*80,
            'Check fastq:',
            '{:>14} : {}'.format('fq1', fq1),
            '{:>14} : {}'.format('fq2', fq2),
            '-'*40,
            'Status',
            '{:>14} : {}'.format('fq1 is list', c1),
            '{:>14} : {}'.format('fq1 exists', c1e),
            '{:>14} : {}'.format('fq2 is list', c2),
            '{:>14} : {}'.format('fq2 is exists', c2e),
            '{:>14} : {}'.format('fq is paired', c2p),
            '-'*40,
            'Status: {}'.format(out),
            '='*80                
        ])
        print(msg)
    return out

        
##########################
## manipulate files     ##
##########################
def check_path(x, **kwargs):
    """Check if x is path
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
    
    create_dirs : bool
        Create the dirs
    """
    show_error = kwargs.get('show_error', False)
    show_log = kwargs.get('show_log', False)
    create_dirs = kwargs.get('create_dirs', True) # default: True
    if isinstance(x, str):
        out = False
        if os.path.isdir(x):
            out = True
        elif os.path.isfile(x):
            if show_error:
                log.error('file exists, not a directory: {}'.format(x))
        else:
            if create_dirs:
                try:
                    os.makedirs(x)
                    out = True
                except:
                    if show_error:
                        log.error('`os.makedirs` failed: {}'.format(x))
        # show log
        flag = 'ok' if out else 'failed'
        if show_log is True:
            log.info('{:<6s} : {}'.format(flag, x))
    elif isinstance(x, list):
        out = all([check_path(i, **kwargs) for i in x])
    else:
        if show_error:
            log.error('x expect str or list, got {}'.format(type(x).__name__))
    return out


def remove_path(x, **kwargs):
    """Remove directory
    
    Parameters
    ----------
    x : str or list
        The file(s) to be removed 

    ask : bool
        Ask the user, before the files removed.

    check_empty : True
        Do not proceed, if directory is not empty
        
    show_log : True
        Display the status of the files
        
    show_error : False
        Display the error messages
    """
    ask = kwargs.get('ask', True)
    check_empty = kwargs.get('check_empty', True)
    show_log = kwargs.get('show_log', True)
    show_error = kwargs.get('show_error', False)
    if isinstance(x, str):
        if os.path.isdir(x):
            x_files = list_dir(x, full_name=True, recursive=False, include_dir=True)
            is_empty = len(x_files) == 0
            if check_empty and not is_empty:
                is_rm = False                
            else:
                is_rm = True
            # rm, ask
            if is_rm:
                ask_msg = input('Remove: {}, [Y|n]: '.format(x)) if ask else 'Y'
            else:
                ask_msg = 'no'
            rm_tag = 'yes' if ask_msg.lower() in ['y', 'yes'] else 'no'
            empty_tag = 'yes' if is_empty else 'no'
            # do-the-thing, removing
            if rm_tag == 'yes':
                try:
                    shutil.rmtree(x)
                except:
                    rm_tag = 'no'
                    if show_error:
                        log.error('failed, remove path: {}'.format(x))
        else:
            rm_tag = 'no'
            empty_tag = 'NA'
            if show_error:
                log.error('x is not path, {}'.format(x))
        if show_log:
            log.info('rm:{:3s}\tis_empty:{}\t{:3s}'.format(rm_tag, empty_tag, x))    
    elif isinstance(x, list):
        [remove_path(i, **kwargs) for i in x]
    elif isinstance(x, dict):
        for k,v in x.items:
            if isinstance(v, str) or isinstance(v, list):
                remove_path(v, **kwargs)
    else:
        log.error('x, str or list or dict expected, got {}'.format(
            type(x).__name__))


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
    if isinstance(x, str):
        if file_exists(x):
            x_size = os.stat(x).st_size
            # empty gzipped file, size=20
            q_size = 20 if x.endswith('.gz') else 0
            out = x_size > q_size if check_empty else True
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
        if show_error:
            log.error('x expect str or list, got {}'.format(type(x).__name__))
        out = False
    return out


def symlink_file(src, dest, absolute_path=False, force=False):
    """Create symlink, 
    
    Parameters
    ----------
    src : str
        The source file
        
    dest : str
        The destinate file, dir exists
        
    absolute_path : bool
        Use abs_path instead

    force : bool
        Copy files, overwrite dest file
    """
    if not isinstance(src, str):
        log.error('src, expect str, got {}'.format(type(src).__name__))
    elif not isinstance(dest, str):
        log.error('dest, expect str, got {}'.format(type(src).__name__))
    elif os.path.isfile(src):
        src = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        if os.path.isdir(dest):
            dest_file = os.path.join(dest, os.path.basename(src))
        else:
            dest_file = dest
        dest_file = os.path.abspath(os.path.expanduser(os.path.expandvars(dest_file)))
        # the relative path of src
        src_dir = os.path.dirname(src)
        src_name = os.path.basename(src)
        dest_dir = os.path.dirname(dest_file)
        src_dir_rel = os.path.relpath(src_dir, dest_dir)
        src_rel = os.path.join(src_dir_rel, src_name)
        src_file = src if absolute_path else src_rel
        # do-the-thing
        if file_exists(dest_file) and not force:
            log.info('symlink_file() skipped, dest exists: {}'.format(dest_file))
        else:
            try:
                os.symlink(src_file, dest_file)
            except:
                log.error('symlink_file() failed, {}'.format(dest_file))
    elif os.path.islink(src):
        pass
    else:
        log.warning('symlink_file() failed, src not vaild: {}'.format(src))


def copy_file(src, dest, force=False):
    """Copy file
    
    Parameters
    ----------
    x : str or list
        The file(s) to be removed 
        
    force : bool
        Copy files, overwrite dest file
    """
    if not isinstance(src, str):
        log.error('src, expect str, got {}'.format(type(src).__name__))
    elif not isinstance(dest, str):
        log.error('dest, expect str, got {}'.format(type(src).__name__))
    elif os.path.isfile(src):
        if os.path.isdir(dest):
            dest_file = os.path.join(dest, os.path.basename(src))
        else:
            dest_file = dest
        # do-the-thing
        if file_exists(dest_file) and not force:
            log.error('copy_file() skipped, dest exists: {}'.format(dest_file))
        else:
            try:
                shutil.copy(src, dest_file)
            except:
                log.error('copy_file() failed, {}'.format(dest_file))
    else:
        log.warning('copy_file() failed, src not vaild: {}'.format(src))


def remove_file(x, **kwargs):
    """Remove files
    
    Parameters
    ----------
    x : str or list
        The file(s) to be removed 
        
    ask : bool
        Ask the user, before the files removed.
    """
    ask = kwargs.get('ask', True)
    show_log = kwargs.get('show_log', True)
    show_error = kwargs.get('show_error', False)
    if isinstance(x, str):
        if os.path.isfile(x):
            file_tag = 'yes'
            ask_msg = input('Remove: {}, [Y|n]: '.format(x)) if ask else 'Y'
        else:
            file_tag = 'no'
            ask_msg = 'no'
        rm_tag = 'yes' if ask_msg.lower() in ['y', 'yes'] else 'no'
        # do-the-thing, removing
        if rm_tag == 'yes':
            try:
                os.remove(x)
            except:
                rm_tag = 'no'
                if show_error:
                    log.error('failed, remove file: {}'.format(x))
        if show_log:
            log.info('rm:{:3s}\tis_file:{}\t{:3s}'.format(rm_tag, file_tag, x))
    elif isinstance(x, list):
        [remove_file(i, **kwargs) for i in x]
    elif isinstance(x, dict):
        for k,v in x.items:
            if isinstance(v, str) or isinstance(v, list):
                remove_file(v, **kwargs)
    else:
        log.error('x, str or list or dict expected, got {}'.format(
            type(x).__name__))


def gzip_file(src, dest=None, decompress=True, **kwargs):
    """Gzip Compress or Decompress files using gzip module in python
    rm, True/False, whether remove old file

    # check the src file by extension: .gz
    """
    a = os.path.exists(src)
    b = file_exists(src)
    show_log = kwargs.get('show_log', True)
    show_error = kwargs.get('show_error', False)
    compresslevel = kwargs.get('compresslevel', 1)
    threads = kwargs.get('threads', 4)
    flag = False
    # input: src
    if not isinstance(src, str):
        if show_error:
            log.error('src expect str, got {}'.format(src))
    if not file_exists(src):
        if show_error:
            log.error('src not exists, {}'.format(src))
    # output: dest
    if dest is None:
        dest = os.path.splitext(src)[0] if decompress else src + '.gz'
    if isinstance(dest, str):
        if file_exists(dest):
            if show_error:
                log.error('dest exists, {}'.format(dest))
        elif os.path.exists(os.path.dirname(dest)):
            flag = True
        else:
            if show_error:
                log.error('dest not valid, {}'.format(dest))
    else:
        if show_error:
            log.error('dest expect str, got {}'.format(dest))
    # do-the-thing
    out = None
    if flag:
        src = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        dest = os.path.abspath(os.path.expanduser(os.path.expandvars(dest)))
        print('!A-1', dest)
        is_gzipped = file_is_gzipped(src)
        if decompress:
            if is_gzipped:
                with xopen(src, 'rb') as r, \
                    xopen(dest, 'wb') as w:
                    shutil.copyfileobj(r, w)
                out = dest
            else:
                if show_error:
                    log.error('src is not gzipped, skipped. {}'.format(src))
        else:
            if is_gzipped:
                if show_error:
                    log.error('src is gzipped, skipped. {}'.format(src))
            else:
                with xopen(src, 'rb') as r, \
                    xopen(dest, 'wb', threads=threads,
                          compresslevel=compresslevel) as w:
                    shutil.copyfileobj(r, w)
                out = dest
    return out


#######################
## 1. search files   ##
#######################
def list_dir(x, full_name=True, recursive=False, include_dir=False):
    """List all the files within the path
    see: list.dir() in R
    see answers on :https://stackoverflow.com/a/3207973
    
    Parameters
    ----------
    x : str
        List files (dirs) in x
        
    full_name : bool
        Return the fullname of the files/dirs 
        
    recursive : bool
        List files/dirs recursively 
        
    include_dir : bool
        Return the dirs
    """
    out = []
    if isinstance(x, str):
        if os.path.isdir(x):
            n = 0 # levels
            for (root, d, f) in os.walk(x):
                dirs = [os.path.join(root, i) for i in d] if full_name else d
                files = [os.path.join(root, i) for i in f] if full_name else f
                out += files
                if include_dir:
                    out += dirs
                if not recursive:
                    break # first level
        else:
            log.error('list_dir() skipped, x not a directory: {}'.format(x))
    else:
        log.error('list_dir() skipped, x expect str, got {}'.format(
            type(x).__name__))
    return sorted(out)


def list_file(path='.', pattern='*', full_name=True, recursive=False,
    include_dir=False):
    """Search files by the pattern, within directory
    fnmatch.fnmatch()
    
    Parameters
    ----------
    x : str
        List files (dirs) in x
        
    pattern : str
        see pattern of fnmatch
        pattern:
        *       matches everything
        ?       matches any single character
        [seq]   matches any character in seq
        [!seq]  matches any char not in seq

        An initial period in FILENAME is not special.
        Both FILENAME and PATTERN are first case-normalized
        if the operating system requires it.
        If you don't want this, use fnmatchcase(FILENAME, PATTERN).
        
    full_name : bool
        Return the fullname of the files/dirs 
        
    recursive : bool
        List files/dirs recursively 
        
    include_dir : bool
        Return the dirs

    example:
    list_file('./', '*.fq')
    """
    files = list_dir(path, full_name, recursive, include_dir)
    files = [f for f in files if fnmatch.fnmatch(os.path.basename(f), pattern)]
    return sorted(files)


def list_fx(x, recursive=False):
    """List the fasta/q files, all
    *.f[aq]
    *.f[aq].gz
    *.fast[aq]
    *.fast[aq].gz
    
    Parameters
    ----------
    x : str
        Path to the directory
        
    recursive : bool
        List files recursively (be carefull, too-much files)
    """
    ext_list = ['.fa', '.fq', '.fasta', '.fastq']
    ext_list += [i + '.gz' for i in ext_list]
    out = list_dir(x, full_name=True, include_dir=False, recursive=recursive)
    if len(out) > 0:
        # out = [i for i in out if os.path.splitext(i)[1].lower() in ext_list]
        p = re.compile('\.f(ast)?(a|q)(\.gz)?', flags=re.IGNORECASE)
        out = [i for i in out if p.search(i)]
    return sorted(list(set(out)))


def list_fx2(x, pattern='*', recursive=False):
    """List the fasta/q files, by name
    *.f[aq]
    *.f[aq].gz
    *.fast[aq]
    *.fast[aq].gz
    
    Parameters
    ----------
    x : str
        Path to the directory
        
    pattern : str
        The pattern of the file
        pattern:
        *       matches everything
        ?       matches any single character
        [seq]   matches any character in seq
        [!seq]  matches any char not in seq

    recursive : bool
        List files recursively (be carefull, too-much files)
    """
    fx = list_fx(x, recursive)
    fx = [f for f in fx if fnmatch.fnmatch(os.path.basename(f), pattern)]
    return sorted(fx)

        
##########################
## information of files ##
##########################
def file_is_gzipped(x):
    """Check if the file is gzipped or not
    see answer: https://stackoverflow.com/a/3703300
    and on wiki: https://en.wikipedia.org/wiki/Gzip 
    see also this blog: https://www.thinbug.com/q/3703276
    
    check the magic number for gzipped file: '1f8b'
    """
    if file_exists(x):
        with open(x, 'rb') as r:
            out = binascii.hexlify(r.read(2)) == b'1f8b'
    else:
        out = False
    return out


def file_prefix(x, with_path=False):
    """Extract the prefix of file, 
    compatiabe for None
    
    Parameters
    ----------
    x : str,list
        Path to a file, or list of files
    
    remove extensions
    .gz, .fq.gz
    """
    if isinstance(x, str):
        if x.endswith('.gz') or x.endswith('.bz2'):
            x = os.path.splitext(x)[0]
        out = os.path.splitext(x)[0]
        if not with_path:
            out = os.path.basename(out)
    elif isinstance(x, list):
        out = [file_prefix(i, with_path) for i in x]
    elif x is None:
        out = None
    else:
        log.error('unknown x, str,list,None expected, got {}'.format(
            type(x).__name__))
        out = None
    return out


def file_abspath(x):
    """Return the absolute path of file
    
    Parameters
    ----------
    x : str,list
        Path to a file, or list of files
    """
    if x is None or x == 'None': # in case toml format?!
        out = None
    elif isinstance(x, str):
        out = os.path.abspath(os.path.expanduser(x))
    elif isinstance(x, list):
        out = [file_abspath(i) for i in x]
    else:
        log.warning('x, expect str,list, got {}'.format(type(x).__name__))
        out = x
    return out


def file_exists(x):
    """Check if file exists or not
    
    Parameters
    ----------
    x : str,list
        Path to a file, or list of files
    """
    if x is None:
        out = False
    elif isinstance(x, str):
        out = os.path.exists(x) # file/dir/link
#         out = os.path.isfile(x)
    elif isinstance(x, list):
        out = [file_exists(i) for i in x]
    else:
        log.warning('x, expect str,list, got {}'.format(type(file).__name__))
        out = False
    return out

        
def file_nrows(x):
    """Count the file rows
    count '\n' 
    from @glglgl on stackoverflow, modified
    https://stackoverflow.com/a/9631635/2530783
    """
    def blocks(files, size = 1024 * 1024):
        while True:
            b = files.read(size)
            if not b: break
            yield b
    if file_exists(x):
        with xopen(x, 'rt') as r:
            out = sum(bl.count('\n') for bl in blocks(r))
    else:
        log.error('x, file not exists, {}'.format(x))
        out = None
    return out


def fx_name(x, fix_pe=False, fix_rep=False, fix_unmap=False):
    """The name of fastx
    fix the pe_suffix, '_1, _2', '_R1, _R2'
    
    Parameters
    ----------
    x : str or list
        Path to the fastx files
        
    fix_pe : bool
        Remove the suffix of Paired-end files, '_1', '_R1'
        
    fix_rep : bool
        Remove suffix for replicate, '_rep1, _rep2'
        
    fix_unmap : bool
        Remove suffix for unmap, '.unmap'
    """
    if isinstance(x, str):
        out = file_prefix(x)
        if fix_pe:
            out = re.sub('[._](r)?[12]$', '', out, flags=re.IGNORECASE)
        if fix_unmap:
            out = re.sub('.unmap$', '', out, flags=re.IGNORECASE)
        if fix_rep:
            out = re.sub('[._](rep|r)[0-9]+$', '', out, flags=re.IGNORECASE)
    elif isinstance(x, list):
        out = [fx_name(i, fix_pe, fix_rep, fix_unmap) for i in x]
    else:
        out = None
    return out


def read_lines(x, nrows=0, skip=0, strip_white=True, comment=''):
    """
    Read plain text file
    save each line as list()
    
    Parameters
    ----------
    x:  str
        Path to a file
        
    nrows:  int
        The maximum number of rows to read 
        default: [0], ignored
        
    skip:  int
        The number of lines of the data to skip before beginning to read data
        default: [0]
        
    strip_white:  bool
        Stipping of leading and trailing white space, default: [True]
    
    comment:  str
        A string of one character, default [''], empty
    
    """
    out = None
    if isinstance(x, str):
        if os.path.exists(x):
            if os.path.isdir(x):
                log.error('read_lines() failed, exptect <file>, got <dir>')
            else:
                l = []
                try:
                    i = 0
                    with open(x) as r:
                        for line in r:
                            i += 1
                            # skip rows
                            if skip > 0 and i <= skip:
                                continue
                            # white-spaces
                            s = line.strip()
                            if strip_white:
                                s = s.lstrip()
                            # comment
                            if len(comment) == 1 and s.startswith(comment):
                                continue
                            # nrows
                            if nrows > 0 and i > nrows:
                                break
                            # save to output
                            l.append(s)
                except IOError as e:
                    log.error(e)
                out = l
        else:
            log.error('read_lines() failed, file not exists')
    else:
        log.error('read_lins() failed, expect str, got {}'.format(
            type(x).__name__))
    return out


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
        if genome_path is None:
            genome_path = os.path.join(str(pathlib.Path.home()),
                'data', 'genome')
        self.genome_path = genome_path


    def get_fa(self):
        """
        Get the fasta file of specific genome
        {genome}/bigZips/{genome}.fa
        also check ".gz" file
        """
        fa = os.path.join(self.genome_path, self.genome, 'bigZips',
            self.genome + '.fa')
        if not file_exists(fa):
            # gencode version
            fa = os.path.join(self.genome_path, self.genome, 'fasta',
                self.genome + '.fa')
        fa_gz = fa + '.gz'
        if not file_exists(fa):
            if file_exists(fa_gz):
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
        if not file_exists(fa_size):
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
        if not file_exists(p):
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
        if not file_exists(g):
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
        if not file_exists(gtf):
            gtf = os.path.join(
            self.genome_path,
            self.genome,
            'gtf',
            self.genome + '.' + version + '.gtf')
        if not file_exists(gtf):
            gtf = None
        return gtf


    def te(self, format='gtf'):
        """Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        # only dm3 supported
        te_gtf = os.path.join(self.genome_path, self.genome,
            self.genome + '_transposon',
            self.genome + '_transposon.' + format)
        if not file_exists(te_gtf):
            te_gtf = None
        return te_gtf


    def piRNA_cluster(self, format='gtf'):
        """Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        # only dm3 supported
        te_gtf = os.path.join(self.genome_path, self.genome,
            self.genome + '_piRNA_clusters',
            self.genome + '_piRNA_clusters.' + format)
        if not file_exists(te_gtf):
            te_gtf = None
        return te_gtf
    
    
    
def bed_to_gtf(file_in, file_out):
    """Convert BED to GTF
    chrom chromStart chromEnd name score strand
    """
    with open(file_in) as r, open(file_out, 'wt') as w:
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
    return file_out
