#!/usr/bin/env python3

"""
General functions for file manipulation

check_file : file exists, check_empty 
symlink_file : symlink file 
copy_file : str 
remove_file : str
gzip_file : str
gunzip_file : str
file_abspath : str
file_is_gzipped : str
file_prefix : str
file_exists : str
file_nrows : str
check_path : str
remove_path : str
copy_path : str
list_dir : str
list_file : str
list_fx : str, all files
list_fx2 : str, pattern
fx_name : str
check_fx : str
check_fx_paired : str
check_fx2 : str,list

read_lines : str
str_distance : str

Genome : 
bed_to_gtf :
"""

import os
from sys import stdout
from re import sub, compile, IGNORECASE
from pathlib import Path # Path.home
from shutil import copy, copytree, copyfileobj, rmtree
from logging import basicConfig, getLogger
from fnmatch import fnmatch
from binascii import hexlify
from pyfastx import Fastx
from pysam import faidx
from xopen import xopen
import Levenshtein as lev # distance


basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=stdout)
log = getLogger(__name__)
log.setLevel('INFO')

# file: only for single file, str
def check_file(x, **kwargs):
    """Check if x is file and exists
    Parameters
    ----------
    x : str
        Path to a file
    
    Parameters
    ------------------
    show_error : bool
        Show the error messages
        
    show_log : bool
        Show the log messages 
        
    check_empty : bool
        Check if the file is empty or not,  gzipped empty file, size=20
    """
    args = {
        'show_error': False,
        'show_log': False,
        'check_empty': False
    }
    args.update(kwargs)
    if isinstance(x, str):
        if os.path.isfile(x): # symlink/file
            out = True
            if args['check_empty']:
                x_size = os.stat(x).st_size
                out = x_size > 20 if x.endswith('.gz') else x_size > 0
        else:
            if args['show_error']:
                log.error('file not exists: {}'.format(x))
            out = False # failed
        # log info
        if args['show_log']:
            flag = 'ok' if out else 'failed'
            log.info('{:<6s} : {}'.format(flag, x))
    # elif isinstance(x, list):
    #     out = all([check_file(i, **kwargs) for i in x])
    else:
        if args['show_error']:
            log.error('x expect str or list, got {}'.format(type(x).__name__))
        out = False
    return out


def symlink_file(src, dest, absolute_path=False, force=False):
    """Create symlink    
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
    if not all(list(map, isinstance, [src, dest], [str, str])):
        log.error('src, dest, expect str, got src: {}, dest: {}'.format(
            type(src).__name__, type(dest).__name__))
    else:
        if check_file(src): # file/link, exists
            src_file = file_abspath(src)
            if os.path.isdir(dest):
                dest_file = os.path.join(dest, os.path.basename(src))
            else:
                dest_file = dest
            dest_file = file_abspath(dest_file)
            if not absolute_path:
                src_dir = os.path.dirname(src_file)
                dest_dir = os.path.dirname(dest_file)
                src_dir_rel = os.apth.relpath(src_dir, dest_dir)
                src_file = os.path.join(src_dir_rel, os.path.basename(src_file))
            if file_exists(dest_file) and not force:
                log.info('symlink_file() skipped, file exists: {}'.format(
                    dest_file))
            else:
                try:
                    os.symlink(src_file, dest_file)
                except:
                    log.error('symlink_file() failed, {}'.format(dest_file))
        else:
            log.error('src not exists, {}'.format(src))


def copy_file(src, dest, force=False):
    """Copy file
    
    Parameters
    ----------
    x : str
        Path to the file
        
    force : bool
        Copy files, overwrite dest file
    """
    if not all(list(map, isinstance, [src, dest], [str, str])):
        log.error('src, dest, expect str, got src: {}, dest: {}'.format(
            type(src).__name__, type(dest).__name__))
    else:
        if check_file(src): # file/link, exists
            src_file = file_abspath(src)
            if os.path.isdir(dest):
                dest_file = os.path.join(dest, os.path.basename(src))
            else:
                dest_file = dest
            if file_exists(dest_file) and not force:
                log.info('copy_file() skipped, file exists: {}'.format(
                    dest_file))
            else:
                try:
                    copy(src_file, dest_file) # shutil.copy
                except:
                    log.error('copy_file() failed, {}'.format(dest_file))
        else:
            log.error('src not exists, {}'.format(src))


def remove_file(x, **kwargs):
    """Remove files
    Parameters
    ----------
    x : str or list
        The file(s) to be removed
        
    ask : bool
        Ask the user, before the files removed.
    """
    args = {
        'ask': True,
        'show_error': False,
        'show_log': False
    }
    args.update(kwargs)
    if check_file(x):
        flag = 'skipped'
        if args['ask']:
            ask_input = input('Remove: {}, [Y|n]: '.format(x))
        else:
            ask_input = 'Y'
        if ask_input.lower() in ['y', 'yes']:
            try:
                os.remove(x)
                flag = 'removed'
            except:
                if args['show_error']:
                    log.error('remove_file() failed: {}'.format(x))
        if args['show_log']:
            log.info('{:<8s}: {}'.format(flag, x))
    else:
        log.error('x is not file, got {}'.format(type(x).__name__))


def gzip_file(src, dest=None, decompress=True, **kwargs):
    """Gzip Compress file using gzip module in python
    rm, True/False, whether remove old file
    action:
    2. compress: file -> gzipped file

    # check the src file by extension: .gz
    """
    args = {
        # 'ask': True,
        'show_error': False,
        'show_log': False,
        'compress_level': 1,
        'threads': 4
    }
    args.update(kwargs)
    if check_file(src):
        if not isinstance(dest, str):
            dest = src + '.gz'
        if check_file(dest):
            if args['show_error']:
                log.warning('dest exists, {}'.format(dest))
        else:
            src = file_abspath(src)
            dest = file_abspath(dest)
            with xopen(src, 'rb') as r, \
                xopen(dest, 'wb', threads=args['threads'],
                      compresslevel=args['compress_level']) as w:
                copyfileobj(r, w) # shutil.copyfileobj
    else:
        if args['show_error']:
            log.error('src not exists: {}'.format(x))


def gunzip_file(src, dest=None, decompress=True, **kwargs):
    """Gzip decompress file using gzip module in python
    rm, True/False, whether remove old file
    action:
    2. decompress: gzipped file -> file

    # check the src file by extension: .gz
    """
    args = {
        # 'ask': True,
        'show_error': False,
        'show_log': False,
        'compress_level': 1,
        'threads': 4
    }
    args.update(kwargs)
    if check_file(src):        
        if file_is_gzipped(x):
            if not isinstance(dest, str):
                dest = os.path.splitext(src)[0]
            if check_file(dest):
                if args['show_error']:
                    log.warning('dest exists, {}'.format(dest))
            else:
                src = file_abspath(src)
                dest = file_abspath(dest)
                with xopen(src, 'rb') as r, \
                    xopen(dest, 'wb', threads=args['threads'],
                          compresslevel=args['compress_level']) as w:
                    shutil.copyfileobj(r, w)
        else:
            if args['show_error']:
                log.error('src is not gzipped file: {}'.format(x))
    else:
        if args['show_error']:
            log.error('src not exists: {}'.format(x))


def file_abspath(x):
    """Expand the absolute path of file
    Parameters
    ----------
    x : str,list
        Path to a file
    """
    if isinstance(x, str):
        out = os.path.abspath(os.path.expanduser(os.path.expandvars(x)))
    # elif isinstance(x, list):
    #     out = [file_abspath(i) for i in x]
    else:
        log.warning('x, expect str, got {}'.format(type(x).__name__))
        out = x
    return out


def file_is_gzipped(x):
    """Check if the file is gzipped or not
    see answer: https://stackoverflow.com/a/3703300
    and on wiki: https://en.wikipedia.org/wiki/Gzip 
    see also this blog: https://www.thinbug.com/q/3703276
    
    check the magic number for gzipped file: '1f8b'
    """
    if check_file(x):
        with open(x, 'rb') as r:
            out = binascii.hexlify(r.read(2)) == b'1f8b'
    else:
        out = False
    return out


def file_prefix(x, with_path=False):
    """Extract the prefix of file, compatiabe for None
    
    Parameters
    ----------
    x : str
        Path to a file, or list of files
    
    remove extensions
    .gz, .fq.gz, tar.gz, ...
    """
    if check_file(x):
        if x.endswith('.gz') or x.endswith('.bz2'):
            x = os.path.splitext(x)[0]
        out = os.path.splitext(x)[0]
        if not with_path:
            out = os.path.basename(out)
    # elif isinstance(x, list):
    #     out = [file_prefix(i, with_path) for i in x]
    elif x is None:
        out = None
    else:
        log.warning('x, expect str, got {}'.format(type(x).__name__))
        out = None
    return out


def file_exists(x):
    """Check if file exists or not    
    Parameters
    ----------
    x : str
        Path to a file, or list of files
    """
    return check_file(x)
    # if check_file(x):
    # if x is None:
    #     out = False
    # elif isinstance(x, str):
    #     out = os.path.exists(x) # file/dir/link
    # elif isinstance(x, list):
    #     out = [file_exists(i) for i in x]
    # else:
    #     log.warning('x, expect str,list, got {}'.format(type(file).__name__))
    #     out = False
    # return out


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
    # if file_exists(x):
    if check_file(x):
        with xopen(x, 'rt') as r:
            out = sum(bl.count('\n') for bl in blocks(r))
    else:
        log.warning('file not exists: {}'.format(x))
        out = None
    return out


# path: only for single path, str
def check_path(x, **kwargs):
    """Check if x is path
    Parameters
    ----------
    x : str or list
        Path to a path
    
    Parameters
    ------------------
    show_error : bool
        Show the error messages
        
    show_log : bool
        Show the log messages 
    
    create_dir : bool
        Create the dirs
    """
    args = {
        'show_error': False,
        'show_log': False,
        'create_dir': True
    }
    args.update(kwargs)
    if isinstance(x, str):
        out = False
        if os.path.isdir(x):
            out = True
        elif os.path.isfile(x):
            if args['show_error']:
                log.error('not a directory: {}'.format(x))
        else:
            if args['create_dir']:
                try:
                    os.makedirs(x)
                    out = True
                except:
                    if args['show_error']:
                        log.error('`os.makedirs` failed: {}'.format(x))
        # show log
        flag = 'ok' if out else 'failed'
        if args['show_log'] is True:
            log.info('{:<6s} : {}'.format(flag, x))
    # elif isinstance(x, list):
    #     out = all([check_path(i, **kwargs) for i in x])
    else:
        if args['show_error']:
            log.error('x expect str, got {}'.format(type(x).__name__))
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
    args = {
        'ask': True,
        'check_empty': True,
        'show_error': False,
        'show_log': False,
        'create_dir': False
    }
    args.update(kwargs)
    if check_path(x, **args):
        flag = 'skipped'
        x_files = list_dir(x, full_name=True, recursive=False, include_dir=True)
        skipped = args['check_empty'] and len(x_files) > 0
        if args['ask']:
            ask_input = input('Remove: {}, [Y|n]: '.format(x))
        else:
            ask_input = 'Y'
        if ask_input.lower() in ['y', 'yes'] and not skipped:
            try:
                rmtree(x) # shutil.rmtree
                flag = 'removed'
            except:
                if args['show_error']:
                    log.error('remove_path() failed: {}'.format(x))
        if args['show_log']:
            log.info('{:<8s}: {}'.format(flag, x))
    else:
        log.error('x is not directory {}'.format(x))

        
def copy_path(src, dest, force=False):
    """Copy the whole directory
    
    Parameters
    ----------
    src : str
        The path to directory

    dest : str
        The path to directory
        
    force : bool
        Copy files, overwrite dest file
    """
    if not all(list(map, isinstance, [src, dest], [str, str])):
        log.error('src, dest, expect str, got src: {}, dest: {}'.format(
            type(src).__name__, type(dest).__name__))
    else:
        if check_path(src, create_dir=False):
            src_path = file_abspath(src)
            dest_path = file_abspath(dest)
            try:
                shutil.copytree(src_path, dest_path) # shutil.copytree
            except:
                log.error('copy_path() failed, {}'.format(dest_path))
        else:
            log.error('src not exists, {}'.format(src_path))
        

# search files
def list_dir(x, full_name=True, recursive=False, include_dir=False):
    """List all the files in path
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
    if check_path(x, create_dir=False):
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
        log.error('list_dir() faild, not a directory: {}'.format(x))
    return out


def list_file(x='.', pattern='*', full_name=True, recursive=False,
    include_dir=False):
    """Search files by the pattern, within directory
    fnmatch()
    
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
    file_list = list_dir(x, full_name, recursive, include_dir)
    file_list = [f for f in file_list if fnmatch(os.path.basename(f), pattern)]
    return sorted(file_list)


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
    # ext_list = ['.fa', '.fq', '.fasta', '.fastq']
    # ext_list += [i + '.gz' for i in ext_list]
    out = list_dir(x, full_name=True, include_dir=False, recursive=recursive)
    if len(out) > 0:
        p = compile('\.f(ast)?(a|q)(\.gz)?', flags=IGNORECASE)
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
    fx = [f for f in fx if fnmatch(os.path.basename(f), pattern)]
    return sorted(fx)


# fastx file: 
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
        if fix_pe: # _1, _2
            out = sub('[._](r)?[12]$', '', out, flags=IGNORECASE) # re.sub
        if fix_rep:
            out = sub('[._](rep|r)[0-9]+$', '', out, flags=IGNORECASE) # re.sub
        if fix_unmap:
            out = sub('.unmap$', '', out, flags=IGNORECASE) # re.sub
    # elif isinstance(x, list):
    #     out = [fx_name(i, fix_pe, fix_rep, fix_unmap) for i in x]
    else:
        out = None
    return out


def check_fx(x, **kwargs):
    """Check if x is fastx file
    1. file exist
    2. fq type: fasta/q, sequence in one-line mode
    
    Keyword Parameters
    ------------------
    check_empty : bool
        Check the file is empty or not
        
    show_error : bool
        Display the error message
    """
    if check_file(x, check_empty=True):
        # determined by first character: @/>
        a = read_lines(x, nrows=1)
        out = a[0] in '>@' # ['@', '>']
    else:
        out = False
    return out


def check_fx_paired(fx1, fx2, **kwargs):
    """Check if fx1 and fx2 are paired or not
    1. file name: of fx1, fx2 are the same
    2. fx name:

    fx1: @the-name/1
    fx2: @the-name/2

    scan the first read from the two files
    
    Keyword Parameters
    ------------------
    check_empty : bool
        Check the file is empty or not
        
    show_error : bool
        Display the error message
    """
    out = False
    if isinstance(fx1, str) and isinstance(fx2, str):
        if check_fx(fx1) and check_fx(fx2):
            # 1. file name
            c1 = fx_name(fx1, fix_pe=True) == fx_name(fx2, fix_pe=True)
            # 2. file name, diff=1
            c2 = str_distance(fx1, fx2, partial=False) == 1
            # 3. fxname
            f1 = read_lines(fx1, nrows=1)
            f2 = read_lines(fx2, nrows=1)
            fn1 = f1.split()[0] # name
            fn2 = f2.split()[0] # name
            c3 = fn1[:-1] == fn2[:-1]
            # the same file?!
            if fx1 == fx2:
                log.warning('fx1 and fx2 are the same file')
                out = False
            out = all([c1, c2, c3])
    # elif isinstance(fx1, list) and isinstance(fx2, list):
    #     if len(fx1) == len(fx2):
    #         out = all([check_fx_paired(i, j) for i,j in zip(fx1, fx2)])
    else:
        pass
    return out


def check_fx2(fx1, fx2=None, **kwargs): # for list, check_fx_args
    """Check the fastx, fx2 could be None.
    
    Parameters
    ----------
    fx1 : str or list
        read1 of PE reads, or SE
        
    fx2 : None, str or list
        read2 of PE reads
        
    check:
    1. file exists
    2. file type
    3. fq paired
    4. check_empty
    """
    args = {
        'show_log': False
    }
    args.update(kwargs)
    if isinstance(fx1, str):
        c1 = check_fx(fx1)
        c2 = check_fx(fx2)
        c3 = check_fx_paired(fx1, fx2):
        if fx2 is None:
            out = c1
        else:
            out = all([c1, c2, c3])
        # show log
        if args['show_log']:
            msg = '\n'.join([
                '-'*30,
                '{:<5s} : {}'.format('ok' if c1 else 'failed', fx1),
                '{:<5s} : {}'.format('ok' if c2 else 'failed', fx2),
                '{:<5s} : is_paired'.format('yes' if c3 else 'no'),
                '-'*30,
                ])
            print(msg)
    elif isinstance(fx1, list):
        if isinstance(fx2, list):
            if len(fx1) == len(fx2):
                out = [check_fx_args(i,j) for i,j in zip(fx1, fx2)]
            else:
                out = False
        elif fx2 is None:
            out = all(check_fx(i) for i in fx1)
        else:
            out = False
    else:
        out = False
    return out


def read_lines(x, nrows=0, skip=0, strip_white=True, comment=''):
    """Read plain text file; warning for large file
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
    if check_file(x, check_empty=True):
        l = []
        try:
            i = 0
            with xopen(x) as r:
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
        log.error('read_lines() failed: {}'.format(x))
    return out


def str_distance(x, y, partial=True):
    """Check distance between a, b
    """
    out = -1
    if isinstance(x, str) and isinstance(y, str):
        try:
            if partial:
                x = x[:len(y)]
                y = y[:len(x)]
            out = lev.distance(x, y)
        except:
            out = -1 # huge number
    return out


#----------------------------#
def save2file(s, f, overwrite=False):
    """Save str to file"""
    if isinstance(s, str):
        pass
    elif isinstance(s, list):
        s = '\n'.join(list(map(str, s))) # to str
    else:
        log.error('s expect str, got {}'.format(type(s).__name__))
        return None
    # write to file
    if file_exists(f) and overwrite is False:
        log.info('save2file() skipped, file exists, {}'.format(f))
    else:
        if isinstance(f, str):
            if not check_path(os.path.dirname(f), create_dir=True):
                log.error('failed write to file: {}'.format(f))
                return None
            try:
                with open(f, 'wt') as w:
                    w.write(s+'\n')
            except IOError as e:
                log.error(e)
        else:
            log.error('f expect str, got {}'.format(type(f).__name__))


#----------------------------#
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

#     default: $HOME/data/genome/{genome}
#     """
#     def __init__(self, genome, genome_path=None,
#         repeat_masked_genome=False, **kwargs):
#         assert isinstance(genome, str)
#         self.genome = genome
#         self.repeat_masked_genome = repeat_masked_genome
#         self.kwargs = kwargs
#         if genome_path is None:
#             genome_path = os.path.join(str(Path.home()),
#                 'data', 'genome')
#         self.genome_path = genome_path


#     def get_fa(self):
#         """
#         Get the fasta file of specific genome
#         {genome}/bigZips/{genome}.fa
#         also check ".gz" file
#         """
#         fa = os.path.join(self.genome_path, self.genome, 'bigZips',
#             self.genome + '.fa')
#         if not file_exists(fa):
#             # gencode version
#             fa = os.path.join(self.genome_path, self.genome, 'fasta',
#                 self.genome + '.fa')
#         fa_gz = fa + '.gz'
#         if not file_exists(fa):
#             if file_exists(fa_gz):
#                 log.error('require to unzip the fasta file: %s' % fa_gz)
#             else:
#                 log.error('fasta file not detected: %s' % fa)
#             return None
#         else:
#             return fa


#     def get_fasize(self):
#         """Get the fasta size file, chromosome size
#         optional, fetch chrom size from ucsc
#         http://hgdownload.cse.ucsc.edu/goldenPath/<db>/bigZips/<db>.chrom.sizes

#         or using UCSC tool: fetchChromSizes
#         fetchChromSizes hg39 > hg38.chrom.sizes
#         """
#         fa = self.get_fa()
#         fa_size = fa + '.chrom.sizes'
#         if not file_exists(fa_size):
#             log.warning('file not exists, run samtools faidx to generate it')
#             pysam.faidx(fa) # create *.fa.fai
#             os.rename(fa + '.fai', fa_size)
#         return fa_size


#     def phylop100(self):
#         """Return the phylop100 bigWig file of hg19, only
#         for conservation analysis
#         """
#         p = os.path.join(self.genome_path, self.genome, 'phyloP100way',
#             self.genome + '.100way.phyloP100way.bw')
#         if not file_exists(p):
#             p = None
#         return p


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
#         if not file_exists(g):
#             g = None
#         return g


#     def gene_gtf(self, version='refseq'):
#         """Return the gene annotation in GTF format
#         support refseq, ensembl, gencode
#         """
#         version = version.lower() #
#         gtf = os.path.join(
#             self.genome_path,
#             self.genome,
#             'annotation_and_repeats',
#             self.genome + '.' + version + '.gtf')
#         if not file_exists(gtf):
#             gtf = os.path.join(
#             self.genome_path,
#             self.genome,
#             'gtf',
#             self.genome + '.' + version + '.gtf')
#         if not file_exists(gtf):
#             gtf = None
#         return gtf


#     def te(self, format='gtf'):
#         """Return TE annotation of the genome
#         or return TE consensus sequence for the genome (dm3)
#         """
#         # only dm3 supported
#         te_gtf = os.path.join(self.genome_path, self.genome,
#             self.genome + '_transposon',
#             self.genome + '_transposon.' + format)
#         if not file_exists(te_gtf):
#             te_gtf = None
#         return te_gtf


#     def piRNA_cluster(self, format='gtf'):
#         """Return TE annotation of the genome
#         or return TE consensus sequence for the genome (dm3)
#         """
#         # only dm3 supported
#         te_gtf = os.path.join(self.genome_path, self.genome,
#             self.genome + '_piRNA_clusters',
#             self.genome + '_piRNA_clusters.' + format)
#         if not file_exists(te_gtf):
#             te_gtf = None
#         return te_gtf
    
    
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

