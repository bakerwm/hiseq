#!/usr/bin/env python3

"""Functions for aligner
check fastq/index, arguments, ...

file existence 
arguments 

Including functions/classes:

AlignReader
check_fx_args
get_args # for Align
"""
import os
import sys
import logging
import argparse
import hiseq
from hiseq.utils.utils import Config, update_obj, log
from hiseq.utils.file import check_file, check_fx, check_fx_paired, file_exists



logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


def default_values(x):
    """
    load default values, saved in package file, in toml format
    hiseq/bin/config.toml
    
    supported_genomes:
    supported_aligner:
    """
    pkg_dir = os.path.dirname(hiseq.__file__)
    f = os.path.join(pkg_dir, 'bin', 'config.toml')
    if file_exists(f):
        df = Config().load(f)
        out = df.get(x, None)
    else:
        out = None
    return out


class AlignReader(object):
    """Read the ALignment directory
    report the details
    config
    files
    map
    """
    def __init__(self, x, **kwargs):
        self.x = x
        self.init_args()
        self.read()
    
    
    def init_args(self):
        """Make sure, x is a dir"""
        self.isdir = False
        if isinstance(self.x, str):
            self.isdir = os.path.isdir(self.x)
    
    
    def list_config(self):
        """List the config.toml file
        in the path:
        x/config/config.toml
        
        old version:
        x/config.pickle
        """
        out = None
        if self.isdir:
            c = [os.path.join(self.x, 'config.toml'),
                 os.path.join(self.x, 'config.pickle'),
                 os.path.join(self.x, 'config', 'config.toml'),
                 os.path.join(self.x, 'config', 'config.pickle')]
            for f in c:
                if os.path.isfile(f):
                    out = f
                    break
        # message
        if out is None:
            log.error('x is not align directory: {}'.format(self.x))
        return out
    
    
    def read(self):
        config = self.list_config()
        p = Config().load(config) if config else {}
        # alignment type
        self.is_align = 'hiseq_type' in p
        self.is_hiseq = self.is_align
        self.hiseq_type = p.get('hiseq_type', None)
        self.is_align_r1 = self.hiseq_type == 'align_r1'
        self.is_align_rn = self.hiseq_type == 'align_rn'
        # for alignment_r1, parse files, stat
        self = update_obj(self, p, force=True)
        # bam, unmap1, unmap
        self.bam = getattr(self, 'bam', None)
        self.unmap1 = getattr(self, 'unmap1', None)
        self.unmap2 = getattr(self, 'unmap2', None)
        self.align_json = getattr(self, 'align_json', None)
        # map, unmap, ...
        q = Config().load(self.align_json) if self.align_json else {}
        self.total = q.get('total', 0)
        self.map = q.get('map', 0)
        self.unique = q.get('unique', 0)
        self.multi = q.get('multi', 0)
        self.unmap = q.get('unmap', 0)
        
    
    def get_output(self):
        """Return the bam, unmap1, unmap2 files"""
        return (self.bam, self.unmap1, self.unmap2)
    

def check_fx_args(fq1, fq2=None, **kwargs):
    """
    Check the fastx in arguments
    
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
    return out # (fq1, fq2) if out else None


## fq files: deprecated: by hiseq.utils.file.check_fx()
# def check_fx(fx, **kwargs):
#     """Check the fastq/a files
#     1. file exist
#     2. fq1 required
#     3. fq type: fasta/q
#     """
#     kwargs['check_empty'] = kwargs.get('check_empty', True)
#     kwargs['show_error'] = kwargs.get('show_error', False)
#     if isinstance(fx, str):
#         if check_file(fx, **kwargs):
#             try:
#                 fx_type = Fastx(fx).format # fasta/q
#                 out = fx_type in ['fasta', 'fastq']
#             except ValueError as err:
#                 log.info('Failed to read file, with error: {}'.format(err))
#                 out = False
# #             finally:
# #                 log.info('check_fx() ...')            
#         else:
#             if kwargs['show_error']:
#                 log.error('fx failed, {}'.format(fx))
#             out = False
#     else:
#         out = False
#     return out


## deprecated: by hiseq.utils.file.check_fx_paired()
# def check_fx_paired(fq1, fq2, **kwargs):
#     """Check the fq1 and fq2 are paired or not
#     the name of fq1, fq2 are the same

#     fq1: @the-name/1
#     fq2: @the-name/2

#     scan the first read from the two files
#     """
#     if isinstance(fq1, str) and isinstance(fq2, str):
#         if check_fx(fq1, **kwargs) and check_fx(fq2, **kwargs):
#             # check fq name:
#             fx1 = pyfastx.Fastx(fq1)
#             fx2 = pyfastx.Fastx(fq2)
#             for a,b in zip(fx1, fx2):
#                 out = a[0][:-1] == b[0][:-1]
#                 break
#         else:
#             out = False
#     else:
#         out = False
#     return out

## deprecated: by hiseq.utils.file.check_file()
# def check_file(x, **kwargs):
#     """Check the x file
#     1. file exists
    
#     Parameters
#     ----------
#     x : str
#         Path to a file
    
#     Keyword Parameters
#     ------------------
#     show_error : bool
#         Show the error messages
        
#     show_log : bool
#         Show the log messages 
        
#     check_empty : bool
#         Check if the file is empty or not,  gzipped empty file, size=20
        
#     emptycheck : bool
#         see check_empty
#     """
#     show_error = kwargs.get('show_error', False)
#     show_log = kwargs.get('show_log', False)
#     check_empty = kwargs.get('check_empty', False)
#     emptycheck = kwargs.get('emptycheck', False) # for old version
#     if isinstance(x, str):
#         if file_exists(x):
#             x_size = os.stat(x).st_size
#             # empty gzipped file, size=20
#             q_size = 20 if x.endswith('.gz') else 0
#             out = x_size > q_size if check_empty or emptycheck else True
#             if show_log:
#                 flag = 'ok' if out else 'failed'
#                 log.info('{:<6s} : {}'.format(flag, x))
#         else:
#             if show_error:
#                 log.error('file not exists: {}'.format(x))
#             out = False # failed
#     elif isinstance(x, list):
#         out = all([check_file(i, **kwargs) for i in x])
#     else:
#         log.error('x expect str or list, got {}'.format(type(x).__name__))
#         out = False
#     return out



def get_args():
    """Parsing arguments for Align
    The main port, aligner
    """
    example = '\n'.join([
        'Examples:',
        '$ python align.py -1 f1.fq -x genome -o output',
        '# add extra para',
        '$ python align.py -1 f1.fq -2 f2.fq -x genome -o output -X "-X 2000"',
        '# unique reads, update index_name',
        '$ python align.py -1 f1.fq -x genome -o output -u -in 01.genome',
    ])    
    parser = argparse.ArgumentParser(
        prog='align',
        description='run algner {bowtie|bowtie2|STAR}',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-a', '--aligner', default='bowtie2', type=str,
                        help='The aligner for alignment, default: [bowtie2]')
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
                        help='Fasta/q file, read1 of PE, or SE read')
    parser.add_argument('-2', '--fq2', nargs='+', required=False, default=None,
                        help='Fasta/q file, read2 of PE, or SE read, optional')
    parser.add_argument('-o', '--outdir', default=None,
                        help='Directory saving results, default: [cwd]')
    parser.add_argument('-g', '--genome', default=None, 
                        help='The name of the genome, [dm6, hg38, mm10]')
    parser.add_argument('-x', '--genome-index', default=None,
                        help='The path to the alignment index')
    parser.add_argument('-n', '--smp-name', nargs='+', default=None, 
                        dest='smp_name',
                        help='The name of the sample')
    parser.add_argument('--spikein', default=None, type=str,
                        help='The genome name of spikein, default: None')
    parser.add_argument('--spikein-index', default=None, dest='spikein_index',
                        help='The alignment index of spikein, default: [None]')
    parser.add_argument('-p', '--threads', default=1, type=int,
                        help='Number of threads, default: [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
                        help='Number of jobs to run in parallel, default: [1]')
    parser.add_argument('-w', '--overwrite', action='store_true',
                        help='Overwrite the exist files')
    parser.add_argument('-u', '--unique-only', action='store_true',
                        dest='unique_only', 
                        help='Report unique mapped reads only')
    parser.add_argument('--n-map', type=int, default=1,
                        help='Number of hits per read')
    parser.add_argument('-l', '--largs-insert', action='store_true',
                        dest='large_insert',
                        help='For large insert, use: -X 1000 --chunkmbs 128')
    parser.add_argument('--clean', dest='keep_tmp', action='store_false',
                        help='Clean temp files')
    parser.add_argument('-X', '--extra-para', dest='extra_para', default=None,
                        help='Add extra parameters, eg: "-X 2000"')
    parser.add_argument('--to-rRNA', action='store_true', dest='to_rRNA', 
                        help='Align reads to rRNA first')
    parser.add_argument('--to-chrM', action='store_true', dest='to_chrM',
                        help='Align reads to mitochromosome first')
    parser.add_argument('--to-MT-trRNA', action='store_true', dest='to_MT_trRNA',
                        help='Align reads to chrM, tRNA and rRNAs first')
    parser.add_argument('--verbose', action='store_true', 
                        help='Show message in details')
    return parser
