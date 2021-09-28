#!/usr/bin/env python3

"""
Align reads to index, using salmon

Basic usage: 

salmon quant -i $index -l A -p 8 --gcBias -o ${outdir}/$prefix -1 $r1 -2 $r2 
# output
# quant.sf

tximport()
txmeta()
"""



import os
import sys
import re
import shutil
import pathlib
import argparse
from Levenshtein import distance
from hiseq.utils.seq import Fastx
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd
from hiseq.utils.file import check_fx, check_file, check_path, file_abspath, \
    file_exists, fx_name, check_fx_args
from hiseq.align.align_index import AlignIndex
# from hiseq.align.utils import check_fx_args


def parse_salmon(x):
    """Wrapper salmon log
    
    # PE
    Only 895980 fragments were mapped, but the number of burn-in fragments was set to 5000000.
    Observed 1000000 total fragments (1000000 in most recent round)
    # SE
    Only 895979 fragments were mapped, but the number of burn-in fragments was set to 5000000.
    Observed 1000000 total fragments (1000000 in most recent round)
    
    ## salmon version 1.4.0
    Observed 10259420 total fragments
    Counted 8,433,361 total reads in the equivalence classes
    
    unique, multiple, unmap, map, total
    """
    total = 0
    mapped = 0
    unmapped = 0
    multi = 0
    unique = 0
    out = 0
    if check_file(x, check_empty=True):
        with open(x) as r:
            for line in r:
                p1 = re.search('Counted ([0-9,]+) total reads in the', line)
                p2 = re.search('Observed (\d+) total fragments', line)
                if p1:
                    mapped = eval(p1.group(1).replace(',', ''))
                elif p2:
                    total = eval(p2.group(1))
                else:
                    pass
    else:
        log.error('file not exists, {}'.format(x))
    return {
        'total': total,
        'map': mapped,
        'unique': mapped,
        'multi': multi,
        'unmap': total - mapped,
        }


class SalmonConfig(object):
    """Check args, prepare files for salmon
    arguments
    output files
    parser ?!
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'salmon',
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'index': None,            
            'index_name': None,
            'smp_name': None,
            'extra_para': None,
            'smp_name': None,
            'threads': 1,
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'salmon_r1'
        self.init_fx()
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        # index name
        if not AlignIndex(self.index, self.aligner).is_valid():
            raise ValueError('index not valid, {}'.format(self.index))
        if self.index_name is None:
            self.index_name = AlignIndex(self.index).index_name()
        # sample
        if self.smp_name is None:
            self.smp_name = fx_name(self.fq1, fix_pe=True)
        # update files
        self.init_files()


    def init_fx(self):
        """Make sure, fx
        1. str
        2. fq1 exists
        3. fq2 None or exists
        """
        if not check_fx(self.fq1):
            raise ValueError('--fq1, not exists, or empty')
        if self.fq2 is not None and not check_fx(self.fq2):
            raise ValueError('--fq2, not exists, or empty')
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('--fq1, --fq2 faild, not properly paired')
        # format
        self.fx_format = Fastx(self.fq1).format # fasta/q
        self.is_paired = file_exists(self.fq2) # fq2
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)


    def init_files(self):
        self.project_dir = os.path.join(self.outdir, self.smp_name, 
            self.index_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        # output files
        prefix = os.path.join(self.project_dir, self.smp_name)
        default_files = {
            'project_dir': self.project_dir,
            'config_yaml': os.path.join(self.config_dir, 'config.yaml'),
            'cmd_shell': os.path.join(self.project_dir, 'cmd.sh'),
            'quant_sf': os.path.join(self.project_dir, 'quant.sf'),
            'bam': prefix + '.bam',
            'sam': prefix + '.sam',
            'unmap': prefix + '.unmap.' + self.fx_format,
            'unmap1': prefix + '.unmap.1.' + self.fx_format, # 
            'unmap2': prefix + '.unmap.2.' + self.fx_format, #
            'align_log': os.path.join(self.project_dir, 'logs', 'salmon_quant.log'),
            'run_log': os.path.join(self.project_dir, 'run.log'),
            'align_stat': prefix + '.align.stat',
            'align_json': prefix + '.align.json',
            'align_flagstat': prefix + '.flagstat',
        }
        self = update_obj(self, default_files, force=True)
        check_path([self.project_dir, self.config_dir], create_dirs=True)


class Salmon(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = SalmonConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'salmon' # force changed
        self.get_cmd()        
        Config().dump(self.__dict__, self.config_yaml)


    def get_cmd(self):
        """
        Example:
        salmon quant -i $index --gcBias -l A -p 8 -o ${outdir}/$prefix -1 $r1 -2 $r2
        """
        # fq 
        cmd_fx = '-1 {} -2 {}'.format(self.fq1, self.fq2) if self.fq2 else '-r {}'.format(self.fq1)
        self.cmd = ' '.join([
            '{} quant'.format(shutil.which('salmon')),
            '--gcBias -l A',
            '-p {}'.format(self.threads),
            '-i {}'.format(self.index),
            '-o {}'.format(self.project_dir),
            cmd_fx
        ])
    

    def run(self):
        if file_exists(self.quant_sf) and not self.overwrite:
            log.info('Salmon() skipped, file exists: {}'.format(self.quant_sf))
        else:
            with open(self.cmd_shell, 'wt') as w:
                w.write(self.cmd + '\n')
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('Salmon() failed, check {}'.format(self.align_log))
            # output file
            if not check_file(self.quant_sf, check_empty=True):
                log.error('Salmon() failed, check {}'.format(self.align_log))
        df = parse_salmon(self.align_log)
        df.update({
            'name': self.smp_name,
            'index': self.index_name
            })
        Config().dump(df, self.align_json)
        return (self.quant_sf, None, None) # bam, unmap1, unmap2


def get_args():
    """Parsing arguments for salmon
    """
    example = '\n'.join([
        'Examples:',
        '$ python salmon.py -1 f1.fq -2 f2.fq -x index -o output -p 8',
    ])    
    parser = argparse.ArgumentParser(
        prog='run_salmon',
        description='run salmon program',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
                        help='Fasta/q file, read1 of PE, or SE read')
    parser.add_argument('-2', '--fq2', required=False, default=None,
                        help='Fasta/q file, read2 of PE, or SE read, optional')
    parser.add_argument('-x', '--index', required=True,
                        help='The alignment index for salmon')
    parser.add_argument('-o', '--outdir', default=None,
                        help='Directory saving results, default: [cwd]')
    parser.add_argument('-in', '--index-name', default=None, dest='index_name',
                        help='The name of the index')
    parser.add_argument('-n', '--smp-name', default=None, dest='smp_name',
                        help='The name of the sample')
    parser.add_argument('-p', '--threads', default=1, 
                        help='Number of threads, default: [1]')
    parser.add_argument('-w', '--overwrite', action='store_true',
                        help='Overwrite the exist files')
    return parser


def main():
    args = vars(get_args().parse_args())
    Salmon(**args).run()


if __name__ == '__main__':
    main()

