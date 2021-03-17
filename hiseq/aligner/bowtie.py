#!/usr/bin/env python3

"""
Align reads to index, using bowtie 

Basic usage: 

1. SE
bowtie -S -x index in.fq > out.sam 2> out.log

2. PE
bowtie -S -x index -1 r1.fq -2 r2.fq > out.sam 2> out.log
"""

import os
import sys
import argparse
from Levenshtein import distance
from hiseq.utils.helper import *
from aligner_index import *
from utils import * # overwrite: check_file


def parse_bowtie(x):
    """Wrapper bowtie directory
    Bowtie:
    # reads processed: 10000
    # reads with at least one alignment: 3332 (33.32%)
    # reads that failed to align: 457 (4.57%)
    # reads with alignments suppressed due to -m: 6211 (62.11%)

    or:

    # reads processed: 10000
    # reads with at least one alignment: 9543 (95.43%)
    # reads that failed to align: 457 (4.57%)

    total map unique multi unmap
    
    Warnings-1
    Warning: Exhausted best-first chunk memory for read ..., skipping read

    solution: --chunkmbs 200 --maxins 1000

    if -k 1, could not tell, unique and multi
    """
    total = 0
    mapped = 0
    unmapped = 0
    multi = 0
    unique = 0
    out = 0
    is_paired = False
    warn_chunkmbs = False
    if check_file(x, check_empty=True):
        # processing
        with open(x) as r:
            for line in r:
                if 'Exhausted best-first chunk memory' in line:
                    warn_chunkmbs = line
                if not line.startswith('#'):
                    continue
                # warnings
                if 'paired-end alignments' in line:
                    is_paired = True
                num = line.strip().split(':')[1]
                value = eval(num.strip().split(' ')[0])
                if 'reads processed' in line:
                    total = value
                elif 'reads with at least one' in line:
                    mapped = value
                elif 'reads that failed to' in line:
                    unmapped = value
                elif 'alignments suppressed due to -m' in line:
                    multi = value
                elif line.startswith('Reported'):
                    out = eval(line.split()[1])
                else:
                    pass
    else:
        log.error('file not exists, {}'.format(x))
    # msg
    if warn_chunkmbs:
        log.warning('{}\nset --chunkmbs 128, to fix the errors'.format(
            warn_chunkmbs))
    return {
        'total': total,
        'map': mapped,
        'unique': mapped - multi,
        'multi': multi,
        'unmap': unmapped,
        }
        


class BowtieConfig(object):
    """Check args, prepare files for bowtie
    arguments
    output files
    parser ?!
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie',
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
            'n_map': 1,
            'unique_only': False,
            'keep_tmp': False,
            'keep_unmap': True,
            'large_insert': False,
        }
        self = update_obj(self, args_init, force=False)
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
            self.smp_name = fq_name(self.fq1, pe_fix=True)
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
        
        
    def init_files(self):
        self.project_dir = os.path.join(self.outdir, self.smp_name, 
            self.index_name)
        # output files
        prefix = os.path.join(self.project_dir, self.smp_name)
        default_files = {
            'project_dir': self.project_dir,
            'config_toml': os.path.join(self.project_dir, 'config.toml'),
            'cmd_shell': os.path.join(self.project_dir, 'cmd.sh'),
            'bam': prefix + '.bam',
            'sam': prefix + '.sam',
            'unmap': prefix + '.unmap.' + self.fx_format,
            'unmap1': prefix + '.unmap.1.' + self.fx_format, # 
            'unmap2': prefix + '.unmap.2.' + self.fx_format, # 
            'unmap1_tmp': prefix + '.unmap_1.' + self.fx_format, # 
            'unmap2_tmp': prefix + '.unmap_2.' + self.fx_format, #
            'align_log': prefix + '.align.log',
            'align_stat': prefix + '.align.stat',
            'align_toml': prefix + '.align.toml',
            'align_flagstat': prefix + '.flagstat',
        }
        self = update_obj(self, default_files, force=True)
        check_path(self.project_dir)


class Bowtie(object):
    """Alignment, using Bowtie
    Single index, SE/PE

    1. SE
    bowtie -S -x index in.fq > out.sam 2> out.log

    2. PE
    bowtie2 -S -x index -1 r1.fq -2 r2.fq > out.sam 2> out.log
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = BowtieConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'bowtie' # force changed
        self.get_cmd()        
        Config().dump(self.__dict__, self.config_toml)


    def get_cmd(self):
        """The command line
        unique: -m 1
        n_map: -k x
        extra_para: 'extra'
        # example:
        bowtie --best -S --no-unal -x index fq1 | \
            samtools view -F 0x4 -Sub - | \
            samtools sort -o out.bam - && \
            samtools index out.bam && \
            samtools flagstat out.bam > out.stat

        # catch errors
        Warning: Exhausted best-first chunk memory   

        1. set "-X 1000" or more, for large insert, or use bowtie2 
        2. set '--chunkmbs 200' or more, default: 64
        """
        if self.n_map < 1:
            self.n_map = 1
        args_extra = self.extra_para if self.extra_para else ''
        args_large_ins = '--chunkmbs 200 -X 1000' if self.large_insert else ''
        args_unique = '-m 1' if self.unique_only else '-v 2 -k {}'.format(
            self.n_map)
        args_fmt = '-f' if self.fx_format == 'fasta' else '-q'
        if self.is_paired:
            args_io = ' '.join([
                '--un {}'.format(self.unmap),
                '-1 {} -2 {}'.format(self.fq1, self.fq2),
                ])
        else:
            args_io = '--un {} {}'.format(self.unmap, self.fq1)
        # command-line
        self.cmd = ' '.join([
            '{}'.format(shutil.which('bowtie')),
            '-p {}'.format(self.threads),
            '--mm --best --sam --no-unal',
            args_unique, 
            args_fmt,
            args_extra,
            args_large_ins,
            '-x {}'.format(self.index),
            args_io,
            '2> {}'.format(self.align_log),
            '| samtools view -Sub -F 0x4 -@ {} -'.format(self.threads),
            '| samtools sort -@ {} -o {} -'.format(self.threads, self.bam),
            '&& samtools index {}'.format(self.bam),
            '&& samtools flagstat {} > {}'.format(
                self.bam, self.align_flagstat),
        ])

    
    def run(self):
        if file_exists(self.bam) and not self.overwrite:
            log.info('bowtie() skipped, file exists: {}'.format(self.bam))
        else:
            # save cmd
            with open(self.cmd_shell, 'wt') as w:
                w.write(self.cmd + '\n')
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('Bowtie() failed, check {}'.format(self.align_log))
            # output file
            if not check_file(self.bam, check_empty=True):
                log.error('Bowtie() failed, check {}'.format(self.align_log))
        # rename unmap files, unmap_1.fastq -> unmap.1.fastq
        if self.is_paired:
            file_symlink(self.unmap1_tmp, self.unmap1)
            file_symlink(self.unmap2_tmp, self.unmap2)
        else:
            self.unmap1, self.unmap2 = (self.unmap, None)
        # log
        df = parse_bowtie(self.align_log)
        df.update({
            'name': self.smp_name,
            'index': self.index_name,
            'unique_only': self.unique_only,
            })
        Config().dump(df, self.align_toml)
        del_list = [
            self.unmap1, self.unmap2, 
            self.unmap1_tmp, self.unmap2_tmp, 
            self.unmap
            ]
        if not self.keep_unmap:
            file_remove(del_list, ask=False)
        return (self.bam, self.unmap1, self.unmap2)


def get_args():
    """Parsing arguments for plotHeatmaps, deeptools
    """
    example = '\n'.join([
        'Examples:',
        '$ python bowtie.py -1 f1.fq -x genome -o output',
        '# add extra para',
        '$ python bowtie.py -1 f1.fq -2 f2.fq -x genome -o output -X "-X 2000"',
        '# unique reads, update index_name',
        '$ python bowtie.py -1 f1.fq -x genome -o output -u -in 01.genome',
    ])    
    parser = argparse.ArgumentParser(
        prog='run_bowtie',
        description='run bowtie program',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
                        help='Fasta/q file, read1 of PE, or SE read')
    parser.add_argument('-2', '--fq2', required=False, default=None,
                        help='Fasta/q file, read2 of PE, or SE read, optional')
    parser.add_argument('-x', '--index', required=True,
                        help='The alignment index for bowtie')
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
    parser.add_argument('-u', '--unique-only', action='store_true',
                        dest='unique_only', 
                        help='Report unique mapped reads only')
    parser.add_argument('-l', '--largs-insert', action='store_true',
                        dest='large_insert',
                        help='For large insert, use: -X 1000 --chunkmbs 128')
    parser.add_argument('--clean', dest='keep_tmp', action='store_false',
                        help='Clean temp files')
    parser.add_argument('-X', '--extra-para', dest='extra_para', default=None,
                        help='Add extra parameters, eg: "-X 2000"')
    return parser


def main():
    args = vars(get_args().parse_args())
    # update: keep_tmp, keep_unmap
    args['keep_unmap'] = args['keep_tmp']
    Bowtie(**args)


if __name__ == '__main__':
    main()