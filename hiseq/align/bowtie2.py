#!/usr/bin/env python3

"""
Align reads to index, using bowtie2

Basic usage: 

1. SE
bowtie2 -S -x index in.fq > out.sam 2> out.log

2. PE
bowtie2 -S -x index -1 r1.fq -2 r2.fq > out.sam 2> out.log
"""

import os
import sys
import re
import shutil
import argparse
from Levenshtein import distance
from hiseq.utils.seq import Fastx
from hiseq.utils.utils import update_obj, Config, log, run_shell_cmd
from hiseq.utils.file import check_fx, check_file, check_path, file_abspath, file_exists, fx_name, symlink_file, copy_file, remove_file
from hiseq.align.align_index import AlignIndex
from hiseq.align.utils import check_fx_args, AlignReader



def parse_bowtie2(x):
    """Wrapper bowtie2 log
    for PE: only proper paired reads

    SE:

    10000 reads; of these:
      10000 (100.00%) were unpaired; of these:
        166 (1.66%) aligned 0 times
        2815 (28.15%) aligned exactly 1 time
        7019 (70.19%) aligned >1 times
    98.34% overall alignment rate

    PE:

    100000 reads; of these:
      100000 (100.00%) were paired; of these:
        92926 (92.93%) aligned concordantly 0 times
        5893 (5.89%) aligned concordantly exactly 1 time
        1181 (1.18%) aligned concordantly >1 times
        ----
        92926 pairs aligned concordantly 0 times; of these:
          1087 (1.17%) aligned discordantly 1 time
        ----
        91839 pairs aligned 0 times concordantly or discordantly; of these:
          183678 mates make up the pairs; of these:
            183215 (99.75%) aligned 0 times
            101 (0.05%) aligned exactly 1 time
            362 (0.20%) aligned >1 times
    8.39% overall alignment rate

    unique, multiple, unmap, map, total

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
                n = line.strip().split()[0]
                if not re.match('^[0-9]+$', n):
                    continue
                n = eval(n)
                # parsing
                if 'were paired; of these:' in line:
                    total = n
                    is_paired = True
                elif 'aligned concordantly 0 times' in line:
                    unmapped = n
                elif 'aligned concordantly exactly 1 time' in line:
                    unique = n
                elif 'aligned concordantly >1 times' in line:
                    multi = n
                elif 'reads; of these' in line and not is_paired:
                    total = n
                elif 'aligned exactly 1 time' in line and not is_paired:
                    unique = n
                elif 'aligned >1 times' in line and not is_paired:
                    multi = n
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
        'map': unique + multi,
        'unique': unique,
        'multi': multi,
        'unmap': unmapped,
        }


class Bowtie2Config(object):
    """Check args, prepare files for bowtie2
    arguments
    output files
    parser ?!
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie2',
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
            'n_map': 0,
            'unique_only': False,
            'keep_tmp': False,
            'keep_unmap': True,
            'large_insert': False,
            'default_bowtie2': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'bowtie2_r1'
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


    def init_files(self):
        self.project_dir = os.path.join(self.outdir, self.smp_name, 
            self.index_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        # output files
        prefix = os.path.join(self.project_dir, self.smp_name)
        default_files = {
#             'project_dir': self.project_dir,
            'config_toml': os.path.join(self.config_dir, 'config.toml'),
            'cmd_shell': os.path.join(self.project_dir, 'cmd.sh'),
            'bam': prefix + '.bam',
            'sam': prefix + '.sam',
            'unmap': prefix + '.unmap.' + self.fx_format,
            'unmap1': prefix + '.unmap.1.' + self.fx_format, # 
            'unmap2': prefix + '.unmap.2.' + self.fx_format, #
            'align_log': prefix + '.align.log',
            'align_stat': prefix + '.align.stat',
            'align_json': prefix + '.align.json',
            'align_flagstat': prefix + '.flagstat',
        }
        self = update_obj(self, default_files, force=True)
        check_path([self.project_dir, self.config_dir], create_dirs=True)
        

class Bowtie2(object):
    """Alignment, using Bowtie2
    Single index, SE/PE

    1. SE
    bowtie2 -S -x index in.fq > out.sam 2> out.log

    2. PE
    bowtie2 -S -x index -1 r1.fq -2 r2.fq > out.sam 2> out.log
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = Bowtie2Config(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'bowtie2' # force changed
        self.get_cmd()        
        Config().dump(self.__dict__, self.config_toml)


    def get_cmd(self):
        """The command line
        ########################
        # Bowtie2 unique reads #
        ########################
        filt by tag: YT:Z:CP        
        YT:Z: String representing alignment type
        CP: Concordant; DP: Discordant; UP: Unpaired Mate; UU: Unpaired.
        answer-1: https://www.biostars.org/p/19283/#19292
        answer-2: https://www.biostars.org/p/133094/#133127
         
        -F 2048: suppress Supplementary alignments
        
        or set: -q 10 

        unique: -m 1
        n_map: -k x
        extra_para: 'extra'
        # example:
        """
        args_fmt = '-f' if self.fx_format == 'fasta' else '-q'
        args_nmap = '-k {}'.format(self.n_map) if self.n_map > 0 else ''
        args_common = '--local --very-sensitive --no-unal --no-mixed --no-discordant'
        args_extra = self.extra_para if self.extra_para else ''
        args_large_ins = '-X 2000' if self.large_insert else ''
        if self.default_bowtie2:
            args_nmap = ''
            args_common = ''
            args_extra = ''
            args_large_ins = ''
        if self.is_paired:
            args_io = ' '.join([
                '--un-conc {}'.format(self.unmap),
                '-1 {} -2 {}'.format(self.fq1, self.fq2),
                ])
        else:
            args_io = '--un {} {}'.format(self.unmap, self.fq1)
        # command-line
        cmd_main = ' '.join([
            '{}'.format(shutil.which('bowtie2')),
            '--mm -p {}'.format(self.threads),
            args_fmt,
            args_common,
            args_nmap,
            args_extra,
            args_large_ins,
            '-x {}'.format(self.index),
            args_io,
            '1> {} 2> {}'.format(self.sam, self.align_log),
        ])
        # for unique
        if self.is_paired:
            if self.unique_only:
                cmd_unique = ' '.join([
                    '&& samtools view -Sub',
                    '<(samtools view -H {};'.format(self.sam),
                    'samtools view -F 2048 {}'.format(self.sam),
                    "| grep 'YT:Z:CP')",
                    ])
            else:
                cmd_unique = ' '.join([
                    '&& samtools view -Sub -F 2048 {}'.format(self.sam),
                    ])
        else:
            if self.unique_only:
                cmd_unique = ' '.join([
                    '&& samtools view -Sub -q 10 -F 2048 {}'.format(self.sam)
                    ])
            else:
                cmd_unique = ' '.join([
                    '&& samtools view -Sub -F 2048 {}'.format(self.sam)
                    ])
        # add cmd
        self.cmd = ' '.join([
            cmd_main, 
            cmd_unique,
            '| samtools sort -@ {} -o {} -'.format(self.threads, self.bam),
            '&& samtools index {}'.format(self.bam),
            '&& samtools flagstat {} > {}'.format(self.bam, self.align_flagstat)
            ])


    def run(self):
        if file_exists(self.bam) and not self.overwrite:
            log.info('Bowtie2() skipped, file exists: {}'.format(self.bam))
        else:
            # save cmd
            with open(self.cmd_shell, 'wt') as w:
                w.write(self.cmd + '\n')
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('Bowtie2() failed, check {}'.format(self.align_log))
            # output file
            if not check_file(self.bam, check_empty=True):
                log.error('Bowtie2() failed, check {}'.format(self.align_log))
        # rename unmap files, unmap_1.fastq -> unmap.1.fastq
        if not self.is_paired:
            self.unmap1, self.unmap2 = (self.unmap, None)
        # log
        df = parse_bowtie2(self.align_log)
        df.update({
            'name': self.smp_name,
            'index': self.index_name,
            'unique_only': self.unique_only,
            })
        Config().dump(df, self.align_json)
        del_list = [self.sam]
        if not self.keep_unmap:
            del_list.extend([self.unmap1, self.unmap2, self.unmap])
        if not self.keep_unmap:
            remove_file(del_list, ask=False)
        return (self.bam, self.unmap1, self.unmap2)


def get_args():
    """Parsing arguments for bowtie2
    """
    example = '\n'.join([
        'Examples:',
        '$ python bowtie2.py -1 f1.fq -x genome -o output',
        '# add extra para',
        '$ python bowtie2.py -1 f1.fq -2 f2.fq -x genome -o output -X "-X 2000"',
        '# unique reads, update index_name',
        '$ python bowtie2.py -1 f1.fq -x genome -o output -u -in 01.genome',
    ])    
    parser = argparse.ArgumentParser(
        prog='run_bowtie',
        description='run bowtie2 program',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
                        help='Fasta/q file, read1 of PE, or SE read')
    parser.add_argument('-2', '--fq2', required=False, default=None,
                        help='Fasta/q file, read2 of PE, or SE read, optional')
    parser.add_argument('-x', '--index', required=True,
                        help='The alignment index for bowtie2')
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
    Bowtie2(**args).run()


if __name__ == '__main__':
    main()