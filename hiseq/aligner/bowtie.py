#!/usr/bin/env python3

"""
Align reads to index, using bowtie 

Basic usage: 

1. SE
bowtie -S -x index in.fq > out.sam 2> out.log

2. PE
bowtie2 -S -x index -1 r1.fq -2 r2.fq > out.sam 2> out.log


Updates:

1. style the code with PEP8: https://www.python.org/dev/peps/pep-0008/

in brief:

1 indentation: 4 spaces  
    1.1 Aligned with opening delimiter
2 Maximum Line Length: 79 characters 
    2.1 Multi lines, operator at end
3 Blank Lines
    3.1 Two blank lines surround top-level function/class  
    3.2 Single blank line, in class  
    3.3 
4 White spaces  
    4.1 no spaces after/before brackets/braces  
    4.2 single space on both sides of binary operator   
    4.3 no spaces around '=' for keyword argument  

5 Naming styles

5.1 lower_case_with_underscores, 

6 Package and Module Names 

6.1 Modules should have short, all-lowercase names

7 Class names should normally use the CapWords convention. 

8 Function and Variable Names,  lowercase, with words separated by underscores as necessary to improve readability.

9 Use `is`, `is not`

10 Use ''.startswith(), ''.endswith()

"""


import os
import sys
from Levenshtein import distance
from utils import *
from aligner_index import *
from hiseq.utils.helper import *




def bowtie_parser(x):
    """Wrapper bowtie directory
    Bowtie:
    # reads processed: 10000
    # reads with at least one reported alignment: 3332 (33.32%)
    # reads that failed to align: 457 (4.57%)
    # reads with alignments suppressed due to -m: 6211 (62.11%)

    or:

    # reads processed: 10000
    # reads with at least one reported alignment: 9543 (95.43%)
    # reads that failed to align: 457 (4.57%)

    unique, multiple, unmap, map, total

    skip: Warning, ...
    """
    dd = {}
    total = 0
    mapped = 0
    unmapped = 0
    multi = 0
    unique = 0
    if check_file(x, check_empty=True):
        with open(x) as r:
            for line in r:
                if not line.startswith('#'):
                    continue
                num = line.strip().split(':')[1]
                value = num.strip().split(' ')[0]
                value = eval(value)
                if 'reads processed' in line:
                    total = value
                elif 'at least one reported alignment' in line:
                    mapped = value
                elif 'failed to align' in line:
                    unmapped = value
                elif 'alignments suppressed due to -m' in line:
                    multi = value
                else:
                    pass
        # unique_only
        if unique_only:
            mapped = unqiue
        else:
            mapped = unique + multi
            
        dd['unique'] = dd['map']
        dd['multiple'] = dd.get('multiple', 0) # default 0

        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']

        # save fqname, indexname,
        dd['fqname'] = self.smp_name
        dd['index_name'] = self.index_name

        # # sort by keys
        self.log_dict = dd

    # save dict to plaintext file
    with open(self.align_stat, 'wt') as w:
        ## version-1
        # for k, v in sorted(dd.items()):
        #     w.write('\t'.join([self.config.fqname, self.config.index_name, k, str(v)]) + '\n')

        # ## version-2
        # w.write('#') # header line
        # w.write('\t'.join(list(map(str, dd.keys()))) + '\n')
        # w.write('\t'.join(list(map(str, dd.values()))) + '\n')

        ## version-3
        groups = ['total', 'map', 'unique', 'multiple', 'unmap', 
            'fqname', 'index_name']
        h = '\t'.join(groups)
        v = '\t'.join([str(dd.get(i, 0)) for i in groups])
        w.write('#' + h + '\n')
        w.write(v + '\n')

    if to_toml:
        Toml(dd).to_toml(self.align_toml)

    return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']




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
        }
        self = update_obj(self, args_init, force=False)
        self.init_fx()
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
            'unmap': prefix + '.unmap.fastq',
            'unmap1': prefix + '.unmap.1.fastq',
            'unmap2': prefix + '.unmap.2.fastq',
            'align_log': prefix + '.align.log',
            'align_stat': prefix + '.align.stat',
            'align_toml': prefix + '.align.toml',
            'align_flagstat': prefix + '.flagstat'
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
        self.prep_cmd()


    def init_args(self):
        """
        check
        """
        args_local = AlignerConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'bowtie' # force changed.
        Toml(self.__dict__).to_toml(self.config_toml)        


    def get_cmd(self):
        """The command line
        unique
        n_map
        extra_para
        """
        # unique/multiple
        if self.extra_para is None:
            self.extra_para = ''
        if self.n_map < 1:
            self.n_map = 1
        if self.unique_only:
            arg_unique = '-m 1'
        else:
            arg_unique = '-v 2 -k {}'.format(self.n_map)
        # file type
        arg_fx = '-f' if self.fq_format == 'fasta' else '-q'
        # input file
        arg_input = '-1 {} -2 {}'.format(self.fq1, self.fq2) if \
            check_fx(self.fq1, self.fq2) else self.fq1
        # command-line
        cmd = ' '.join([
            '{}'.format(shutil.which('bowtie')),
            '--mm --best --sam --no-unal',
            arg_unique, 
            arg_fx,
            self.extra_para,
            '--un {}'.format(self.unmap),
            '{}'.format(self.index),
            arg_input,
            '2> {}'.format(self.align_log),
            '&& samtools view -@ {}'.format(self.threads),
            '-Sub -F 0x4 -',
            '| samtools sort -@ {}'.format(self.threads),
            '-o {} -'.format(self.bam),
            '&& samtools index {}'.format(self.bam),
            '&& samtools flagstat {} > {}'.format(
                self.bam, 
                self.align_flagstat),
        ])
        # save cmd
        with open(self.cmd_shell, 'wt') as w:
            w.write(self.cmd + '\n')
        return cmd

    
    def run(self):
        if file_exists(self.bam) and not self.overwrite:
            log.info('bowtie() skipped, file exists: {}'.format(self.bam))
        else:
            cmd = self.get_cmd()
            try:
                run_shell_cmd(cmd)
            except:
                log.error('align() failed, check {}'.format(self.align_log))
        # rename unmap files
        unmap1 = self.subdir + '/' + self.smp_name + '.unmap_1.fastq'
        unmap2 = self.subdir + '/' + self.smp_name + '.unmap_2.fastq'
        if self.is_paired:
            file_copy(unmap1, self.unmap1)
            file_copy(unmap2, self.unmap2)
            file_remove([unmap1, unmap2], ask=False)
        else:
            self.unmap1 = self.unmap
            self.unmap2 = None
        # parse log file
        if file_exists(self.align_log):
            self.parse_align(to_toml=True)
        # temp files
        del_list = [self.sam, self.unmap1, self.unmap2]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)
        return (self.bam, self.unmap1, self.unmap2)





