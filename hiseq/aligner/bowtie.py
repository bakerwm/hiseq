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
from hiseq.utils.helper import *




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


    def prep_cmd(self):
        """
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
        arg_input = '{}'.format(self.fq1) if not self.is_paired \
            else '-1 {} -2 {}'.format(self.fq1, self.fq2)

        self.cmd = ' '.join([
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


    def parse_align(self, to_toml=True):
        """Wrapper bowtie log

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
        with open(self.align_log) as r:
            for line in r:
                # if not ':' in line or line.startswith('Warning'):
                #     continue
                if not line.startswith('#'):
                    continue
                num = line.strip().split(':')[1]
                value = num.strip().split(' ')[0]
                value = eval(value)
                if 'reads processed' in line:
                    dd['total'] = value
                elif 'at least one reported alignment' in line:
                    dd['map'] = value
                elif 'failed to align' in line:
                    dd['unmap'] = value
                elif 'alignments suppressed due to -m' in line:
                    dd['multiple'] = value
                else:
                    pass

        # unique_only
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


    def run(self):
        # run cmd
        if file_exists(self.bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(self.bam))
        else:
            try:
                run_shell_cmd(self.cmd)
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





