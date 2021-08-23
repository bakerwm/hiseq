#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert bam to bigWig

Solution-1: BAM -> bedGraph -> bigWig (bedtools, bedGraphToBigWig)
Solution-2: BAM -> bigWig (bamCoverage)

How to normalize

1. total reads to 1M (CPM) 
2. RPGC (reads per genomic content, 1x)
3. RPKM per bin: 
4. CPM  per bin:
"""

import os
import pathlib
import argparse
from shutil import which
from hiseq.utils.file import file_exists, check_path, file_prefix
from hiseq.utils.utils import log, Config, run_shell_cmd, get_date, update_obj
from hiseq.utils.bam import Bam



class Bam2bwRn(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
    
    def init_args(self):
        args_init = {
            'bam': None,
            'scaleFactor': 1.0,
        }
        self = update_obj(self, args_init, force=False)
        if isinstance(self.bam, str):
            self.bam = [self.bam]
        if isinstance(self.scaleFactor, float):
            self.scaleFactor = [self.scaleFactor]*len(self.bam)
        elif isinstance(self.scaleFactor, list) and len(self.scaleFactor) == len(self.bam):
            pass # all float ?
        else:
            log.warning('illegal scaleFactor: {}'.format(self.scaleFactor))
            self.scaleFactor = [1.0]*len(self.bam)
        
        
    def run(self):
        for bam, sf in zip(self.bam, self.scaleFactor):
            bname = file_prefix(bam)
            subdir = os.path.join(self.outdir, bname)
            args_local = self.__dict__.copy()
            args_local.update({'bam': bam, 'scaleFactor': sf, 'outdir': subdir})
            Bam2bw(**args_local).run()


# solution-1: BAM -> bg -> bw
class Bam2bw(object):
    """
    Convert bam to bigWig
    version-1: bedtools genomecov -bdg , bedGraphToBigWig, ...
    version-2: deeptools, genomeCoverage

    Convert bam to bigWig using deeptools
    https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    history:
    1. Mappable sequence of a genome, see Table 1 in 
       url: https://www.nature.com/articles/nbt.1518.pdf
    2. effective genome size:
        - non-N bases
        - regions (of some size) uniquely mappable
    3. UCSC
    http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics
    http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics

    !!! strandness
    default: dUTP-based library (read2 is sense strand, read1 is anti-sense strand)
    general RNA library: (NSR, small-RNA-library), read1 is sense, read2 is antisense
    """
    def __init__(self, **kwargs):
        """
        bam, file to bam file
        outdir, str, path to save all files
        genome, dm3,dm6,...
        binsize, int, default: 50
        strandness, 0, [0|1|2|12], 1=fwd, 2=rev, 12=rev,rev, 0=non-strandness
        reference, str, fasta file
        overwrite, bool, False

        Required:
        bam, outdir, genome, binsize, strandness, overwrite, reference
        """
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        Config().dump(self.__dict__, self.config_yaml)


    def init_args(self):
        # genome size
        effsize = {
            'dm3': 162367812,
            'dm6': 142573017,
            'mm9': 2620345972,
            'mm10': 2652783500,
            'hg19': 2451960000,
            'hg38': 2913022398,
            'GRCh38': 2913022398}
        # gsize = effsize.get(genome, 0) 
        args_init = {
            'bam': None,
            'prefix': None,
            'outdir': str(pathlib.Path.cwd()),
            'binsize': 50,
            'strandness': 0,
            'overwrite': False,
            'reference': None,
            'genome': None,
            'genome_size': None,
            'scaleFactor': 1.0,
            'normalizeUsing': 'RPGC',
            'threads': 4,
            'config_txt': os.path.join(self.outdir, 'config.txt'),
            'config_yaml': os.path.join(self.outdir, 'config.yaml'),
            'flag': True # whether run/not
        }
        self = update_obj(self, args_init, force=False)
        self.bamcoverage = which('bamCoverage')
        if not file_exists(self.bam):
            self.flag = False
        Bam(self.bam).index()
#         self.prefix = os.path.splitext(os.path.basename(self.bam))[0]
        if not isinstance(self.prefix, str):
            self.prefix = file_prefix(self.bam)
        self.bw_fwd = os.path.join(self.outdir, self.prefix + '.fwd.bigWig')
        self.bw_rev = os.path.join(self.outdir, self.prefix + '.rev.bigWig')
        self.bw = os.path.join(self.outdir, self.prefix + '.bigWig')
        self.bw_fwd_log = os.path.join(self.outdir, self.prefix + '.fwd.deeptools.log')
        self.bw_rev_log = os.path.join(self.outdir, self.prefix + '.rev.deeptools.log')
        self.bw_log = os.path.join(self.outdir, self.prefix + '.deeptools.log')
        check_path(self.outdir, create_dirs=True)
        # binsize
        if not isinstance(self.binsize, int):
            log.warning('binsize: int expected. {} got, auto-set: 50'.format(self.binsize))
            self.binsize = 50
        # strand
        if not self.strandness in [0, 1, 2, 12]:
            log.warning('strandness: [0|1|2|12] expected, {} got, auto-set: 0'.format(self.strandness))
            self.strandness = 0
        # overwrite
        if not isinstance(self.overwrite, bool):
            log.warning('overwrite: bool expected, {} got, auto-set: False'.format(self.overwrite))
            self.overwrite = False
        # genome/genome_size/reference
        if self.genome in effsize:
            self.genome_size = effsize.get(self.genome, None)
        elif isinstance(self.genome_size, int):
            pass
        elif not self.reference is None:
            self.genome_size = fasize(self.reference, sum_all=True)
        else:
            log.error('genome, genome_size, reference; one of arg is requred')
            self.flag = False # do not run


    def run_cmd(self, cmd, bw, bw_log):
        if os.path.exists(bw) and not self.overwrite:
            log.info('Bam2bw() skipped, file exists, ..., {}'.format(bw))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(bw_log, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')
                # p1 = subprocess.run(shlex.split(cmd), stdout=w, stderr=w)
        if not os.path.exists(bw):
            log.error('Bam2bw() failed, bw not found: {}'.format(self.bw))
            # raise ValueError('output file is missing, check log file: %s' % bw_log)


    def run_fwd(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_fwd),
            '--binSize {}'.format(self.binsize), 
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--filterRNAstrand forward',
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing {}'.format(self.normalizeUsing),
            '2> {}'.format(self.bw_log)
            ])
        cmd_txt = os.path.join(self.outdir, 'cmd_fwd.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        self.run_cmd(cmd, self.bw_fwd, self.bw_fwd_log)


    def run_rev(self):
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_rev),
            '--binSize {}'.format(self.binsize), 
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--filterRNAstrand reverse',
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing {}'.format(self.normalizeUsing)])
        cmd_txt = os.path.join(self.outdir, 'cmd_rev.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + 'n')

        self.run_cmd(cmd, self.bw_rev, self.bw_rev_log)
        

    def run_non(self):
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw),
            '--binSize {}'.format(self.binsize),
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing  {}'.format(self.normalizeUsing)])
        cmd_txt = os.path.join(self.outdir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        self.run_cmd(cmd, self.bw, self.bw_log)
        

    def run(self):
        """
        determined by: strandness
        exit by: flag
        """
        if self.flag is False:
            log.error('Bam2bw() failed, ...')
        elif self.strandness == 1:
            self.run_fwd()
        elif self.strandness == 2:
            self.run_rev()
        elif self.strandness == 12:
            self.run_fwd()
            self.run_rev()
        elif self.strandness == 0:
            self.run_non()
        else:
            log.error('Bam2bw() failed, unknown strandness')


def bw_compare(bw1, bw2, bw_out, operation='log2', **kwargs):
    """
    Compare two bigWig files: ip over input

    example:
    bigwigCompare -b1 bw1 -b2 bw2 --operation 
    {log2, ratio, subtract, add, mean, reciprocal_ratio, 
    first, second}
    -o out.bw
    """
    binsize = kwargs.get('binsize', 50)
    threads = kwargs.get('threads', 1)
    overwrite = kwargs.get('overwrite', False)
    # outdir
    bw_out_dir = os.path.dirname(bw_out)
    bw_stdout = os.path.join(bw_out_dir, 'bw_compare.stdout')
    bw_stderr = os.path.join(bw_out_dir, 'bw_compare.stderr')
    check_path(bw_out_dir, create_dirs=True)
    # cmd
    cmd = ' '.join([
        '{}'.format(which('bigwigCompare')),
        '--bigwig1 {} --bigwig2 {}'.format(bw1, bw2),
        '--operation {}'.format(operation),
        '--skipZeroOverZero',
        '--skipNAs',
        '--binSize {}'.format(binsize),
        '-p {}'.format(threads),
        '-o {}'.format(bw_out),
        '1> {}'.format(bw_stdout),
        '2> {}'.format(bw_stderr),
    ])
    cmd_txt = os.path.join(bw_out_dir, 'cmd.sh')
    with open(cmd_txt, 'wt') as w:
        w.write(cmd + '\n')
    # run
    if os.path.exists(bw_out) and not overwrite:
        log.info('bwCompare() skipped, file exists: {}'.format(bw_out))
    else:
        run_shell_cmd(cmd)



def get_args():
    parser = argparse.ArgumentParser(description='hiseq bam2bw -i bam -g dm6 -o outdir')
    parser.add_argument('-i', '--bam', nargs='+', required=True,
        help='BAM files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-s', '--strandness', default=0, type=int,
        help='the strandness, 0=no, 1=fwd, 2=rev, 12=fwd&rev, default: [0]')
    parser.add_argument('-b', '--binsize', default=50, type=int,
        help='set binSize for bigWig, default: [50]')
    parser.add_argument('-g', '--genome', default=None,
        help='choose genome for the bam file, default: [None]')
    parser.add_argument('--scaleFactor', nargs='+', type=float, default=1.0,
        help='The scaling factor, default: [1.0]')
    parser.add_argument('--normalizeUsing', default='None',
        choices=['RPGC', 'RPKM', 'CPM', 'BPM', 'None'],
        help='Possible choices: RPKM, CPM, BPM, RPGC, None, default: [None] \
        see https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html, \
        for details from deeptools documentation.')
    parser.add_argument('-gs', '--genome-size', dest='genome_size', default=None,
        help='set the genome size for input genome, default: [None]' )
    parser.add_argument('-r', '--reference', default=None,
        help='the reference genome in fasta format, default: [None]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Whether overwrite exists files')
    parser.add_argument('-p', '--threads', type=int, default=4,
        help='Number of threads, default: [4]')
    return parser


def main():
    args = vars(get_args().parse_args())
    Bam2bwRn(**args).run()


if __name__ == '__main__':
    main()

#
        
        
        
        
# class Bam2bw2(object):
#     """
#     Convert bam to bigWig
#     version-2: bedtools genomecov -bdg , bedGraphToBigWig, ...

#     !!! strandness
#     default: dUTP-based library (read2 is sense strand, read1 is anti-sense strand)
#     general RNA library: (NSR, small-RNA-library), read1 is sense, read2 is antisense
#     """
#     def __init__(self, **kwargs):
#         """
#         bam, file to bam file
#         outdir, str, path to save all files
#         genome, dm3,dm6,...
#         binsize, int, default: 50
#         strandness, 0, [0|1|2|12], 1=fwd, 2=rev, 12=rev,rev, 0=non-strandness
#         reference, str, fasta file
#         overwrite, bool, False

#         Required:
#         bam, outdir, genome, binsize, strandness, overwrite, reference
#         """
#         # update self (obj)
#         for k, v in kwargs.items():
#             setattr(self, k, v)

#         # update, defaults
#         self.init_args()
#         Config().dump(self.__dict__, self.config_yaml)


#     def init_args(self):
#         args_init = {
#             'bam': None,
#             'outdir': str(pathlib.Path.cwd()),
#             'binsize': 50,
#             'strandness': 0,
#             'overwrite': False,
#             'reference': None,
#             'genome': None,
#             'genome_size': None,
#             'norm_size': 1000000, # 1 million
#             'fragment': True,
#             'threads': 4,
#             'config_txt': os.path.join(self.outdir, 'config.txt'),
#             'config_yaml': os.path.join(self.outdir, 'config.yaml'),
#             'flag': True # whether run/not
#         }
#         self = update_obj(self, args_init, force=False)
#         self.bamcoverage = which('bamCoverage')
#         if not file_exists(self.bam):
#             self.flag = False
#         Bam(self.bam).index()
#         self.prefix = os.path.splitext(os.path.basename(self.bam))[0]
#         self.bw_fwd = os.path.join(self.outdir, self.prefix + '.fwd.bigWig')
#         self.bw_rev = os.path.join(self.outdir, self.prefix + '.rev.bigWig')
#         self.bw = os.path.join(self.outdir, self.prefix + '.bigWig')
#         self.bw_fwd_log = os.path.join(self.outdir, self.prefix + '.fwd.deeptools.log')
#         self.bw_rev_log = os.path.join(self.outdir, self.prefix + '.rev.deeptools.log')
#         self.bw_log = os.path.join(self.outdir, self.prefix + '.deeptools.log')
#         check_path(self.outdir, create_dirs=True)
#         if not isinstance(self.binsize, int):
#             log.warning('binsize: int expected. {} got, auto-set: 50'.format(self.binsize))
#             self.binsize = 50
#         # norm size
#         # default: 1 million
#         if isinstance(self.norm_size, int):
#             if self.norm_size < 1:
#                 self.norm_size = 1000000 # 1 million
#         else:
#             self.norm_size = 1000000 # 1 million
#         # strand
#         if not self.strandness in [0, 1, 2, 12]:
#             log.warning('strandness: [0|1|2|12] expected, {} got, auto-set: 0'.format(self.strandness))
#             self.strandness = 0
#         # overwrite
#         if not isinstance(self.overwrite, bool):
#             log.warning('overwrite: bool expected, {} got, auto-set: False'.format(self.overwrite))
#             self.overwrite = False
#         # genome/genome_size/reference
#         if self.genome in effsize:
#             self.genome_size = effsize.get(self.genome, None)
#         elif isinstance(self.genome_size, int):
#             pass
#         elif not self.reference is None:
#             self.genome_size = fasize(self.reference, sum_all=True)
#         else:
#             log.error('genome, genome_size, reference; one of arg is requred')
#             self.flag = False # do not run


#     def get_reads(self, fragment=False):
#         """
#         fragment, pair of PE reads
#         normalize to 1 million mapped reads
#         RPM, CPM
#         """
#         b = Bam(self.bam)
#         b_is_paired = b.isPaired()
#         b_reads = b.getNumReads()
#         if fragment and b_is_paired:
#             b_frag = b_reads / 2
#         else:
#             b_frag = b_reads

#         return b_frag


#     def get_scale(self):
#         """
#         Norm to 1 million reads
#         """
#         r = self.get_reads(fragment=self.fragment)
#         # norm by reads/fragment
#         return r/self.norm_size


#     def run_cmd(self, cmd, bw, bw_log):
#         if os.path.exists(bw) and not self.overwrite:
#             log.info('file exists, bigWig skipped ..., {}'.format(bw))
#         else:
#             _, stdout, stderr = run_shell_cmd(cmd)
#             with open(bw_log, 'wt') as w:
#                 w.write(stdout + '\n' + stderr + '\n')
#                 # p1 = subprocess.run(shlex.split(cmd), stdout=w, stderr=w)
#         if not os.path.exists(bw):
#             log.error('Bam2bw() failed, bw not found: {}'.format(self.bw))
#             # raise ValueError('output file is missing, check log file: %s' % bw_log)


#     def run_fwd(self):
#         # bedtools cmd
#         cmd = ' '.join([
#             '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_fwd),
#             '--binSize {}'.format(self.binsize), 
#             '--effectiveGenomeSize {}'.format(self.genome_size),
#             '--filterRNAstrand forward',
#             '--numberOfProcessors {}'.format(self.threads), 
#             '--scaleFactor {}'.format(self.scaleFactor),
#             '--normalizeUsing {}'.format(self.normalizeUsing)])
#         cmd_txt = os.path.join(self.outdir, 'cmd_fwd.txt')
#         with open(cmd_txt, 'wt') as w:
#             w.write(cmd + '\n')
#         self.run_cmd(cmd, self.bw_fwd, self.bw_fwd_log)


#     def run_rev(self):
#         # bedtools cmd
#         cmd = ' '.join([
#             '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_rev),
#             '--binSize {}'.format(self.binsize), 
#             '--effectiveGenomeSize {}'.format(self.genome_size),
#             '--filterRNAstrand reverse',
#             '--numberOfProcessors {}'.format(self.threads), 
#             '--scaleFactor {}'.format(self.scaleFactor),
#             '--normalizeUsing {}'.format(self.normalizeUsing)])
#         cmd_txt = os.path.join(self.outdir, 'cmd_rev.txt')
#         with open(cmd_txt, 'wt') as w:
#             w.write(cmd + 'n')
#         self.run_cmd(cmd, self.bw_rev, self.bw_rev_log)
        

#     def run_non(self):
#         # bedtools cmd
#         cmd = ' '.join([
#             '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw),
#             '--binSize {}'.format(self.binsize),
#             '--effectiveGenomeSize {}'.format(self.genome_size),
#             '--numberOfProcessors {}'.format(self.threads), 
#             '--scaleFactor {}'.format(self.scaleFactor),
#             '--normalizeUsing  {}'.format(self.normalizeUsing)])
#         cmd_txt = os.path.join(self.outdir, 'cmd.txt')
#         with open(cmd_txt, 'wt') as w:
#             w.write(cmd + '\n')
#         self.run_cmd(cmd, self.bw, self.bw_log)
        

#     def run(self):
#         """
#         determined by: strandness
#         exit by: flag
#         """
#         if self.flag is False:
#             log.error('bam2bw skipped, ...')
#         elif self.strandness == 1:
#             self.run_fwd()
#         elif self.strandness == 2:
#             self.run_rev()
#         elif self.strandness == 12:
#             self.run_fwd()
#             self.run_rev()
#         elif self.strandness == 0:
#             self.run_non()
#         else:
#             log.error('unknown strandness')