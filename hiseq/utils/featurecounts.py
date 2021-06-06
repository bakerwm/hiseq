#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for FeatureCounts

1. run shell command 
2. read output
"""

import os
import tempfile
import pybedtools
import pandas as pd
from shutil import which
from hiseq.utils.bam import Bam
from hiseq.utils.file import file_exists, file_abspath, file_prefix, \
    check_path, remove_path
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd


def bed_to_gtf(file_in, file_out):
    """Convert BED to GTF
    chrom chromStart chromEnd name score strand
    """
    if file_exists(file_out) and not overwrite:
        log.info('bed_to_gtf() skipped, file exists: {}'.format(file_out))
    else:
        try:
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
        except:
            log.error('bed_to_gtf() failed, see: {}'.format(file_out))
    return file_out


def bed_to_saf(file_in, file_out, overwrite=False):
    """Convert BED to SAF format, for featureCounts
    GeneID Chr Start End Strand

    see: https://www.biostars.org/p/228636/#319624
    """
    if file_exists(file_out) and not overwrite:
        log.info('bed_to_saf() skipped, file exists: {}'.format(file_out))
    else:
        try:
            with open(file_in, 'rt') as r, open(file_out, 'wt') as w:
                for line in r:
                    tabs = line.strip().split('\t')
                    chr, start, end = tabs[:3]
                    # strand
                    strand = '.'
                    if len(tabs) > 5:
                        s = tabs[5]
                        if s in ['+', '-', '.']:
                            strand = s
                    # id
                    if len(tabs) > 3:
                        name = os.path.basename(tabs[3])
                    else:
                        name = '_'.join([chr, start, end, strand])
                    # output
                    w.write('\t'.join([name, chr, start, end, strand])+'\n')
        except:
            log.error('bed_to_saf() failed, see: {}'.format(file_out))
    return file_out


def read_fc_txt(x, fix_name=False):
    """Read the count.txt
    fix the name of bam file
    """
    try:
        df = pd.read_csv(x, '\t', comment='#')
        if fix_name:
            bam_list = df.columns.to_list()[6:] # from 7-column to end
            bam_list = [os.path.splitext(i)[0] for i in bam_list]
            out = [os.path.basename(i) for i in bam_list]
        else:
            out = df
    except:
        log.error('reading file failed, {}'.format(x))
        out = None
    return out


def read_fc_summary(x):
    """
    Read the count.txt.summary
    output: dict
    {index:, total:, map:, pct:}
    """
    if not file_exists(x):
        log.error('file not exists, {}'.format(x))
        return None
    # parse the strand from count.txt file, line-1
    count_txt = os.path.splitext(x)[0]
    if not file_exists(count_txt):
        log.error('count.txt not found, {}'.format(count_txt))
        return None
    # parse arguments
    with open(count_txt) as r:
        cmd_line = r.readline()
    args = cmd_line.replace('"', '').split()
    if '-s' in args:
        s = args[args.index('-s')+1]
    else:
        s = '0' # default
    # read summary
    try:
        df = pd.read_csv(x, '\t', index_col=0)
        df.columns = list(map(os.path.basename, df.columns.to_list()))
        total = df.sum(axis=0, skipna=True)
        assign = df.loc['Assigned', ]
        assign_pct = assign / total
        assign_pct = assign_pct.round(decimals=4)
        assign_df = assign_pct.to_frame('assigned')
        assign_df['strandness'] = s
        # mimimal value
        assign_min = assign_pct.min()
        if assign_min < 0.50:
            log.warning('Caution: -s {}, {:.2f}% assigned, see {}'.format(
                s, assign_min, x))
        log.info(assign_df)
    except:
        log.warning('reading file failed: {}'.format(self.summary))
        total, assign, assign_pct = [1, 0, 0]
    df = pd.DataFrame([total, assign, assign_pct]).T
    df.columns = ['total', 'map', 'pct']
    return df


class FeatureCounts(object):
    """
    Run featureCounts for {BED|GTF} and BAM(s)
    
    Keyword parameters
    ------------------
    gtf : str
        Path to the gtf, require (features:exon, gene)
        
    bam_list : str or list
        List of bam files
    
    outdir : str
        The directory saving the results
        
    strandness : int or None
        The strandness, 0=no, 1=sens, 2=anti, None=no, default: [None]
        
    prefix : str
        The prefix of the output file, default: [count.txt] 
        
    threads : int
        Number of threads, default: [4] 
        
    overwrite : bool
        Overwrite the exists file,
        
    feature_type : str
        Specify the feature type in GTF annotation, 'exon', 'gene', 
        
    example:
    FeatureCounts(gtf=a, bam_list=b, outdir=c).run()
    gene file format: bed, gtf, saf
    
    bed_to_saf
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'gtf': None,
            'bam_list': None,
            'outdir': None,
            'prefix': None,
            'strandness': 0,
            'threads': 4,
            'overwrite': False,
            'feature_type': 'exon',
        }
        self = update_obj(self, args_init, force=False)
        if not isinstance(self.outdir, str):
            self.outdir = self._tmp(dir=True)
        self.outdir = file_abspath(self.outdir)
        if not isinstance(self.prefix, str):
            self.prefix = 'count.txt'
        self.gtf = file_abspath(self.gtf)
        self.bam_list = file_abspath(self.bam_list)
        if not self.strandness in [0, 1, 2]:
            raise ValueError('strandness=, not valid, expect [0, 1, 2], \
                got {}'.format(strandness))
        # update files
        self.init_files()
        self.init_gtf()
        # dirs
        check_path(self.config_dir, create_dirs=True)
        Config().dump(self.__dict__, self.config_yaml)


    def init_files(self):
        self.config_dir = os.path.join(self.outdir, 'config')
        prefix = os.path.join(self.outdir, self.prefix)
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
            'count_txt': prefix,
            'summary': prefix + '.summary',
            'summary_json': prefix + '.summary.json',
            'log_stdout': prefix + '.featureCounts.stdout',
            'log_stderr': prefix + '.featureCounts.stderr',
            'stat': prefix + '.featureCounts.stat',
            'cmd_shell': prefix + '.cmd.sh',
            'saf': os.path.join(self.outdir, file_prefix(self.gtf) + '.saf')
        }
        self = update_obj(self, default_files, force=True) # key


    def init_gtf(self):
        if isinstance(self.bam_list, str):
            self.bam_list = [self.bam_list]
        msg = '\n'.join([
            '='*80,
            'Check files for FeatureCounts',
            '{:>14s} : {}'.format('gtf', self.gtf),
            '{:>14s} : {}'.format('bam_list', self.bam_list),
            '='*80,
        ])
        print(msg)
        if not all([
            isinstance(self.gtf, str),
            isinstance(self.bam_list, list),
            file_exists(self.gtf),
            all(file_exists(self.bam_list)),
        ]):
            raise ValueError('FeatureCounts() failed, check above message.')
        # check gtf
        self.gtf_ext = os.path.splitext(self.gtf)[1].lower()
        self.is_paired = self.is_pe()


    def _tmp(self, dir=False):
        """Create a tmp file
        """
        if dir:
            tmp = tempfile.TemporaryDirectory()
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=False)
        return tmp.name


    def is_pe(self):
        """Check whether the input bam files are Paired or Single file
        Bam().isPaired()
        """
        return all([Bam(i).is_paired() for i in self.bam_list])


    def get_args_gtf(self):
        """
        gtf
        saf
        bed -> saf
        """
        # check gtf or saf or bed
        if self.gtf_ext in ['.gtf']:
            s = '-a {} -F GTF -t {} -g gene_id '.format(
                self.gtf, self.feature_type)
        elif self.gtf_ext in ['.bed', '.narrowpeak', '.broadpeak']:
            bed_to_saf(self.gtf, self.saf)
            s = '-a {} -F SAF'.format(self.saf) # SAF
        elif self.gtf_ext in ['.saf']:
            s = '-a {} -F SAF'.format(self.gtf)
        else:
            raise ValueError('unknown file, expect *.gtf, got {}'.format(
                self.gtf))
            s = ''
        return s


    def show_msg(self):
        msg = '\n'.join([
            '-'*80,
            'Run featureCounts:',
            '{:>10s}: {}'.format('bam_files', self.bam_list),
            '{:>10s}: {}'.format('gtf file', self.gtf),
            '{:>10s}: {}'.format('feature', self.feature_type),
            '{:>10s}: {}'.format('strand', self.strandness),
            '{:>10s}: {}'.format('count.txt', self.count_txt),
            '{:>10s}: {}'.format('log_stdout', self.log_stdout),
            '{:>10s}: {}'.format('log_stderr', self.log_stderr),
            '-'*80,
        ])
        print(msg)


    def run_featurecounts(self):
        """
        prepare args for featureCounts
        """
        self.args_gtf = self.get_args_gtf()
        self.args_bam = ' '.join(self.bam_list)
        self.args_pe = '-p -C -B' if self.is_paired else ''
        cmd = ' '.join([
            '{}'.format(which('featureCounts')),
            '-s {}'.format(self.strandness),
            '-o {}'.format(self.count_txt),
            '-T {}'.format(self.threads),
            '-M -O --fraction',
            self.args_pe,
            self.args_gtf,
            self.args_bam,
            '1> {}'.format(self.log_stdout),
            '2> {}'.format(self.log_stderr),
            ])
        # save command
        with open(self.cmd_shell, 'wt') as w:
            w.write(cmd + '\n')
        if file_exists(self.count_txt) and not self.overwrite:
            log.info("FeatureCounts() skippped, file exists: {}".format(
                self.count_txt))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('FeatureCounts() failed, see: {}'.format(
                    self.log_stderr))


    def wrap_summary(self):
        # !!! single or multiple bam files ?! multiple
        df = read_fc_summary(self.summary) # pd.DataFrame
        df['filename'] = df.index # add filename to column
        d = df.transpose().to_dict('dict') # to dict
        # df2 = df.reset_index().to_dict('dict')
        Config().dump(d, self.summary_json)
        return d
        

    def run(self):
        """Run featureCounts
        """
        self.show_msg()
        self.run_featurecounts()
        return self.wrap_summary()

    
class LibStrand(object):
    """
    Run featureCounts for {BED|GTF} and BAM(s)
    
    Parameters
    ----------
    bam : str
        Single bam file
        
    gtf : str
        Path to the gene annotation, GTF/BED format
    
    outdir : str
        The directory saving the results
        
    size : int
        Number of reads subset from BAM file, default [200000]
        
    clean_up : bool
        Remove the temp files, default [True]
        
    Guess the library type of the HiSeq data
    Strandness
    strandness: 1 ++, 1 --, / 2 +-, 2 -+ :
    dUTP, NSR: 1 +-, 1 -+, / 2 ++, 2 -- :

    see: infer_experiment.py from RSeQC package.
    
    Examples:
    >>> s = LibStrand(bam='in.bam', gtf='gene.gtf')
    >>> s.sense_strand # 
    >>> s.is_strand_specific # 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.run()
        
 
    def init_args(self):
        args_init = {
            'bam': None,
            'gtf': None,
            'size': 200000,
            'outdir': None,
            'clean_up': True,
        }
        self = update_obj(self, args_init, force=False)
        # output, tmp dir
        if self.outdir is None:
            self.outdir = self._tmp(is_dir=True)
        check_path(self.outdir)
        self.init_files()


    def _tmp(self, is_dir=False, suffix='.txt'):
        """
        Create a tmp file to save json object
        """
        if is_dir:
            tmp = tempfile.TemporaryDirectory(prefix='tmp')
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
                delete=False)
        return tmp.name

    
    def init_files(self):
        c1 = isinstance(self.bam, str)
        c1e = file_exists(self.bam)
        c2 = isinstance(self.gtf, str)
        c2e = file_exists(self.gtf)
        self.file_is_ok = all([c1, c1e, c2, c2e])
        if self.file_is_ok:
            self.bam_sub = Bam(self.bam).subset(self.size, self.outdir)
            self.bam_name = os.path.basename(self.bam)
        else:
            self.bam_sub = None
            self.bam_name = None
            log.error('bam, expect str, got {}'.format(self.bam))
            log.error('gtf, expect str, got {}'.format(self.gtf))


    def run_fc(self, strandness=1):
        """
        strandness; 1 or 2; for featureCounts
        using -s 1, -s 2
        and return the mapping reads
        """
        prefix = 'count.s_{}.txt'.format(strandness)
        args = {
            'gtf': self.gtf,
            'bam_list': self.bam_sub,
            'outdir': self.outdir,
            'prefix': prefix,
            'strandness': strandness
        }
        return FeatureCounts(**args).run()
        

    def guess_sense(self, with_status=False):
        """
        Return the sense strand, s=
        """
        msg = '\n'.join([
            '='*80,
            'Guess strandness of BAM file, by featureCounts',
            'Library strandness: -s=1 or -s=2',
            'bam: {}'.format(self.bam),
            'gtf: {}'.format(self.gtf),
            '='*80,
        ])
        print(msg)
        
        if not self.file_is_ok:
            return None # skipped
        # run fc: s=1
        df1 = self.run_fc(strandness=1) # dict {}
        df1b = df1.get(self.bam_name, {})
        pct1 = df1b.get('pct', 0)
        # run fc: s=2
        df2 = self.run_fc(strandness=2) # dict
        df2b = df2.get(self.bam_name, {})
        pct2 = df2b.get('pct', 0)
        # check strand:
        if pct1 + pct2 < 0.2:
            log.warning('Less than 20% reads assigned, {}'.format(self.gtf))
            s = 0
        elif pct1 == pct2:
            log.info('Library type (s=0): non-stranded')
            s = 0
        elif pct1 > pct2:
            log.info('Library type (s=1): 1 ++, 1 --, / 2 +-, 2 -+')
            s = 1
        else:
            log.info('Library type (s=2): 1 +-, 1 -+, / 2 ++, 2 --')
            s = 2
        # clean up
        if self.clean_up:
            remove_path(self.outdir, ask=False)
        # log
        msg = '\n'.join([
            '='*80,
            'Status:',
            'forward stranded, -s=1, {:>6.2f}%'.format(pct1*100),
            'reverse stranded, -s=2, {:>6.2f}%'.format(pct2*100),
            '='*80,
        ])
        print(msg)
        return s # sense
    

    def run(self):
        s = self.guess_sense() # 0, 1, 2, None
        self.is_strand_specific = False
        self.read1_is_sense = True # non-stranded
        self.read2_is_sense = True # non-stranded
        self.sense_strand = s #
        if s == 1:
            self.anti_strand = 2
        elif s == 2:
            self.anti_strand = 1
        else:
            self.anti_strand = s
        if s is None:
            log.error('LibStrand() failed, check bam/gtf')
        elif s == 0:
            log.info('non-stranded BAM file: {}'.format(self.bam))
        elif s > 0:
            self.is_strand_specific = True
            self.read1_is_sense = s == 1
            self.read2_is_sense = s == 2
