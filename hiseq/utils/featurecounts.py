#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for FeatureCounts

1. run shell command 
2. read output
"""

import os
import pandas as pd
import pybedtools
from hiseq.utils.helper import *
from .bed import *
from .bam import *


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
    """Read the count.txt.summary
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
        print(assign_df)
    except:
        log.warning('reading file failed: {}'.format(self.summary))
        total, assign, assign_pct = [1, 0, 0]
    df = pd.DataFrame([total, assign, assign_pct]).T
    df.columns = ['total', 'map', 'pct']
    return df


class FeatureCounts(object):
    """Run featureCounts for {BED|GTF} and BAM(s)
    
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
        self.init_files()


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
        if not isinstance(self.gtf, str):
            raise ValueError('gtf=, expect str, got {}'.format(
                type(self.gtf).__name__))
        if not file_exists(self.gtf):
            raise ValueError('gtf=, file not exists: {}'.format(self.gtf))
        if isinstance(self.bam_list, str):
            self.bam_list = [self.bam_list]
        elif isinstance(self.bam_list, list):
            pass
        else:
            raise ValueError('bam_list, expect str or list, got {}'.format(
                type(self.bam_list).__name__))
        if not all(file_exists(self.bam_list)):
            raise ValueError('bam-list, file not exists:')
        # [Bam(i).index() for i in self.bam_list]
        if not isinstance(self.outdir, str):
            self.outdir = self._tmp(dir=True)
        # check_path(self.outdir)
        if not isinstance(self.prefix, str):
            self.prefix = 'count.txt'
        self.gtf = file_abspath(self.gtf)
        self.bam_list = file_abspath(self.bam_list)
        self.outdir = file_abspath(self.outdir)
        # strand
        if not self.strandness in [0, 1, 2]:
            raise ValueError('strandness=, not valid, expect [0, 1, 2], \
                got {}'.format(strandness))


    def init_files(self):
        config_dir = os.path.join(self.outdir, 'config')
        prefix = os.path.join(self.outdir, self.prefix)
        default_files = {
            'config_toml': config_dir + '/config.toml',
            'count_txt': prefix,
            'summary': prefix + '.summary',
            'log': prefix + '.featureCounts.log',
            'stat': prefix + '.featureCounts.stat',
            'cmd_shell': prefix + '.cmd.sh'
        }
        self.config_dir = config_dir
        self = update_obj(self, default_files, force=True) # key
        # check_path(self.config_dir)
        self.gtf_ext = os.path.splitext(self.gtf)[1].lower()
        if self.gtf_ext in ['.bed', '.narrowpeak', '.broadpeak']:
            self.saf = os.path.join(self.outdir, file_prefix(self.gtf)[0] + '.saf')
            bed_to_saf(self.gtf, self.saf)


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
            arg = '-a {} -F GTF -t {} -g gene_id '.format(self.gtf, self.feature_type)
        elif self.gtf_ext in ['.bed', '.narrowpeak', '.broadpeak']:
            arg = '-a {} -F SAF'.format(self.saf) # SAF
        elif self.gtf_ext in ['.saf']:
            arg = '-a {} -F SAF'.format(self.gtf)
        else:
            arg = ''
        return arg


    def get_cmd(self):
        """
        prepare args for featureCounts
        """
        args_gtf = self.get_args_gtf()
        args_pe = '-p -C -B' if self.is_pe() else ''
        return ' '.join([
            '{}'.format(shutil.which('featureCounts')),
            '-s {}'.format(self.strandness),
            args_gtf,
            args_pe,
            '-o {}'.format(self.count_txt),
            '-T {}'.format(self.threads),
            '-M -O --fraction',
            ' '.join(self.bam_list),
            '2> {}'.format(self.log)
            ])


    def pre_run(self):
        # save config
        check_path([self.config_dir])
        Config().to_toml(self.__dict__, self.config_toml)
        # show message
        msg = '\n'.join([
            '-'*80,
            'Run featureCounts:',
            '{:>10s}: {}'.format('bam_files', self.bam_list),
            '{:>10s}: {}'.format('gtf file', self.gtf),
            '{:>10s}: {}'.format('feature', self.feature_type),
            '{:>10s}: {}'.format('strand', self.strandness),
            '{:>10s}: {}'.format('count.txt', self.count_txt),
            '{:>10s}: {}'.format('log', self.log),
            '-'*80,
        ])
        print(msg)


    def run(self):
        """Run featureCounts
        """
        self.pre_run()
        # run all
        cmd = self.get_cmd()
        with open(self.cmd_shell, 'wt') as w:
            w.write(cmd + '\n')
        if file_exists(self.count_txt) and not self.overwrite:
            log.info("FeatureCounts() skippped, file exists: {}".format(
                self.count_txt))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('FeatureCounts() failed, see: {}'.format(self.log))
        return read_fc_summary(self.summary)

