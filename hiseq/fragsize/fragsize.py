#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate the fragment size of PE reads,

Available tools: CollectInsertSizeMetrics (Picard), bamPEFragmentSize (deeptools)

Here are the samtools method:

Statistics:
1. length count
2. mean/medium/qrt.
"""

import os
import sys
import pathlib
import numpy as np
import pandas as pd
import pysam
import tempfile
import logging
import argparse
import hiseq
from multiprocessing import Pool
from hiseq.utils.utils import update_obj, log, get_date, run_shell_cmd
from hiseq.utils.file import file_exists, file_prefix, check_path
from hiseq.utils.bam import Bam


class BamFragSize(object):
    """Calculate the read size of SE, multiple BAM
    multiple BAM files
    rename (labels)
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args() # update, default args


    def init_args(self):
        args_init = {
            'bam': None,
            'outdir': None,
            'labels': None,
            'as_se': False,
            'max_count': 0,
            'strandness': False,
            'csv_file': None,
            'plot_file': None,
            'threads': 4,
            'parallel_jobs': 1,
            'overwrite': False
        }
        self = update_obj(self, args_init, force=False)
        # bam
        if isinstance(self.bam, str):
            self.bam = [self.bam]
        elif isinstance(self.bam, list):
            pass
        else:
            raise ValueError('bam, str or list expected, got {}'.format(
                type(self.bam).__name__))
        if isinstance(self.labels, str):
            self.labels = [self.labels]
        elif isinstance(self.labels, list):
            pass
        else:
            self.labels = None
        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        # labels
        bam_names = file_prefix(self.bam)
        if self.labels is None:
            self.labels = bam_names
        else:
            self.labels = self.labels[:len(bam_names)] + bam_names[len(self.labels):]


    def run_single_bam(self, bam):
        """Run for single bam
        stat
        pdf
        require outdir
        """
        bam_i = self.bam.index(bam) # index for bam file
        labels = self.labels[bam_i] # labels for bam
        csv_file = os.path.join(self.outdir, labels + '.fragsize.csv')
        if os.path.exists(csv_file) and not self.overwrite:
            log.info('BamFragSize() skipped, file exists: {}'.format(csv_file))
        else:
            args_local = {
                'bam': bam,
                'outdir': self.outdir,
                'labels': labels,
                'as_se': self.as_se,
                'strandness': self.strandness,
            }
            BamFragSizeR1(**args_local).run()
            # BamPEFragSize(**args_local).plot(pdf_file)


    def run_multiple_bam(self):
        if len(self.bam) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.threads) as pool:
                pool.map(self.run_single_bam, self.bam)
        else:
            [self.run_single_bam(i) for i in self.bam]
        # plot merge
        csv_file_list = [os.path.join(self.outdir, i + '.fragsize.csv') for i in self.labels]
        csv_file_list = [i for i in csv_file_list if os.path.exists(i)]
        pdf_merge = os.path.join(self.outdir, 'frag_size.pdf')
        # save to plot
        # to-do: matplotlib function
        hiseq_dir = os.path.dirname(hiseq.__file__)
        frag_plot_r = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
        stdout = os.path.join(self.outdir, 'frag_plot.plot.stdout')
        stderr = os.path.join(self.outdir, 'frag_plot.plot.stderr')
        cmd_file = os.path.join(self.outdir, 'frag_plot.plot.cmd.sh')
        cmd = ' '.join([
            'Rscript',
            frag_plot_r,
            pdf_merge,
            ' '.join(csv_file_list),
            '1> {}'.format(stdout),
            '2> {}'.format(stderr)
        ])
        with open(cmd_file, 'wt') as w:
            w.write(cmd+'\n')
        try:
            run_shell_cmd(cmd)
        except:
            log.warning('fragsize.py failed')


    def run(self):
        self.run_multiple_bam()



class BamFragSizeR1(object):
    """Calculate the read size of SE, single BAM

    Parameters
    ----------
    bam: str
        bam file, single

    labels: str
        label

    as_se: bool
        Treat BAM as SE reads

    max_count: int
        maximum number of reads to process, default: [0], all

    strandness: bool
        Check the strandness, default: [False]

    csv_file: str
        File to save the results, in csv format

    sample size = 1000 (SE or PE)

    > Table.csv
    length strand count
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'bam': None,
            'outdir': None,
            'labels': None,
            'as_se': False,
            'max_count': 0,
            'strandness': False,
            'csv_file': None,
            'plot_file': None,
        }
        self = update_obj(self, args_init, force=False)
        # only for single bam file
        if not isinstance(self.bam, str):
            raise ValueError('str expected, {} got'.format(
                type(self.bam).__name__))
        if not self.is_bam(self.bam):
            raise Exception('bam={}, is not bam file'.format(self.bam))
        # bai
        if not file_exists(self.bam+'.bai'):
            Bam(self.bam).index()
        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        check_path(self.outdir)
        # check strandness
        # as_se: paired end only
        if not (self.is_paired(self.bam) and not self.as_se):
            self.as_se = True
        if not self.labels:
            self.labels = file_prefix(self.bam)
        # output files
        if self.csv_file is None:
            self.csv_file = os.path.join(self.outdir, self.labels + '.fragsize.csv')
        if self.plot_file is None:
            self.plot_file = os.path.join(self.outdir, self.labels + '.fragsize.pdf')


    def is_empty(self, bam):
        """Check input file is empty or not
        pysam.Samfile().count
        pysam.Samfile().mapped
        """
        pysam.index(bam)
        sam = pysam.Samfile(bam)
        return sam.count() == 0 and sam.mapped == 0


    def is_bam(self, bam):
        """Check input is BAM file
        update: self.bam

        string
        *.bam
        """
        if bam is None:
            bam = self.bam
        out = False
        if isinstance(self.bam, str):
            out = self.bam.endswith('.bam') and os.path.exists(self.bam)
        return out


    def is_paired(self, bam, topn=1000):
        """Check input bam is Paired end alignment"""
        out = False
        if self.is_bam(bam):
            samfile = pysam.AlignmentFile(bam)
            out = all([read.is_paired for read in samfile.head(topn)])
            samfile.close()
        return out


    def cal_freq(self, x, return_dataframe=True):
        """Calculate the frequency of list
        ['length', 'strand']
        return dataframe
        """
        header = ['length', 'strand', 'count']
        if isinstance(x, list):
            df = pd.DataFrame(x, columns=header).groupby(['length', 'strand']).count().reset_index()
        else:
            df = pd.DataFrame(columns=header)
        if not self.strandness:
            df['strand'] = '*'
        return df


    def cal_frag_size(self, bam=None, chunk=1000000):
        """Extract the read length
        length count id
        """
        if bam is None:
            bam = self.bam
        # for SE or PE
        if self.is_paired(bam):
            pass
        else:
            pass
        # empty check
        if self.is_empty(bam):
            log.error('bam is empty: {}'.format(bam))
            return pd.DataFrame(columns=['length','strand','count','id']) #!!
        # read sam/bam file
        sam = pysam.AlignmentFile(bam)
        counter  = 0
        frag_size_list = []
        frames = []
        for read in sam:
            if self.as_se:
                # reads sizes
                if not read.is_unmapped \
                and not read.is_duplicate > 0:
                    counter += 1
                    strand = '-' if read.is_reverse else '+'
                    frag_size_list.append([read.infer_query_length(), strand, 1])
            else:
                strand = '*'
                # fragment sizes
                if read.is_proper_pair \
                and not read.is_unmapped \
                and not read.mate_is_unmapped \
                and not read.is_read1 \
                and not read.is_duplicate \
                and read.template_length > 0:
                    counter += 1
                    frag_size_list.append([read.template_length, strand, 1])
            # sample size
            if self.max_count > 0 and counter  > self.max_count:
                log.info('stop at: {}'.format(counter))
                break # stop
            # chunk
            if counter > 0 and counter%chunk == 0:
                frames.append(self.cal_freq(frag_size_list))
                frag_size_list = [] # empty
                log.info('{} : {} {}'.format('Processed', counter , self.labels))
        # last chunk
        if len(frag_size_list) > 0:
            frames.append(self.cal_freq(frag_size_list))
            frag_size_list = [] # empty
            log.info('{} : {} {}'.format('Processed', counter , self.labels))
        # overall
        df = pd.concat(frames, axis=0).groupby(['length', 'strand']).sum().reset_index()
        df['id'] = self.labels
        return df


    def distribution(self):
        """Basic statistics values
        value + freq
        mean, medium, mode, std, min, max, Q1, Q2, Q3
        """
        if self.freq_table.shape[0] == 0:
            out = pd.DataFrame(
                columns=['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1',
                         'Q2', 'Q3'])
        else:
            val = self.freq_table['length']
            freq = self.freq_table['count']
            inserts = np.repeat(val, freq)
            # statistics
            q_mean = np.mean(inserts)
            q_median = np.median(inserts)
            q_median_dev = np.median(np.absolute(inserts - q_median))
            q_mode = val[np.argmax(freq)]
            q_std = np.std(inserts)
            q_min = np.min(inserts)
            q_max = np.max(inserts)
            q_qual = np.quantile(inserts, [0.25, 0.5, 0.75], axis=0)
            # core distribution
            s = np.array([q_mean, q_median, q_mode, q_std, q_min, q_max]).round(2)
            s = np.append(s, q_qual)
            # DataFrame
            out = pd.DataFrame(s).T
            out.columns = ['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1',
                           'Q2', 'Q3']
        return out


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.csv',
            delete=False)
        return tmp.name


    def save_as(self, csv_file=None):
        """Save to file"""
        if csv_file is None:
            csv_file = self.csv_file # default
        if csv_file is None:
            csv_file = self._tmp()
        log.info('saving to file: {}'.format(csv_file))
        try:
            self.freq_table.to_csv(csv_file, index=False)
            # save statistics
            stat_file = csv_file + '.stat'
            da = self.distribution()
            da.to_csv(stat_file, sep='\t', index=False)
        except:
            log.warning('failed saving file: {}'.format(csv_file))


    def plot(self, plot_file=None):
        """Generate freq table plot
        line plot, hist
        """
        if plot_file is None:
            plot_file = self.plot_file
        hiseq_dir = os.path.dirname(hiseq.__file__)
        frag_plot_r = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
        stdout = os.path.join(self.outdir, self.labels+'.fragsize.plot.stdout')
        stderr = os.path.join(self.outdir, self.labels+'.fragsize.plot.stderr')
        cmd_file = os.path.join(self.outdir, self.labels+'.fragsize.plot.cmd.sh')
        cmd = ' '.join([
            'Rscript',
            frag_plot_r,
            plot_file,
            self.csv_file,
            '1> {}'.format(stdout),
            '2> {}'.format(stderr)
        ])
        with open(cmd_file, 'wt') as w:
            w.write(cmd+'\n')
        try:
            run_shell_cmd(cmd)
        except:
            log.warning('fragsize.py failed')
        return(plot_file)


    def run(self):
        self.freq_table = self.cal_frag_size(self.bam) # dataframe
        self.save_as() # to csv
        self.plot() # to pdf



def frag_size_picard(bam, outdir=None, smp_name=None):
    """Calculate PE fragment size, insertion size using Picard
    
    see: https://gatk.broadinstitute.org/hc/en-us/articles/360056968712-CollectInsertSizeMetrics-Picard-
    
    Example:

    java -jar picard.jar CollectInsertSizeMetrics \
      I=input.bam \
      O=insert_size_metrics.txt \
      H=insert_size_histogram.pdf \
      M=0.5
      
    Output: 
    ## htsjdk.samtools.metrics.StringHeader
    # CollectInsertSizeMetrics HISTOGRAM_FILE=insert_histogram.pdf MINIMUM_PCT=0.5 INPUT=rep1.bam OUTPUT=insert_metrics.txt    DEVIATIONS=10.0 METRIC_ACCUMULATION_LEVEL=[ALL_READS] INCLUDE_DUPLICATES=false ASSUME_SORTED=true STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
    ## htsjdk.samtools.metrics.StringHeader
    # Started on: Tue Jun 01 18:01:14 CST 2021

    ## METRICS CLASS        picard.analysis.InsertSizeMetrics
    MEDIAN_INSERT_SIZE      MODE_INSERT_SIZE        MEDIAN_ABSOLUTE_DEVIATION       MIN_INSERT_SIZE MAX_INSERT_SIZE MEAN_INSERT_SIZE        STANDARD_DEVIATION    READ_PAIRS       PAIR_ORIENTATION        WIDTH_OF_10_PERCENT     WIDTH_OF_20_PERCENT     WIDTH_OF_30_PERCENT     WIDTH_OF_40_PERCENT     WIDTH_OF_50_PERCENT   WIDTH_OF_60_PERCENT      WIDTH_OF_70_PERCENT     WIDTH_OF_80_PERCENT     WIDTH_OF_90_PERCENT     WIDTH_OF_95_PERCENT     WIDTH_OF_99_PERCENT     SAMPLE  LIBRARYREAD_GROUP
    177     64      79      22      2101    184.303587      104.531204      9495677 FR      31      63      95      127     159     187     213     239     315   385      1071

    ## HISTOGRAM    java.lang.Integer
    insert_size     All_Reads.fr_count
    22      7
    """
    # check: picard
    if isinstance(bam, str) and file_exists(bam):
        # outdir
        if outdir is None:
            outdir = os.path.dirname(bam)
        check_path(outdir)
        # smp_name
        if smp_name is None:
            smp_name = file_prefix(bam)
        # files
        out_metrics = os.path.join(outdir, smp_name+'.insert_meterics.txt')
        cmd_txt = os.path.join(outdir, smp_name+'.insert_mertics.cmd.sh')
        stdout = os.path.join(outdir, smp_name+'.insert_mertics.stdout')
        stderr = os.path.join(outdir, smp_name+'.insert_mertics.stderr')
        # command
        cmd = ' '.join([
            '{}'.format(which('picard')), # 
            'CollectInsertSizeMetrics',
            'I={}'.format(bam),
            'O={}'.format(out_metrics),
            'M=0.05',
            '1> {}'.format(stdout),
            '2> {}'.format(stderr),
        ])
        with open(cmd_txt, 'wt') as w:
            w.write(cmd+'\n')
        try:
            run_shell_cmd(cmd)
        except:
            log.error('frag_size_picard() failed')
    else:
        log.error('frag_size_picard() failed, expect str, got {}'.format(
            type(bam).__name__))
        

def frag_size_deeptools(bam, out):
    """Calculate PE fragment size, using deeptools
    Example:

    deepTools2.0/bin/bamPEFragmentSize \
    -hist fragmentSize.png \
    -T "Fragment size of PE RNA-seq data" \
    --maxFragmentLength 1000 \
    -b testFiles/RNAseq_sample1.bam testFiles/RNAseq_sample2.bam \
    testFiles/RNAseq_sample3.bam testFiles/RNAseq_sample4.bam \
    -samplesLabel sample1 sample2 sample3 sample4
    """
    pass


def frag_size_samtools(infile, outfile=None):
    """
    Calculate Fragment length from BAM file
    BAM: column-9

    BAM is PE
    BAM is sorted
    BAM is indexed
    """
    d = {}
    for read in pysam.AlignmentFile(infile):
        if not read.is_proper_pair:
            continue
        if not read.is_read1:
            continue
        size = abs(read.tlen)
        d[size] = d.get(size, 0) + 1
    d = sorted(d.items(), key=lambda x: x[0])

    if not outfile is None:
        with open(outfile, 'wt') as w:
            for k, v in d:
                w.write(str(k) + '\t' + str(v) + '\n')
    return d


def get_args():
    parser = argparse.ArgumentParser(description='hiseq fragsize -i bam -o outdir')
    parser.add_argument('-i', '--bam', nargs='+', required=True,
        help='BAM files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-l', '--labels', nargs='+', default=None,
        help='label of the bam files')
    parser.add_argument('--se', dest='as_se', action='store_true',
        help='Treat BAM file as SE, calculate the read size')
    parser.add_argument('-p', '--threads', default=4, type=int,
        help='number of processes, default: [4]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    return parser


def main():
    args = vars(get_args().parse_args())
    BamFragSize(**args).run()


if __name__ == '__main__':
    main()

#



# class BamPEFragSize(object):
#     """
#     Calculate insert size of PE, single BAM

#     sample size = 1000

#     > Table.csv
#     length strand count
#     """
#     def __init__(self, bam, labels=None, max_count=0):
#         self.bam = bam
#         self.max_count = max_count
#         self.labels = labels if labels else file_prefix(bam)
#         self.freq_table = self.cal_frag_size() # dataframe


#     def is_bam(self):
#         """Check input is BAM file
#         update: self.bam

#         string
#         *.bam
#         """
#         out = False
#         if isinstance(self.bam, str):
#             out = self.bam.endswith('.bam') and os.path.exists(self.bam)
#         if not out:
#             raise ValueError('bam file expected, got {}'.format(self.bam))
#         return out


#     def is_paired(self, topn=10000):
#         """Check input bam is Paired end alignment"""
#         out = False
#         if self.is_bam():
#             samfile = pysam.AlignmentFile(self.bam)
#             out = all([read.is_paired for read in samfile.head(topn)])
#             samfile.close()
#         return out


#     def is_empty(self, bam):
#         """Check input file is empty or not
#         pysam.Samfile().count
#         pysam.Samfile().mapped
#         """
#         pysam.index(bam)
#         sam = pysam.Samfile(bam)
#         return sam.count() == 0 and sam.mapped == 0


#     def cal_freq(self, x, return_dataframe=True):
#         """Calculate the frequency of list
#         return dataframe

#         index count
#         """
#         if isinstance(x, list):
#             var, freq = np.unique(x, return_counts=True)
#             df = pd.DataFrame(data=freq, index=var, columns=['count'])
#         else:
#             df = pd.DataFrame(columns=['count'])
#         return df


#     def cal_frag_size(self, bam=None, chunk=1000000):
#         """Extract the fragment length of paired end alignment
#         return
#         length count id
#         """
#         out = pd.DataFrame(columns = ['length','count','id']) # empty
#         if bam is None:
#             bam = self.bam
#         if not self.is_paired():
#             log.warning('not a PE bam, cal_frag_size() skipped ...')
#             return out
#         if self.is_empty(bam):
#             log.error('bam is empty: {}'.format(bam))
#             return out
#         # is_paired
#         sam = pysam.AlignmentFile(bam)
#         flag = 0
#         frag_size_list = []
#         df = pd.DataFrame(columns=['count'])
#         for read in sam:
#             if read.is_proper_pair \
#             and not read.is_unmapped \
#             and not read.mate_is_unmapped \
#             and not read.is_read1 \
#             and not read.is_duplicate \
#             and read.template_length > 0:
#                 flag += 1
#                 frag_size_list.append(read.template_length)
#             # get sample size
#             if self.max_count > 0 and flag > self.max_count:
#                 log.info('stop iteration at: {}'.format(flag))
#                 break # stop
#             # chunk
#             if flag > 0 and flag % chunk == 0:
#                 dfn = self.cal_freq(frag_size_list)
#                 df = pd.concat([df, dfn], axis=1).sum(axis=1)
#                 frag_size_list = [] # empty
#                 log.info('{} : {} {}'.format('Processed', flag, self.labels))
#         # last chunk
#         if len(frag_size_list) > 0:
#             dfn = self.cal_freq(frag_size_list)
#             df = pd.concat([df, dfn], axis=1).sum(axis=1)
#             frag_size_list = [] # empty
#             log.info('{} : {} {}'.format('Processed', flag, self.labels))
#         # convert to data.frame
#         df = df.reset_index()
#         df.columns = ['length', 'count']
#         df['id'] = self.labels
#         return df


#     def distribution(self):
#         """Basic statistics values
#         value + freq

#         mean, medium, mode, std, min, max, Q1, Q2, Q3
#         """
#         if self.freq_table.shape[0] == 0:
#             out = pd.DataFrame(
#                 columns=['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1',
#                          'Q2', 'Q3'])
#         else:
#             val = self.freq_table['length']
#             freq = self.freq_table['count']
#             inserts = np.repeat(val, freq)
#             # inserts = np.repeat(self.freq_table['length'], self.freq_table['count'])
#             # statistics
#             q_mean = np.mean(inserts)
#             q_median = np.median(inserts)
#             q_median_dev = np.median(np.absolute(inserts - q_median))
#             q_mode = val[np.argmax(freq)]
#             q_std = np.std(inserts)
#             q_min = np.min(inserts)
#             q_max = np.max(inserts)
#             q_qual = np.quantile(inserts, [0.25, 0.5, 0.75], axis=0)
#             # core distribution
#             s = np.array([q_mean, q_median, q_mode, q_std, q_min, q_max])
#             s = np.append(s, q_qual)
#             # DataFrame
#             out = pd.DataFrame(s).T
#             out.columns = ['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1',
#                            'Q2', 'Q3']
#         return out


#     def _tmp(self):
#         """
#         Create a tmp file to save json object
#         """
#         tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.csv',
#             delete=False)
#         return tmp.name


#     def save_as(self, csv_file=None):
#         """Save to file"""
#         if csv_file is None:
#             csv_file = self._tmp()
#         log.info('saving to file: {}'.format(csv_file))
#         try:
#             # save table
#             self.freq_table.to_csv(csv_file, index=False)
#             # save statistics
#             stat_file = csv_file + '.stat'
#             da = self.distribution()
#             da.to_csv(stat_file, sep='\t', index=False)
#         except:
#             log.warning('failed saving file: {}'.format(csv_file))


#     def plot(self, plot_file):
#         """Generate freq table plot
#         line plot, hist
#         """
#         # save to csv file
#         csv_file = os.path.splitext(plot_file)[0] + '.csv'
#         self.save_as(csv_file) # save csv, stat
#         # save to plot
#         # to-do: matplotlib function
#         hiseq_dir = os.path.dirname(hiseq.__file__)
#         frag_plot_r = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
#         stdout = os.path.join(self.outdir, 'frag_plot.plot.stdout')
#         stderr = os.path.join(self.outdir, 'frag_plot.plot.stderr')
#         cmd_file = os.path.join(self.outdir, 'frag_plot.plot.cmd.sh')
#         cmd = ' '.join([
#             'Rscript',
#             frag_plot_r,
#             plot_file,
#             csv_file,
#             '1> {}'.format(stdout),
#             '2> {}'.format(stderr)
#         ])
#         with open(cmd_file, 'wt') as w:
#             w.write(cmd+'\n')
#         try:
#             os.system(cmd)
#         except:
#             log.warning('fragsize.py failed')


#     def run(self):
#         """Generate dataframe, pdf
#         bam, str or list
#         """
#         if isinstance(self.bam, str):
#             bam_list = [self.bam]
#         elif isinstance(self.bam, list):
#             bam_list = self.bam
#         else:
#             raise Exception('str and list expected, {} got'.format(type(self.bam).__name__))
#         # create df
#         frames = [self.cal_frag_size(b) for b in bam_list]
#         df = pd.concat(frames, axis=1)
#         # save to txt
#         df.to_csv(frag_stat, index=False)
#         # save to plot
#         hiseq_dir = os.path.dirname(hiseq.__file__)
#         frag_plot_r = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
#         stdout = os.path.join(self.outdir, 'frag_plot.plot.stdout')
#         stderr = os.path.join(self.outdir, 'frag_plot.plot.stderr')
#         cmd_file = os.path.join(self.outdir, 'frag_plot.plot.cmd.sh')
#         cmd = ' '.join([
#             'Rscript',
#             frag_plot_r,
#             plot_file,
#             csv_file,
#             '1> {}'.format(stdout),
#             '2> {}'.format(stderr)
#         ])
#         with open(cmd_file, 'wt') as w:
#             w.write(cmd+'\n')
#         try:
#             os.system(cmd)
#         except:
#             log.warning('fragsize.py failed')