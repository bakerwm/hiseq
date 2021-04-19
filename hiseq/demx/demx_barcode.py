#!/usr/bin/env python3


"""
Function:
1. Demultiplex barcode

## others
1. defq, ultra fasta multi-threaded fq demultiplexing
https://github.com/OpenGene/defq

2. deML, maxlikelihood demultiplexing
https://github.com/grenaud/deML

## test 
step1: split by p7-index: index1_index2.r1.fq
step2: split by inline barcode:

## search index
index1:index2:barcode
"""


import os
import sys
import re
import pathlib
import argparse
from xopen import xopen
from contextlib import ExitStack
import Levenshtein as lev # distance
from hiseq.utils.helper import update_obj, check_file, check_path, file_abspath, combinations, log, Config
from hiseq.utils.seq import Fastx


class DemxBarcode(object):
    """Demultiplex barcode, locate in the head of read1/2
    
    Example:
    
    >>> args = {
        'fq1': '1m.r1.fq.gz',
        'fq2': '1m.r2.fq.gz',
        'outdir': 'aaa',
        'index_table': 'index.csv',
        'mismatch': 0,
        'demo': False,
        'gzipped': True,
    }

    >>> DemxBarcode(**args).run()

    --------------------------------------------------------------------------------
    Demultiplex:
                     fq1 | /data/1m.r1.fq.gz
                     fq2 | /data/1m.r2.fq.gz
                 PE_mode | yes     
             bc_in_read2 | yes     
               barcode_n_left | 0       
              barcode_n_right | 1       
                bc_width | 6       
             index_table | /data/index.csv
                  outdir | /data/aaa
                mismatch | 0       
               overwrite | no      
                run_demo | no      
    --------------------------------------------------------------------------------

    --------------------------------------------------------------------------------
    Demultiplex report
    Input reads :    1000000    1.0M
    order filename                                      count percent million  to_400M 
      1 CLIPseq_PTB_UV254_RIP_rep1                       89   0.01%      0.0      0.0
      2 CLIPseq_PTB_UV254_input_rep1                     48   0.00%      0.0      0.0
      3 CLIPseq_PTB_noUV_input_rep1                       2   0.00%      0.0      0.0
      4 CLIPseq_YFP_UV254_RIP_rep1                       35   0.00%      0.0      0.0
      5 CLIPseq_YFP_UV254_input_rep1                     48   0.00%      0.0      0.0
      6 RNAseq_I_R_Hybrid_Dys_6h_rep1                   178   0.02%      0.0      0.1
      7 RNAseq_I_R_Hybrid_Dys_6h_rep2                    32   0.00%      0.0      0.0
      8 RNAseq_I_R_Hybrid_None_Dys_6h_rep1               78   0.01%      0.0      0.0
      9 RNAseq_I_R_Hybrid_None_Dys_6h_rep2               47   0.00%      0.0      0.0
     10 undemx                                      999,443  99.94%      1.0    399.8
    --------------------------------------------------------------------------------
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
    
    
    def init_args(self):
        args_init = {
            'fq1': None,
            'fq2': None,
            'in_read2': True,
            'outdir': None,
            'barcode_n_left': 0,
            'barcode_n_right': 1,
            'index_table': None, #name,index
            'mismatch': 0,
            'overwrite': False,
            'demo': False,
            'gzipped': False,
        }
        self = update_obj(self, args_init, force=False)        
        # fastq file
        if not check_file(self.fq1, emptycheck=True):
            log.error('fq1, file not exists')
            raise ValueError('fq1, fastq file required')
        self.is_pe = isinstance(self.fq2, str) # os.path.exists(self.fq2)
        if not check_file(self.index_table, emptycheck=True):
            raise ValueError('index_table not exists or empty')
        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        if not self.mismatch in range(4):
            raise ValueError('illegal mimatche: [{}], expect [0,1,2,3]'.format(
                self.mismatch))
        # index table
        self.idx = self.load_index()
        # compatiable
        if not self.check_index():
            raise ValueError('index not compatiable with mimatch={}'.format(
            self.mismatch))
        # absolute path
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.index_table = file_abspath(self.index_table)
        self.outdir = file_abspath(self.outdir)

    
    def check_index(self):
        """Make sure the index compatiable
        support the mismatches
        """
        idx_list = list(self.idx.keys())
        print('!AAAA-2', idx_list)
        d = []
        for i, j in combinations(idx_list, 2):
            d.append(not self.is_valid_index(i, j, self.mismatch))
        out1 = all(d)
        # check barcode width
        idx_width = [len(i) for i in idx_list]
        out2 = len(set(idx_width)) == 1 # identical width
        self.bc_width = idx_width.pop()
        return all([out1, out2])
    

    def sample_fq(self, smp='demo', sample_size=1000000):
        """
        Create a subsample of input fastq files, default: 1M reads
        Run the whole process for demostration        
        update: fq1, fq2, outdir
        """
        log.info('Running demo with subsample: {} reads'.format(sample_size))
        # update args
        self.outdir = os.path.join(self.outdir, smp)
        self.data_dir = os.path.join(self.outdir, 'data')
        check_path(self.data_dir)
        # subset
        demo_fq1 = os.path.join(self.data_dir, os.path.basename(self.fq1))
        Fastx(self.fq1).sample(demo_fq1, sample_size)
        self.fq1 = demo_fq1
        if self.fq2:
            demo_fq2 = os.path.join(self.data_dir, os.path.basename(self.fq2))
            Fastx(self.fq2).sample(demo_fq2, sample_size)
            self.fq2 = demo_fq2
    
    
    def show_msg(self):
        msg = '\n'.join([
            '-'*80,
            'Demultiplex:',
            '{:>20} | {:<8}'.format('fq1', self.fq1),
            '{:>20} | {:<8}'.format('fq2', self.fq2 if self.fq2 else 'None'),
            '{:>20} | {:<8}'.format('PE_mode', 'yes' if self.is_pe else 'no'),
            '{:>20} | {:<8}'.format('bc_in_read2', 'yes' if self.in_read2 else 'no'),
            '{:>20} | {:<8}'.format('barcode_n_left', self.barcode_n_left),
            '{:>20} | {:<8}'.format('barcode_n_right', self.barcode_n_right),
            '{:>20} | {:<8}'.format('bc_width', self.bc_width),
            '{:>20} | {:<8}'.format('index_table', self.index_table),
            '{:>20} | {:<8}'.format('outdir', self.outdir),
            '{:>20} | {:<8}'.format('mismatch', self.mismatch),
            '{:>20} | {:<8}'.format('overwrite', 'yes' if self.overwrite else 'no'),
            '{:>20} | {:<8}'.format('run_demo', 'yes' if self.demo else 'no'),
            '-'*80,
        ])
        print(msg)
    
    
    def load_index(self):
        d = {}
        with open(self.index_table) as r:
            for line in r:
                if line.startswith('#'):
                    continue
                s = re.split('[,\s\t]', line.strip())
                if not len(s) == 2:
                    log.error('illegal format, expect:name,index')
                    raise ValueError('illegal index_table: {}'.format(line))
                name,bc = s
                if not re.match('^[ACGTN$]+$', bc.strip()):
                    continue
                if bc in d:
                    raise ValueError('index not unique, {}'.format(
                        self.index_table))
                d[bc] = name
        return d

    
    def str_distance(self, x, y, partial=True):
        """
        Check distance between a, b
        """
        try:
            if partial:
                x = x[:len(y)]
                y = y[:len(x)]
            out = lev.distance(x, y)
        except:
            out = 10 # huge number
        return out
        
    
    def is_valid_index(self, x, y, mm=0):
        """
        Check if the index x and y, compatiable
        """
        return self.str_distance(x, y) <= mm
    

    def match_index(self, x, y, mm=0):
        """
        Search x(str) in y(dict)
        """
        for k,v in y.items():
            if self.is_valid_index(x, k, mm):
                out = v
                # To speed up, exit the loop when meeting the first hit.
                # if mm>0, this strategy might skip the best hit, when
                # mismatches hit at first.
                break # 
            else:
                out = None
        return out


    def readfq(self, fh): # this is a generator function
        """
        source: https://github.com/lh3/readfq/blob/master/readfq.py
        processing fastq file
        """
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fh: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            [name, _, comment], seqs, last = last[1:].partition(" "), [], None
            for l in fh: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None, comment # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fh: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs), comment; # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None, comment # yield a fasta record instead
                    break


    def index_se(self, fh, outdir):
        """
        barcode saved in read1/read2
        extract: {left}+{right} to comment
        """
        # prepare output files, named by index
        fnames = list(self.idx.values())
        fnames.append('undemx')
        f_files = [os.path.join(outdir, i + '.fq') for i in fnames]
        if self.gzipped:
            f_files = [i+'.gz' for i in f_files]
        # check status
        self.fn_toml = os.path.join(outdir, 'read_count.toml')
        if os.path.exists(self.fn_toml):
            log.info('Skipped {}, read_count.json exists'.format(outdir))
            return f_files
        # open multiple files
        bl = self.barcode_n_left
        br = self.barcode_n_right
        bw = self.bc_width
        fn = {}
        with ExitStack() as stack:
            fws = [stack.enter_context(xopen(f, 'wt')) for f in f_files]
            n = 0 # counter
            for name, seq, qual, comment in self.readfq(fh):
                n += 1
                if n%1e6 == 0:
                    log.info('Processed reads: {}'.format(n))
                fq = '\n'.join(['@' + name + ' ' + comment, seq, '+', qual])
                i = seq[bl:(bl+bw)] # left:(left+width)
                hit = self.match_index(i, self.idx, self.mismatch)
                if hit is None:
                    hit = 'undemx'
                fw = fws[fnames.index(hit)]
                fw.write(fq+'\n')
                fn[hit] = fn.get(hit, 0) + 1
        # save dict to json
        Config().dump(fn, self.fn_toml)
        return f_files
    
    
    def index_pe(self, fh1, fh2, outdir):
        """
        barcode saved in read1/read2
        """
        # prepare output files, named by index
        fnames = list(self.idx.values())
        fnames.append('undemx')
        f1_files = [os.path.join(outdir, i + '_1.fq') for i in fnames]
        f2_files = [os.path.join(outdir, i + '_2.fq') for i in fnames]
        if self.gzipped:
            f1_files = [i+'.gz' for i in f1_files]
            f2_files = [i+'.gz' for i in f2_files]
        # check status
        self.fn_toml = os.path.join(outdir, 'read_count.toml')
        if os.path.exists(self.fn_toml):
            log.info('Skipped {}, read_count.json exists'.format(outdir))
            return (f1_files, f2_files)
        # open multiple files
        bl = self.barcode_n_left
        br = self.barcode_n_right
        bw = self.bc_width
        fn = {}
        with ExitStack() as stack:
            fws1 = [stack.enter_context(xopen(f, 'wt')) for f in f1_files]
            fws2 = [stack.enter_context(xopen(f, 'wt')) for f in f2_files]
            n = 0 # counter            
            for r1, r2 in zip(self.readfq(fh1), self.readfq(fh2)):
                n += 1
                if n%1e6 == 0:
                    log.info('Processed reads: {}'.format(n))
                r1_name, r1_seq, r1_qual, r1_comment = r1
                r2_name, r2_seq, r2_qual, r2_comment = r2
                fq1 = '\n'.join(['@' + r1_name + ' ' + r1_comment, r1_seq, '+', r1_qual])
                fq2 = '\n'.join(['@' + r2_name + ' ' + r2_comment, r2_seq, '+', r2_qual])
                # check output
                seq = r2_seq if self.in_read2 else r1_seq
                i = seq[bl:(bl+bw)] # left:(left+width)
                hit = self.match_index(i, self.idx, self.mismatch)
                if hit is None:
                    hit = 'undemx'                
                fw1 = fws1[fnames.index(hit)]
                fw2 = fws2[fnames.index(hit)]
                fw1.write(fq1 + '\n')
                fw2.write(fq2 + '\n')
                fn[hit] = fn.get(hit, 0) + 1
        # save dict to json
        Config().dump(fn, self.fn_toml)
        return (f1_files, f2_files)

    
    def report(self):
        """Organize the report
        name, count, pct, to_400M
        """
        df = Config().load(self.fn_toml)
        total = sum(df.values())
        scale = 4e8/total # to_400M
        i = 0
        f_stat = []
        for k,v in df.items():
            i += 1
            s = '{:>3} {:<40s} {:>10,} {:6.2f}% {:8.1f} {:8.1f}'.format(
                i, k, v, v/total*100, v/1e6, v/1e6*scale)
            f_stat.append(s)
        # output
        msg = '\n'.join([
            '-'*80,
            'Demultiplex report',
            '{} : {:>10} {:6.1f}M'.format('Input reads', total, total/1e6),
            '{:>3} {:<40s} {:>10} {:6} {:8s} {:8s}'.format(
                'order', 'filename', 'count', 'percent', 'million', 'to_400M'),
            '\n'.join(f_stat),
            '-'*80,
        ])
        # save to file
        self.report_txt = os.path.join(self.outdir, 'report.txt')
        with open(self.report_txt, 'wt') as w:
            w.write(msg+'\n')
        print(msg)
    

    def run(self):
        check_path(self.outdir)
        self.show_msg()
        if self.demo:
            self.sample_fq()
        if self.is_pe:
            with xopen(self.fq1) as r1, xopen(self.fq2) as r2:
                self.index_pe(r1, r2, self.outdir)
        else:
            with xopen(self.fq1) as r1:
                self.index_se(r1, self.outdir)
        # report
        self.report()
