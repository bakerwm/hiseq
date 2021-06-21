#!/usr/bin/env python

"""
Mission: guess the hiseq library

1. TruSeq 
  - i7 
  - barcode (----A{bc}-adapter)
  - NSR (read1: CT----TC-adapter
         read2: GA----AG-adapter)

2. Nextera
  - i7

Statistics:
1. Percentage of P7 adapter (p7a)  
2. Percentage of i7 index 
3. Percentage of barcode

output:
is_truseq
is_nextera
is_nsr
has_barcode

Parse the i7 and inline barcode from sequencing reads

support: TruSeq library (122nt+INS)
# p5 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' # 58nt
# p7a = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' # 34nt 
# p7b = 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
# p7 # 6nt

support: Nextera library (128nt+INS)
# p5 = 'AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG' # 62nt
# p7a = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' # 34nt
# p7b = 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
# p7 # 8nt


output format:
1. P7 structure: P7a-i7-P7b
"group", "total", "count", "pct", "sample"
#P7a_i7 17107874        16313570        95.4    RNAseq_shWhite_0-2h_rep1_1

2. i7 and bc pct: 
seq  i7_rc   i7_name   demx  demx_pct group total i7 i7_pct sample
>CAACTA CAACTA  TruSeq_Index29  16210525         99.4   i7      17107874        16313570        95.4    RNAseq_shWhite_0-2h_rep1_1

## ChangeLog

## Version-3 - 2021-06-15
1. guess library, 
2. support NSR

## Version-2 - 2021-05-11
1. support Nextera 
2. support large insert, (partial P7)

## Version-1 - 2021-05-11 
1. check full version of p7a, p7b
2. search barcode and i7
"""

import os
import sys
import re
import argparse
import pathlib
import shutil
import pyfastx
import hiseq
from multiprocessing import Pool
from hiseq.utils.fastx import Fastx
from hiseq.utils.utils import update_obj, log, run_shell_cmd, Config
from hiseq.utils.file import fx_name, check_path, list_file, list_fx, \
    file_abspath, file_exists, file_abspath
from hiseq.demx.sample_sheet import HiSeqIndex # check index name/seq
from collections import Counter


class HiseqLib(object):
    """"
    check HiSeq library for multiple fq 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args() # update


    def init_args(self):
        args_init = {
            'fq': None,
            'outdir': None,
            'parallel_jobs': 1,
            'overwrite': False,
            'n_max': 100000,
            }
        self = update_obj(self, args_init, force=False)
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.fq_list = self.init_fq(self.fq)


    def init_fq(self, fq):
        """
        Make sure fq is read1 of PE reads
        """
        out = []
        if isinstance(fq, str):
            fq = os.path.abspath(fq)
            if os.path.isdir(fq):
                out = list_fx(fq)
                out = [i for i in out if '_1.' in i] # rep1
            elif os.path.isfile(fq):
                fname = fx_name(fq, fix_pe=False)
                if fname.endswith('_1') and os.path.exists(fq):
                    out = [fq]
            else:
                log.warning('illegal fq, str expect, got {}'.format(
                    type(fq).__name__))
        elif isinstance(fq, list):
            out = [i for k in fq for i in self.init_fq(k)]
        else:
            log.warning('illegal fq, str or list expect, got {}'.format(
                    type(fq).__name__))
        return out


    def run_single(self, x):
        p = HiseqLibR1(x, outdir=self.outdir, n_max=self.n_max)
        if file_exists(p.p7_json) and not self.overwrite:
            log.info('parse_i7 skipped, file exists: {}'.format(p.p7_json))
        else:
            p.parse_i7()


    def run_multi(self, x):
        if self.parallel_jobs > 1 and len(x) > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single, x)
        else:
            for i in x:
                self.run_single(i)

    
    def report(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_p7_report.R')
        report_html = os.path.join(self.outdir, 'HiSeq_P7_report.html')
        stdout = os.path.join(self.outdir, 'report.stdout')
        stderr = os.path.join(self.outdir, 'report.stderr')
        cmd_file = os.path.join(self.outdir, 'cmd.sh')
        cmd = ' '.join([
            shutil.which('Rscript'),
            qc_reportR,
            self.outdir,
            '1> {}'.format(stdout),
            '2> {}'.format(stderr)
        ])
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(report_html) and self.overwrite is False:
            log.info('file exists, skip generating html.')
        else:
            run_shell_cmd(cmd)
        if not os.path.exists(report_html):
            log.error('failed, generating html file')


    def run(self):
        check_path(self.outdir)
        self.run_multi(self.fq_list)
        self.report()

        
class HiseqLibR1(object):
    """    
    guess the fastq file is NSR or not: single fq
    
    1. is nextera
    2. read1: CT...; read2: GA...    
    """
    def __init__(self, fq, outdir=None, n_max=1000000):
        self.fq = fq
        self.outdir = outdir
        self.n_max = n_max
        if not isinstance(fq, str):
            raise ValueError('str expect, got {}'.format(
                type(fq).__name__))
        if not file_exists(fq):
            raise ValueError('file not exists: {}'.format(fq))
        self.init_args()
        self.lib = guess_adapter(fq, show_log=False)
        if self.lib is None:
            log.warning('no adapters found for: {}'.format(fq))
            self.lib = 'truseq'
        self.hiseq_type = self.lib
        self.is_nextera = self.lib == 'nextera'
        self.is_truseq = self.lib == 'truseq'
        self.is_smallrna = self.lib == 'smallrna'
        self.is_nsr = self.guess_nsr()
        Config().dump(self.__dict__, self.config_yaml)        


    def init_args(self):
        self.fq = file_abspath(self.fq)
        self.fname = fx_name(self.fq, fix_pe=False)
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.config_yaml = os.path.join(self.outdir, self.fname+'.config.yaml')
        self.p7_json = os.path.join(self.outdir, self.fname+'.p7.json')
        self.i7_json = os.path.join(self.outdir, self.fname+'.i7.json')
        self.barcode_json = os.path.join(self.outdir, self.fname+'.barcode.json')
        check_path(self.outdir)        


    def guess_nsr(self, cutoff=0.9):
        """
        cutoff : float
            Percentage of reads contains the prefix, default: [0.9]
        read1: CT....
        read2: GA....
        """
        n_max = 100000 # 100k reads
        i = 0
        n_ct = 0
        n_ga = 0
        for _,s,_,_ in pyfastx.Fastx(self.fq):
            i += 1 # counter
            n_ct += int(s.startswith('CT'))
            n_ga += int(s.startswith('GA'))
            if i >= n_max:
                break
        # check percentage
        is_read1 = self.fname.endswith('1')
        is_read2 = self.fname.endswith('2')
        pct_ct = n_ct/i
        pct_ga = n_ga/i
        # check filename
        if is_read1:
            out = pct_ct > cutoff
        elif is_read2:
            out = pct_ga > cutoff
        else:
            r = 'rea1' if pct_ct > cutoff else \
                'read2' if pct_ga > cutoff else None
            if isinstance(r, str):
                log.error('possible a NSR {}, but filename not matched: {}'.format(
                    r, self.fname))
            out = False
        # output
        return out

    
    def parse_i7(self):
        """
        Guess the i7 index
        located between p7a and p7b
        
        # i7.stat
        index_1, index_2, ...
        
        # p7.stat
        p7a, p7a-p7b, p7b        
        """
        # get the p7a, p7b sequence
        a = HiseqAdapter(self.lib) # p7a, p7b, p5a, p5b, ...
        # index length: defaut: nextera:8, truseq:6
        isizes = {
            'truseq': 6,
            'nextera': 8
        }
        index_size = isizes.get(self.lib, 6) # default
        # patterns
        p1 = re.compile('([ACGTN]{6})'+a.p7a[:6]) # p7a; barcode
        p2 = re.compile(a.p7a[-6:]+'([ACGTN]{5,10})'+a.p7b[:6]) # p7a_i7_p7b
        p3 = re.compile(a.p7b[:6]) # p7b;
        # parse file
        i = 0
        n1 = 0
        n2 = 0
        n3 = 0
        bc = []
        i7 = []
        for _,s,_,_ in pyfastx.Fastx(self.fq):
            i += 1
            if i%1e6 == 0:
                print('{} processed'.format(i))
            # for bc+p7a
            s1 = p1.search(s)
            if s1:
                n1 += 1
                bc_seq = s1.group(1)
                bc.append(bc_seq)
                # for i7
                s2 = p2.search(s)
                s3 = p3.search(s)
                if s2:
                    n2 += 1
                    i7_seq = s2.group(1)
                    i7.append(i7_seq)
                elif s3:
                    n3 += 1
                    i7_seq = '-'
                else:
                    i7_seq = '-'
            else:
                bc_seq = '-'
                i7_seq = '-'
            # skipped
            if i >= self.n_max and self.n_max > 0: # topN seq
                break
        # wrap log
        # for p7; p7a, p7a-p7b, p7b
        p7_d = {
            'p7a': n1-n2-n3,
            'p7a_p7b': n2,
            'p7b': n3,
            'no_ad': i-n1
        }
        Config().dump(p7_d, self.p7_json)
        # show msg
        p7_fmt = '#{}\t{}\t{}\t{:.1f}\t{}'
        msg = '\n'.join([
            p7_fmt.format('p7a', i, p7_d['p7a'], p7_d['p7a']/i*100, self.fname),
            p7_fmt.format('p7a_p7b', i, n2, n2/i*100, self.fname),
            p7_fmt.format('p7b', i, n3, n3/i*100, self.fname),
            p7_fmt.format('no_ad', i, i-n1, (i-n1)/i*100, self.fname)
        ])
        print(msg)
        # for i7: freq
        top_i7 = self.top_seq(i7, n=5)
        # print('!AAAA-1', 'i7', len(i7))
        Config().dump(top_i7, self.i7_json)
        # for barcode: freq
        top_bc = self.top_seq(bc, n=5, revcmp=True)
        Config().dump(top_bc, self.barcode_json)


    def top_seq(self, x, n=5, revcmp=False):
        """
        Top N seq in list
        """
        out = {}
        if isinstance(x, list):
            t = len(x)
            a = Counter(x).most_common(n)
            for i in a:
                k, v = i
                rp = self.rev_comp(k) if revcmp else k
                rn = HiSeqIndex(rp).name # NULL if not found
                if rn is None or rn == 'NULL': # in case, toml, 'NULL'
                    rn = rp
                # seq, seq-rev, seq-name, pct, count
                out.update({
                    k: {
                    'seq': k,
                    'seq_revcomp': rp,
                    'name': rn,
                    'count': v,
                    'pct': v/t*100,
                }})
                # out.append('{}\t{}\t{}\t{}\t{:5.1f}'.format(k, rp, rn, v, v/t*100))
            # sorted
            out = dict(sorted(out.items(), key=lambda x: x[1]['count'], reverse=True))
        return out
            

    def rev_comp(self, x):
        """Reverse complement DNAseq"""
        t = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'N': 'N',
        }
        s = x[::-1] # reverse
        return ''.join([t.get(i, i) for i in s]) # complement


class HiseqAdapter(object):
    """
    support: TruSeq library (122nt+INS)
    # p5 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' # 58nt
    # p7a = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' # 34nt 
    # p7b = 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
    # p7 # 6nt

    support: Nextera library (128nt+INS)
    # p5 = 'AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG' # 62nt
    # p7a = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' # 34nt
    # p7b = 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
    # p7 # 8nt
    
    Sturcture
    P5a{i5}P5b--insert--P7a{i7}P7b
    
    TruSeq: P5a{58nt} + P7a{i7}P7b{64nt} = 122nt
    Nextera: P5a{62nt} + p7a{i7}P7b{66nt} = 128nt
    
    Current: single-index library (i7)
    """
    def __init__(self, lib='TruSeq'):
        lib = lib.lower()
        if not lib in ['truseq', 'nextera']:
            raise ValueError('unknown lib, expect [truseq, nextera], got: {}'.format(lib))
        self.hiseq_type = lib
        self.p5a = self.get_p5a()
        self.p5b = self.get_p5b(lib)
        self.p7a = self.get_p7a(lib)
        self.p7b = self.get_p7b()
        
    
    def get_p7(self, lib='truseq'):
        d = {
            'truseq': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG', # 64nt
            'nextera': 'CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG', # 66nt
        }
        return d.get(d, '')
    
    
    def get_p7a(self, lib='truseq'):
        d = {
            'truseq': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', # 34nt
            'nextera': 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC', # 34nt
        }
        return d.get(lib, '')
        

    def get_p7b(self):
        return 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
        

    def get_p5(self, lib='truseq'):
        d = {
            'truseq': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
            'nextera': 'AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
        }
        return d.get(lib, '')
        

    def get_p5a(self):
        return 'AATGATACGGCGACCACCGAGATCTACAC', # 29nt
        

    def get_p5b(self, lib='truseq'):
        d = {
            'truseq': 'TCTTTCCCTACACGACGCTCTTCCGATCT', # 29nt
            'nextera': 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG' # 33nt
        }
        return d.get(lib, '')


def guess_adapter(fq, show_log=True):
    """
    see: hiseq.utils.fastx.Fastx().detect_adapter()
    Guess the hiseq library type, by the p7a adapter:
    'truseq': 'AGATCGGAAGAGC',
    'nextera': 'CTGTCTCTTATA',
    'smallrna': 'TGGAATTCTCGG'
    """
    
    d = Fastx(fq).detect_adapter(show_log) # dict: {seqtype:count}
    if len(d) > 1:
        dx = d.copy()
        dx.pop('unknown', 0) # remove unknown
        n_total = sum(d.values())
        n_max = max(dx.values())
        i = []
        for k,v in dx.items():
            if v == n_max:
                i.append(k)
        # check unique
        if len(set(i)) > 1:
            log.error('multiple hiseq libraries detected: {}'.format(
                ','.join(i)))
        # log
        if show_log:
            msg = ','.join([
                '{}:{:.2f}%'.format(j, d.get(j, 0)/n_total*100) for j in i
            ])
            print(msg)        
        # output
        out = i.pop()
    else:
        log.info('unknown library')
        out = None
    return out


def get_args():
    parser = argparse.ArgumentParser(
        description='hiseq i7')
    parser.add_argument('-i', '--fq', nargs='+', required=True,
        help='FASTQ files, or directory contains fastq files')
    parser.add_argument('-o', '--outdir', default=None, required=True,
        help='The directory to save results.')
    parser.add_argument('-n', '--n-max', dest='n_max', type=int, default=100000,
        help='Maximum number of reads to process, default: [100000]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='overwrite the exists files')
    return parser
        

def main():
    args = vars(get_args().parse_args())
    HiseqLib(**args).run()

    
if __name__ == '__main__':
    main()
        
# 