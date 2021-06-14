#!/usr/bin/env python

"""
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
from hiseq.utils.utils import update_obj, log, run_shell_cmd
from hiseq.utils.file import fx_name, check_path, list_file, list_fx
from hiseq.demx.sample_sheet import HiSeqIndex # check index name/seq
from collections import Counter



def top_seq(x, n=5, revcmp=False):
    """
    Top N seq in list
    """
    t = len(x)
    a = Counter(x).most_common(n)
    out = []
    for i in a:
        k, v = i
        # format:
        rp = rev_comp(k) if revcmp else k
        rn = HiSeqIndex(rp).name # NULL if not found
        if rn is None or rn == 'NULL':
            rn = rp
        # seq, seq-rev, seq-name, pct, count
        out.append('{}\t{}\t{}\t{}\t{:5.1f}'.format(k, rp, rn, v, v/t*100))
    return out


def rev_comp(x):
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


## guess adapter
def guess_adapter(fx):
    """
    Guess the adapters:
    TruSeq    AGATCGGAAGAGC
    Nextera   CTGTCTCTTATACACATCT
    smallRNA  TGGAATTCTCGG
    """
    ad = {
        'truseq': 'AGATCGGAAGAGC',
        'nextera': 'CTGTCTCTTATA',
        'smallrna': 'TGGAATTCTCGG'
        }
    d = {} # key:count
    n_max = 10000
    n = 0
    for _,s,_,_ in pyfastx.Fastx(fx):
        n += 1
        if n > n_max:
            break
        # count reads
        for k,v in ad.items():
            if v in s:
                d[k] = d.get(k, 0) + 1
    # sort d by values
    d = dict(sorted(d.items(), key=lambda x: x[1], reverse=True))
    out = None
    if len(d) > 0:
        if max(d.values()) == 0:
            log.warning('no adapters detected')
        else:
            m = []
            for k,v in d.items():
                m.append('{}:{:.1f}%'.format(k, v/n_max*100))
                if v == max(d.values()):
                    out = k # incase, top, equal, ...
            print('{}: {} [{}]  {}'.format('Guess_adapter', out, '; '.join(m), fx))
    else:
        log.warning('no adapters detected')
    # output
    return out
        
    
def get_p7_seq(x):
    """
    Return the sequence of P7: TruSeq, Nextera, smallRNA
    """
    p7 = {
        'truseq': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 6, 'ATCTCGTATGCCGTCTTCTGCTTG'],
        'nextera': ['CTGTCTCTTATACACATCTCCGAGCCCACGAGAC', 8, 'ATCTCGTATGCCGTCTTCTGCTTG'],
        'smallrna': [],
    }
    return p7.get(x.lower(), None)
    

## patterns
def parse_i7(fx, outdir, save_seq=False):
    fname = fx_name(fx, fix_pe=False)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    out_bc = os.path.join(outdir, fname + '.barcode.txt')
    out_msg = os.path.join(outdir, fname + '.stat.log')
    ## adapters
    lib_type = guess_adapter(fx)
    if lib_type is None:
        print('unknown lib type: {}'.format(fx))
        sys.exit(1)
    p7a, i7_size, p7b = get_p7_seq(lib_type) 
    ## patterns
    p1 = re.compile('([ACGTN]{6})'+p7a[:6]) # p7a; gather barcode
    p2 = re.compile(p7a[-6:]+'([ACGTN]{5,10})'+p7b[:6]) # p7a_i7_p7b
    p3 = re.compile(p7b[:6]) # p7b;
    ## parse file
    n0 = 0
    n1 = 0
    n2 = 0
    n3 = 0
    bc = []
    i7 = []
    with open(out_bc, 'wt') as w:
        for _,s,_,_ in pyfastx.Fastx(fx):
            n0 += 1
            if n0%1e6 == 0:
                print('{} processed'.format(n0))
            ## for bc+p7a
            s1 = p1.search(s)
            if s1:
                n1 += 1
                bc_seq = s1.group(1)
                bc.append(bc_seq)
                ## for i7
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
            t = '{}\t{}'.format(bc_seq, i7_seq, s)
            if save_seq and not t.startswith('-'):
                w.write(t+'\n')
    # wrap log
    ## p7 component
    p7c_fmt = '#{}\t{}\t{}\t{:.1f}\t{}'
    p7c_1 = p7c_fmt.format('p7a', n0, n1-n2-n3, (n1-n2-n3)/n0*100, fname)
    p7c_2 = p7c_fmt.format('p7a_p7b', n0, n2, n2/n0*100, fname)
    p7c_3 = p7c_fmt.format('p7b', n0, n3, n3/n0*100, fname)
    p7c_4 = p7c_fmt.format('no_ad', n0, n0-n1, (n0-n1)/n0*100, fname)
    ## top bc, i7
    idx_fmt = '>{}\t{}\t{}\t{}\t{:.1f}\t{}'
    top_bc = top_seq(bc, 3, revcmp=True)
    top_bc_list = [idx_fmt.format(i, 'bc', n0, n1, n1/n0*100, fname) for i in top_bc]
    top_i7 = top_seq(i7, 3)
    top_i7_list = [idx_fmt.format(i, 'i7', n0, n1, n1/n0*100, fname) for i in top_i7]
    ## message
    msg = '\n'.join([
        '-'*80,
        '{:<6}: {:>12}'.format('fastq', fx),
        '{:<6}: {:12,}'.format('total', n0),
        '-'*80,
        '\t'.join(["group", "total", "count", "pct", "sample"]),
        p7c_1,
        p7c_2,
        p7c_3,
        p7c_4,
        '-'*80,
        '\t'.join(['demx', 'demx_rc', 'demx_name', 'demx_n', 'demx_pct', 'group', 'total', 'i7', 'i7_pct', 'sample']),
        '\n'.join(top_bc_list),
        '\n'.join(top_i7_list),
        '-'*80,
    ])
    print(msg)
    with open(out_msg, 'wt') as w:
        w.write(msg+'\n')



class HiSeqP7(object):
    """"
    Check the P7 of HiSeq library
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
            'save_seq': False,
            'sub_sample': 0,
            'subdir': None,
            }
        self = update_obj(self, args_init, force=False)
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        if not isinstance(self.subdir, str):
            self.subdir = os.path.join(self.outdir, 'subset')
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
        # check target
        xname = fx_name(x, fix_pe=False)
        x_stat = os.path.join(self.outdir, xname+'.stat.log')
        if os.path.isfile(x_stat) and not self.overwrite:
            log.info('parse_i7 skipped, file exists: {}'.format(x_stat))
            with open(x_stat) as r:
                msg = r.readlines()
            print(msg)
        else:
            sub_fx = os.path.join(self.subdir, os.path.basename(x))
            if self.sub_sample > 0:
                Fastx(x).sample(out=sub_fx, n=self.sub_sample)
            else:
                symlink_file(x, sub_fx)
            parse_i7(sub_fx, self.outdir, self.save_seq)


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
        cmd_file = os.path.join(self.outdir, 'cmd.sh')
        cmd = ' '.join([
            shutil.which('Rscript'),
            qc_reportR,
            self.outdir,
            self.outdir])
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(report_html) and self.overwrite is False:
            log.info('file exists, skip generating html.')
        else:
            run_shell_cmd(cmd)
        if not os.path.exists(report_html):
            log.error('failed, generating html file')


    def run(self):
        check_path(self.subdir)
        self.run_multi(self.fq_list)
        self.report()

        
def get_args():
    parser = argparse.ArgumentParser(
        description='hiseq i7')
    parser.add_argument('-i', '--fq', nargs='+', required=True,
        help='FASTQ files, or directory contains fastq files')
    parser.add_argument('-o', '--outdir', default=None, required=True,
        help='The directory to save results.')
    parser.add_argument('-s', '--save-seq', dest='save_seq', 
        action='store_true',
        help='Saving the barcode and index sequences')
    parser.add_argument('--sub-sample', dest='sub_sample', type=int, 
        default=0,
        help='Sub sample the input fastq, default: [0], total reads')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    return parser
        

def main():
    args = vars(get_args().parse_args())
    HiSeqP7(**args).run()

    
if __name__ == '__main__':
    main()
        
# 



# ## patterns
# def parse_i7(fx, outdir, save_seq=False):
#     fname = os.path.basename(fx)
#     fname = fname.replace('.fq.gz', '')
#     if not os.path.isdir(outdir):
#         os.makedirs(outdir)
#     ## output
#     out_bc = os.path.join(outdir, fname + '.barcode.txt')
#     out_msg = os.path.join(outdir, fname + '.stat.log')
#     ## adapters
#     # p1 = re.compile('AGATCG[ACGTN]{22}AGTCAC')  # p7a
#     p1 = re.compile('([ACGTN]{6})AGATCG[ACGTN]{22}AGTCAC([ACGTN]{6})(.*)') # p7a+p7 
#     p2 = re.compile('AGATCG[ACGTN]{22}AGTCAC([ACGTN]{6})ATCTCG[ACGT]{12}TGCTTG') # p7a+p7+p7b
#     p3 = re.compile('([ACGTN]{6})ATCTCG[ACGT]{12}TGCTTG') # p7+p7b
#     ## count reads
#     n0 = 0
#     n1 = 0
#     n2 = 0
#     n3 = 0
#     bc = []
#     i7 = []
#     with open(out_bc, 'wt') as w:
#         for _,s,_,_ in pyfastx.Fastx(fx):
#             t = []
#             n0 += 1
#             if n0%1e6 == 0:
#                 print('{} processed'.format(n0))            
#             s1 = p1.search(s)
#             s3 = p3.search(s)
#             if s1:
#                 bc.append(s1.group(1))
#                 i7.append(s1.group(2))
#                 t.extend([s1.group(1), s1.group(2), s1.group(3)])
#                 n1 += 1
#                 s2 = p2.search(s)        
#                 if s2:
#                     n2 += 1
#                     n3 += 1
#             elif s3:
#                 t.append(s3.group(1))
#                 n3 += 1
#             if len(t) > 0 and save_seq:
#                 w.write('\t'.join(t)+'\n')
#     ## top i7
#     kbc = '\n'.join(['>{}\tbc\t{}\t{}\t{:.1f}\t{}'.format(i, n0, n1, n1/n0*100, fname) for i in top_seq(bc, 3, revcmp=True)])
#     ki7 = '\n'.join(['>{}\ti7\t{}\t{}\t{:.1f}\t{}'.format(i, n0, n1, n1/n0*100, fname) for i in top_seq(i7, 3)])
#     ## message
#     msg = '\n'.join([
#         '-'*80,
#         '{:>11}: {:>12}'.format('fastq', fx),
#         '{:>11}: {:12,}'.format('total', n0),
#         '\t'.join(["group", "total", "count", "pct", "sample"]),
#         '#{}\t{}\t{}\t{:.1f}\t{}'.format('P7a_i7', n0, n1, n1/n0*100, fname),
#         '#{}\t{}\t{}\t{:.1f}\t{}'.format('P7a_7b', n0, n2, n2/n0*100, fname),
#         '#{}\t{}\t{}\t{:.1f}\t{}'.format('i7_P7b', n0, n3, n3/n0*100, fname),
#         '-'*80,
#         '\t'.join(['demx', 'demx_rc', 'demx_name', 'demx_n', 'demx_pct', 'group', 'total', 'i7', 'i7_pct', 'sample']),
#         ki7,
#         kbc,
# #         '[{:<}]: \n{}'.format('i7', ki7),
# #         '[{:<}]: \n{}'.format('barcode', kbc),
#         '-'*80,
#     ])
#     print(msg)
#     with open(out_msg, 'wt') as w:
#         w.write(msg+'\n')

