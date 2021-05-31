#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# to-do: in-case, 0-count files
"""
Description:
1. P7 index from comment section in fastq: @fqname  comment
2. inline barcode in the sequence, position1, 2
3. SE or PE reads

mission-1:
1. P7 barcode:
2. inline-barcode:


utils:
1. bc reader
2. errors
3. errors: levenshtein distance


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
from multiprocessing import Pool
from contextlib import ExitStack
import Levenshtein as lev # distance
from hiseq.utils.utils import update_obj, log, Config
from hiseq.utils.file import check_file, check_path, file_abspath, \
    list_file, fx_name, symlink_file
from hiseq.utils.seq import Fastx
from hiseq.demx.demx import Demx
from hiseq.demx.demx_index import DemxIndex
from hiseq.demx.demx_barcode import DemxBarcode
from hiseq.demx.read_index import IndexTable
from hiseq.demx.sample_sheet import SampleSheet


class Demx2(object):
    """
    This script is designed for barcode demx

    Arguments
    ---------
    x, str
        The table of sample list, could be: .csv, or .xlsx
        required columns:
        1. xlsx file, required columns
        ['Sample_name*', 'P7_index_id*', 'Barcode_id*', 'Reads, M']
        2. csv file
        ['sample_name', 'i7_id', 'i5_id', 'bc_id', 'reads']

    datadir, str
        The path to the fastq files

    outdir, str
        The path to the dir, final output

    Description
    -----------
    The fastq files in datadir are named by the i7_index_name; in case, some
    of the i7_index file contains multiple sub_files, distinguished by
    in-line barcode
    This function is designed to do:
    1. rename i7_only files, (retrieve sample_name from table, by i7_index_name)
    2. demultiplex the i7 files, contains barcode
    3. organize the report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'x': None,
            'datadir': 'from_illumina',
            'outdir': 'results',
            'mismatch': 0,
            'in_read2': True,
            'barcode_n_left': 0,
            'barcode_n_right': 1,
            'overwrite': False,
            'demo': False,
            'gzipped': True,
        }
        self = update_obj(self, args_init, force=False)
        if not os.path.isfile(self.x):
            raise ValueError('file not exists x={}'.format(self.x))
        if not os.path.isdir(self.datadir):
            raise ValueError('not a directory: {}'.format(self.datadir))
        if not isinstance(self.outdir, str):
            self.outdir = 'results'
        # check input fastq files
        f1 = list_file(self.datadir, '*.fq.gz', recursive=True)
        f2 = list_file(self.datadir, '*.fastq.gz', recursive=True)
        self.raw_fq_list = f1 + f2
        if len(self.raw_fq_list) < 2:
            raise ValueError('no fastq files in: {}'.format(self.datadir))
        # read table
        self.load_table()
        # fix files
        # self.demx_bc_shell = os.path.join(self.outdir, 'demx_bc.sh')
        self.read_count_toml = os.path.join(self.outdir, 'read_count.toml')
        self.demx_report = os.path.join(self.outdir, 'report.txt')


    def load_table(self):
        self.sheet = SampleSheet(x=self.x, outdir=self.outdir) # load table
        df = self.sheet.df # 'name', 'i7', 'bc', 'reads'
        dfx = df.loc[:, ['name', 'i7']].set_index('i7').to_dict('dict')
        self.d_smp = dfx['name'] # dict, {index:sample_name}
        s = df.groupby('i7').size()
        self.i7_ids = s[s==1] # i7 only
        self.bc_ids = s[s>1]  # i7_index, with barcode
        self.bc_fq = self.get_bc_raw_files() # fastq files for barcode


    def get_bc_raw_files(self):
        bc_fq = {}
        for fq in self.raw_fq_list:
            prefix = fx_name(fq, fix_pe=True)
            prefix = re.sub('_raw$', '', prefix)
            s_name = self.d_smp.get(prefix, None) # i7_index -> sample_name
            if s_name and prefix in self.bc_ids:
                bc_fq[prefix] = bc_fq.get(prefix, [])
                bc_fq[prefix].append(fq)
        return bc_fq

    
    def demx_barcode(self):
        bc_list = self.sheet.to_barcode_table() # barcode.csv
        threads = 8 if len(bc_list) > 8 else len(bc_list)
        if len(bc_list) > 0:
            with Pool(processes=threads) as pool:
                pool.map(self.demx_barcode_single, bc_list)


    def demx_barcode_single(self, k):
        k = os.path.abspath(k)
        k_id = os.path.basename(k)
        k_id = re.sub('^i7_index.|.csv$', '', k_id)
        k_fq = self.bc_fq.get(k_id, []) # [r1, r2]
        k_fq = sorted(k_fq)
        k_outdir = os.path.join(self.outdir, k_id)
        k_log = os.path.join(k_outdir, 'demx.log')
        s_name = self.d_smp.get(k_id, None)
        s_file = os.path.join(self.outdir, s_name+'_1.fq.gz')
        args = {
            'fq1': k_fq[0],
            'fq2': k_fq[1],
            'outdir': k_outdir,
            'index_table': k,
            'in_read2': self.in_read2,
            'barcode_n_left': self.barcode_n_left,
            'barcode_n_right': self.barcode_n_right,
            'overwrite': self.overwrite,
            'demo': self.demo,
            'gzipped': self.gzipped,
        }
        if os.path.isfile(s_file):
            log.info('barcode done: {}'.format(k_id))
        else:
            Demx(**args).run()


    def rename_i7_files(self):
        """
        in case:
        i7_name were sanitized (by "-"):
        eg: 2.1 -> 2-1
        """
        q_size = {}
        for fq in self.raw_fq_list:
            is_r1 = re.search('_1.f(ast)?q+.gz', fq, re.IGNORECASE)
            suffix = '_1.fq.gz' if is_r1 else '_2.fq.gz'
            prefix = fx_name(fq, fix_pe=True)
            ###################################################################
            # Fix sample names: !!!!
            # fix for MGI, Next_Ad2.1 <- Next_Ad2-1_raw
            # fix for MGI, Next-Ad2.1 <- Next-Ad2_1_raw
            prefix = re.sub('_raw$', '', prefix)
            prefix = re.sub('TruSeq.Index', 'TruSeq_Index', prefix, flags=re.IGNORECASE)
            prefix = re.sub('Next.Ad', 'Next_Ad', prefix, flags=re.IGNORECASE)
            # prefix = re.sub('Next.Ad', 'Next_Ad', prefix)
            prefix = re.sub('Ad2[^0-9]', 'Ad2.', prefix)
            s_name = self.d_smp.get(prefix, None) # i7_index -> sample_name
            ###################################################################
            s_file = os.path.join(self.outdir, s_name+suffix) # target file
            if s_name and prefix in self.i7_ids:
                if os.path.isfile(s_file):
                    log.warning('renaming skipped, file exists: {}'.format(s_file))
                else:
                    # os.rename(fq, s_file)
                    symlink_file(fq, s_file)
                # count fq
                if is_r1:
                    try:
                        if os.path.isfile(self.read_count_toml):
                            n_size = Config().load(self.read_count_toml)
                            n_fq = n_size.get(s_name, 0)
                        else:
                            n_fq = Fastx(s_file).number_of_seq()
                    except OSError as e:
                        print(e)
                        n_fq = 0
                    q_size.update({
                        s_name: n_fq
                    })
                    log.info('counting reads, {}, {}: {}'.format(
                            prefix, s_name, n_fq))
        # return q_size
        self.i7_size = q_size


    def rename_bc_files(self):
        q_size = {}
        n_undemx = 0
        if len(self.bc_ids) > 0:
            for i in self.bc_ids.index.to_list():
                bc_dir = os.path.join(self.outdir, i)
                bc_files = list_file(bc_dir, '*.fq.gz')
                for q in bc_files:
                    q_name = fx_name(q)
                    q_new = os.path.join(self.outdir, os.path.basename(q))
                    if q_name.startswith('undemx_'):
                        continue
                    else:
                        if os.path.isfile(q_new):
                            log.info('renaming skipped, file exists: {}'.format(q_new))
                        else:
                            # os.rename(q, q_new)
                            symlink_file(q, q_new)
                # read count file
                t = os.path.join(bc_dir, 'read_count.toml')
                try:
                    d = Config().load(t)
                    n_undemx += d.get('undemx', 0)
                    d.pop('undemx')
                    q_size.update(d)
                except:
                    log.warning('file not exists: {}'.format(t))
        # return q_size #
        self.bc_size = q_size
        self.n_undemx = n_undemx
        self.bc_size.update({'undemx': n_undemx})


    def wrap_dir(self):
        # Expect reads
        exp_df = self.sheet.df.loc[:, ['name', 'reads']].set_index('name')
        exp_size = exp_df.to_dict('dict')['reads'] # sample_name:reads
        # to toml
        self.rename_i7_files() #
        self.rename_bc_files() #
        q_size = self.i7_size
        q_size.update(self.bc_size) # sample_name:reads
        q_size = dict(sorted(q_size.items(), key=lambda x:x[0])) # sort by key
        Config().dump(q_size, self.read_count_toml)
        # to txt
        total = sum(q_size.values())
        if total < 1:
            total = 1
        f_stat = []
        i = 0
        for k,v in q_size.items():
            i += 1
            e = exp_size.get(k, 0) # million
            s = '{:>3} {:<40s} {:>10,} {:6.2f}% {:8.1f} {:8.1f}'.format(
                i, k, v, v/total*100, v/1e6, e)
            f_stat.append(s)
        # message
        msg = '\n'.join([
            '-'*80,
            'Demultiplex report:',
            '{} : {:>10} {:6.1f}M'.format('Input reads', total, total/1e6),
            '{:>3} {:<40s} {:>10} {:6} {:8s} {:8s}'.format(
                'order', 'filename', 'count', 'percent', 'million', 'design'),
            '\n'.join(f_stat),
            '-'*80,
        ])
        with open(self.demx_report, 'wt') as w:
            w.write(msg+'\n')
        print(msg)


    def run(self):
        log.info('Demulplexing starting')
        self.demx_barcode()
        self.wrap_dir()
        log.info('Demulplexing finish')

        

def get_args():
    """
    Demultiplexing, multi barcode files
    """
    parser = argparse.ArgumentParser(description='hiseq demx2')
    parser.add_argument('-s', '--xlsx-table', dest='x', required=True,
        help="Sample table in (xlsx|csv) format; xlsx: require the columns\
        ['Sample_name*', 'P7_index_id*', 'Barcode_id*', 'Reads, M']; \
        csv: require the columns: ['name', 'i7', 'i5', 'bc', 'reads'] \
        the csv file could be `hiseq sheet -s a.xlsx -o data` output: *.demx.csv")
    parser.add_argument('-d', '--datadir', dest='datadir', required=True,
        help='Directory saving the fastq files')
    parser.add_argument('-o', '--outdir', dest='outdir',
        help='directory to save the reulsts')
    parser.add_argument('--demo', action='store_true',
        help='run demo (1M reads) for demostration, default: off')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]')
    parser.add_argument('-x', '--barcode-in-read2', dest='barcode_in_read2',
        action='store_true', help='barcode in read2')
    parser.add_argument('-l', '--barcode-n-left', type=int, dest='barcode_n_left',
        default=0, help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int, dest='barcode_n_right',
        default=0, help='bases locate on the right of barcode')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    return parser


def main():
    args = vars(get_args().parse_args())
    Demx2(**args).run()


if __name__ == '__main__':
    main()

#