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
import shutil
import pathlib
import argparse
from xopen import xopen
from multiprocessing import Pool
from contextlib import ExitStack
import Levenshtein as lev # distance

import hiseq
from hiseq.utils.utils import update_obj, log, Config, run_shell_cmd
from hiseq.utils.file import check_file, check_path, file_abspath, \
    list_file, list_fx, fx_name, symlink_file, file_exists
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
            'sample_sheet': None,
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
        if not os.path.isfile(self.sample_sheet):
            raise ValueError('file not exists x={}'.format(self.sample_sheet))
        if not os.path.isdir(self.datadir):
            raise ValueError('not a directory: {}'.format(self.datadir))
        if not isinstance(self.outdir, str):
            self.outdir = 'results'
        self.outdir = file_abspath(self.outdir)
        self.datadir = file_abspath(self.datadir)
        self.sample_sheet = file_abspath(self.sample_sheet)
        # check input fastq files
        self.raw_fq_list = list_fx(self.datadir, recursive=True)
        if len(self.raw_fq_list) < 2:
            raise ValueError('no fastq files in: {}'.format(self.datadir))
        self.init_files()
        # read table
        self.load_table()
        # fix files
        self.read_count_json = os.path.join(self.outdir, 'read_count.json')
        self.demx_report = os.path.join(self.outdir, 'report.txt')
        # save part of the objects
        # bc_ids, i7_ids, sheet, # remove the classes
        dd = self.__dict__.copy()
        for i in ['bc_ids', 'i7_ids', 'd_smp', 'sheet']:
            dd.pop(i)
        Config().dump(dd, self.config_yaml)


    def init_files(self):
        self.config_dir = os.path.join(self.outdir, 'config')
        self.config_yaml = os.path.join(self.config_dir, 'config.yaml')
        self.report_dir = os.path.join(self.outdir, 'report')
        self.hiseq_type = 'demx_r1'
        check_path([self.config_dir, self.report_dir])

        
    def load_table(self):
        self.sheet = SampleSheet(x=self.sample_sheet, outdir=self.outdir) # load table
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
        # save read count
        self.count_dir = os.path.join(self.outdir, 'read_count')
        check_path(self.count_dir)
        q_size = {}
        for fq in self.raw_fq_list:
            is_r1 = re.search('_(R)?1.f(ast)?q+.gz', fq, re.IGNORECASE)
            suffix = '_1.fq.gz' if is_r1 else '_2.fq.gz'
            prefix = fx_name(fq, fix_pe=True)
            ###################################################################
            # fix Index name: 
            # TruSeq_Index[1-48]
            # Next_Ad2.[1-24]
            # P7_[1-8][AB]
            #
            # Fix sample names: !!!!
            # fix for MGI, Next_Ad2.1 <- Next_Ad2-1_raw
            # fix for MGI, Next-Ad2.1 <- Next-Ad2_1_raw
            p1 = re.compile('TruSeq.Index(\d+)', flags=re.IGNORECASE)
            p2 = re.compile('Next.Ad2.(\d+)', flags=re.IGNORECASE)
            p3 = re.compile('P7_([\d][AB])', flags=re.IGNORECASE)
            m1 = p1.match(prefix)
            m2 = p2.match(prefix)
            m3 = p3.match(prefix)
            if m1:
                prefix = 'TruSeq_Index{}'.format(m1.group(1))
            elif m2:
                prefix = 'Next_Ad2.{}'.format(m2.group(1))
            elif m3:
                prefix = 'P7_{}'.format(m3.group(1))
            else:
                pass
            # output filenames
            s_name = self.d_smp.get(prefix, prefix)
            ###################################################################
            s_file = os.path.join(self.outdir, s_name+suffix) # target file
            if s_name and prefix in self.i7_ids:
                if os.path.isfile(s_file):
                    log.warning('renaming skipped, file exists: {}'.format(s_file))
                else:
                    symlink_file(fq, s_file)
                # count fq
                # save read count to file: read_count/s_name.count.json
                s_count_json = os.path.join(self.count_dir, s_name + '.count.json')
                if is_r1:
                    try:
                        if os.path.isfile(self.read_count_json):
                            n_size = Config().load(self.read_count_json)
                            n_fq = n_size.get(s_name, 0)
                        elif file_exists(s_count_json):
                            n_fq = Config().load(s_count_json).get(s_name, 0)
                        else:
                            n_fq = Fastx(s_file).number_of_seq()
                    except OSError as e:
                        print(e)
                        n_fq = 0
                    # save to file
                    s_count_d = {s_name: n_fq}
                    Config().dump(s_count_d, s_count_json)
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
                if self.demo:
                    t = os.path.join(bc_dir, 'demo', 'read_count.json')
                else:
                    t = os.path.join(bc_dir, 'read_count.json')
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
        n_exp = exp_df['reads'].sum()
        # to json
        self.rename_i7_files() #
        self.rename_bc_files() #
        q_size = self.i7_size
        q_size.update(self.bc_size) # sample_name:reads
        q_size = dict(sorted(q_size.items(), key=lambda x:x[0])) # sort by key
        Config().dump(q_size, self.read_count_json)
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
            '{} : {:>10} {:6.1f}M'.format('Expect reads', '', n_exp),
            '{} : {:>10} {:6.1f}M (pct: {:5.1f}%)'.format(
                'output reads', total, total/1e6, total/1e6/n_exp*100),
            '{:>3} {:<40s} {:>10} {:6} {:8s} {:8s}'.format(
                'order', 'filename', 'count', 'percent', 'million', 'design'),
            '\n'.join(f_stat),
            '-'*80,
        ])
        with open(self.demx_report, 'wt') as w:
            w.write(msg+'\n')
        print(msg)

        
    def report(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        hiseq_report_R = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        report_stdout = os.path.join(self.report_dir, 'report.stdout')
        report_stderr = os.path.join(self.report_dir, 'report.stderr')
        hiseq_report_html = os.path.join(
            self.report_dir,
            'HiSeq_report.html')
        cmd = ' '.join([
            '{}'.format(shutil.which('Rscript')),
            hiseq_report_R,
            self.outdir,
            self.report_dir,
            '1>{}'.format(report_stdout),
            '2>{}'.format(report_stderr),
            ])
        # save command
        cmd_txt = os.path.join(self.report_dir, 'report.cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # report_html
        if file_exists(hiseq_report_html):
            log.info('report() skipped, file exists: {}'.format(
                hiseq_report_html))
        else:
            run_shell_cmd(cmd)


    def run(self):
        log.info('Demulplexing starting')
        self.demx_barcode()
        self.wrap_dir()
        self.report()
        log.info('Demulplexing finish')

        

def get_args():
    """
    Demultiplexing, multi barcode files
    """
    parser = argparse.ArgumentParser(description='hiseq demx2')
    parser.add_argument('-s', '--xlsx-table', dest='sample_sheet', required=True,
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
    parser.add_argument('-x', '--in-read2', dest='in_read2',
        action='store_true', help='barcode in read2, default in read1, ...... demx2')
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