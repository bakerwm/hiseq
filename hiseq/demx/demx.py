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
from hiseq.demx.demx_index import DemxIndex
from hiseq.demx.demx_barcode import DemxBarcode
from hiseq.demx.read_index import IndexTable
from hiseq.demx.sample_sheet import SampleSheet


class Demx(object):
    """
    This script is designed for i7+barcode demx

    Arguments
    ---------
    index_table, str
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
        self.mission = self.get_mission()
        check_path(self.outdir)


    def init_args(self):
        args_init = {
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'index_table': None, #name,index
            'mismatch': 0,
            'in_read2': True,
            'barcode_n_left': 0,
            'barcode_n_right': 1,
            'overwrite': False,
            'demo': False,
            'gzipped': True,
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
        # absolute path
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.index_table = file_abspath(self.index_table)
        self.outdir = file_abspath(self.outdir)
        # index
        self.check_index()


    def check_index(self):
        p = IndexTable(index_table=self.index_table, mismatch=self.mismatch)
        self.samples = p.samples # name:(p7,p5,bc)
        self.p7_list = p.p7_list
        self.p5_list = p.p5_list
        self.bc_list = p.bc_list
        out = p.run()
        if not out:
            log.warning('index unknown: {}'.format(self.index_table))
#             raise ValueError('index illegal: {}'.format(self.index_table))


    def get_mission(self):
        """
        Determine the mission type: p7:4, p5:2, bc:1
        1. p7 only
        2. bc only
        3. p7 + bc
        4. p5 (optional)
        """
        has_p7 = len(self.p7_list) > 0
        has_bc = len(self.bc_list) > 0
        if has_p7 and has_bc:
            m = 5
        elif has_p7:
            m = 4
        elif has_bc:
            m = 1
        else:
            m = 0
        return m


    def demx_p7(self):
        """
        Demultiplex p7 only
        mission: 4
        """
        s = ['{},{}'.format(k, v[0]) for k,v in self.samples.items()]
        p7_idx_table = os.path.join(self.outdir, 'index_table.csv')
        with open(p7_idx_table, 'wt') as w:
            w.write('\n'.join(s)+'\n')
        if self.mission == 4:
            req_args = [
                'fq1', 'fq2', 'outdir', 'mismatch', 'demo', 'gzipped',
                'overwrite']
            args = {i:getattr(self, i, None) for i in req_args if hasattr(self, i)}
            args.update({
                'index_type': 'i7',
                'index_table': p7_idx_table,
            })
            DemxIndex(**args).run()


    def demx_bc(self):
        """
        Demultiplex bc only
        mission: 1
        """
        s = ['{},{}'.format(k, v[2]) for k,v in self.samples.items()]
        bc_idx_table = os.path.join(self.outdir, 'index_table.csv')
        with open(bc_idx_table, 'wt') as w:
            w.write('\n'.join(s)+'\n')
        if self.mission == 1:
            req_args = [
                'fq1', 'fq2', 'outdir', 'in_read2', 'mismatch',
                'barcode_n_left', 'barcode_n_right', 'demo', 'gzipped',
                'overwrite']
            args = {i:getattr(self, i, None) for i in req_args if hasattr(self, i)}
            args.update({'index_table': bc_idx_table})
            DemxBarcode(**args).run()


    def split_p7_bc(self):
        """
        For p7+bc mode
        split table into
        1. p7 only
        2. bc only
        """
        if self.mission < 5:
            log.warning('no need to split index_table')
            return None
        # full version of p7_list
        p7_list = [i[0] for i in self.samples.values()]
        # split by p7-index
        d1 = {} # level-2
        s2 = [] # level-1
        for k,v in self.samples.items():
            p7,p5,bc = v
            if p7_list.count(p7) > 1:
                if p7 in d1:
                    d1[p7].append('{},{}'.format(k, bc))
                else:
                    d1[p7] = ['{},{}'.format(k, bc)]
                    # add to level-1
                    s2.append('{},{}'.format(p7, p7))
            else:
                s2.append('{},{}'.format(k, p7))
        # save as file
        p7_dir = os.path.join(self.outdir, 'p7_index')
        p7_idx_table = os.path.join(p7_dir, 'index_table.csv')
        check_path(p7_dir)
        with open(p7_idx_table, 'wt') as w:
            w.write('\n'.join(s2)+'\n')
        # save bc table
        bc_idx_tables = []
        for k,v in d1.items():
            bc_dir = os.path.join(self.outdir, 'bc_index', k)
            check_path(bc_dir)
            bc_idx_table = os.path.join(bc_dir, 'index_table.csv')
            with open(bc_idx_table, 'wt') as w:
                w.write('\n'.join(v)+'\n')
            bc_idx_tables.append(bc_idx_table)
        return (p7_idx_table, bc_idx_tables)


    def demx_p7_bc(self):
        """
        Demultiplex both p7 and bc
        default order: p7->bc
        """
        if self.mission == 5:
            # for p7
            p7_table, bc_tables = self.split_p7_bc()
            p7_dir = os.path.join(self.outdir, 'p7_index')
            req_args = ['fq1', 'fq2', 'outdir', 'mismatch', 'demo', 'gzipped', 'overwrite']
            args = {i:getattr(self, i, None) for i in req_args if hasattr(self, i)}
            args.update({
                'outdir': p7_dir,
                'index_type': 'i7',
                'index_table': p7_table,
            })
            DemxIndex(**args).run()

            # for bc
            bc_args_list = []
            req_args = [
                'in_read2', 'mismatch', 'barcode_n_left', 'barcode_n_right', 'demo',
                'gzipped', 'overwrite']
            bc_args = {i:getattr(self, i, None) for i in req_args if hasattr(self, i)}
            for bc_table in bc_tables:
                idir = os.path.dirname(bc_table)
                iname = os.path.basename(idir)
                bc_args2 = bc_args.copy()
                if self.is_pe:
                    fq1 = os.path.join(p7_dir, iname+'_1.fq.gz')
                    fq2 = os.path.join(p7_dir, iname+'_2.fq.gz')
                else:
                    fq1 = os.path.join(p7_dir, iname+'_1.fq.gz')
                    fq2 = None
                bc_args2.update({
                    'fq1': fq1,
                    'fq2': fq2,
                    'outdir': os.path.dirname(bc_table),
                    'index_table': bc_table,
                })
                bc_args_list.append(bc_args2)
                DemxBarcode(**bc_args2).run()

            # organize files and stat
            self.wrap_p7_bc()


    def wrap_p7_bc(self):
        """
        1. organize fastq files
        2. report.txt
        """
        undemx = 0
        p7_table, bc_tables = self.split_p7_bc()
        # level-1: p7 files
        p7_dir = os.path.join(self.outdir, 'p7_index')
        f1_list = list_file(p7_dir, '*.fq.gz')
        f1_list = [i for i in f1_list if fx_name(i, fix_pe=True) in self.samples]
        f1x_list = [os.path.join(self.outdir, os.path.basename(i)) for i in f1_list]
        for f1, f1x in zip(f1_list, f1x_list):
            if os.path.exists(f1) and not os.path.exists(f1x):
                os.rename(f1, f1x)
        # for read count
        p7_stat_toml = os.path.join(p7_dir, 'read_count.toml')
        p7_stat = Config().load(p7_stat_toml)
        df_stat = {k:v for k,v in p7_stat.items() if k in self.samples}
        undemx += p7_stat.get('undemx', 0)
        # level-2: bc files
        for bc_table in bc_tables:
            bc_dir = os.path.dirname(bc_table)
            f2_list = list_file(bc_dir, '*.fq.gz')
            f2_list = [i for i in f2_list if fx_name(i, fix_pe=True) in self.samples]
            f2x_list = [os.path.join(self.outdir, os.path.basename(i)) for i in f2_list]
            for f2, f2x in zip(f2_list, f2x_list):
                if os.path.exists(f2) and not os.path.exists(f2x):
                    os.rename(f2, f2x)
            # for read count
            bc_stat_toml = os.path.join(bc_dir, 'read_count.toml')
            bc_stat = Config().load(bc_stat_toml)
            df_stat.update({k:v for k,v in bc_stat.items() if k in self.samples})
            undemx += bc_stat.get('undemx', 0)
        # update undemx
        df_stat['undemx'] = undemx
        # save to file
        df_stat = dict(sorted(df_stat.items(), key=lambda x:x[0])) # sort by key
        stat_toml = os.path.join(self.outdir, 'read_count.toml')
        Config().dump(df_stat, stat_toml)
        # report
        stat_report = os.path.join(self.outdir, 'report.txt')
        total = sum(df_stat.values())
        scale = 4e8/total # to_400M
        f_stat = []
        i = 0
        for k,v in df_stat.items():
            i += 1
            s = '{:>3} {:<40s} {:>10,} {:6.2f}% {:8.1f} {:8.1f}'.format(
                i, k, v, v/total*100, v/1e6, v/1e6*scale)
            f_stat.append(s)
        # message
        msg = '\n'.join([
            '-'*80,
            'Demultiplex report:',
            '{} : {:>10} {:6.1f}M'.format('Input reads', total, total/1e6),
            '{:>3} {:<40s} {:>10} {:6} {:8s} {:8s}'.format(
                'order', 'filename', 'count', 'percent', 'million', 'to_400M'),
            '\n'.join(f_stat),
            '-'*80,
        ])
        with open(stat_report, 'wt') as w:
            w.write(msg+'\n')
        print(msg)


    def run(self):
        log.info('Demulplexing starting')
        if self.mission == 1:
            self.demx_bc()
        elif self.mission == 4:
            self.demx_p7()
        elif self.mission == 5:
            self.demx_p7_bc()
        else:
            pass
        log.info('Demulplexing finish')


def get_args():
    parser = argparse.ArgumentParser(description='hiseq demx')
    parser.add_argument('-1', '--fq1', required=True,
        help='read1 in fastq format, gzipped')
    parser.add_argument('-2', '--fq2',
        help='read2 in fastq format, gzipped, (optional)')
    parser.add_argument('-o', '--outdir', required=True,
        help='directory to save the reulsts')
    parser.add_argument('-s', '--index-table', dest='index_table', required=True,
        help='index list in csv format, [filename,i7,i5,barcode]')
    parser.add_argument('--demo', action='store_true',
        help='run demo (1M reads) for demostration, default: off')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]')
    parser.add_argument('-x', '--barcode-in-read2', action='store_true',
        help='barcode in read2')
    parser.add_argument('-l', '--barcode-n-left', type=int, dest='barcode_n_left',
        default=0, help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int, dest='barcode_n_right',
        default=0, help='bases locate on the right of barcode')
    parser.add_argument('-p', '--threads', type=int, default=1,
        help='number of threads, default: [1]')
    parser.add_argument('-j', '--parallel-jobs', type=int, dest='parallel_jobs',
        default=1, help='number of josb run in parallel, default: [1]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    return parser


def main():
    args = vars(get_args().parse_args())
    Demx(**args).run()


if __name__ == '__main__':
    main()

#