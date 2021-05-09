#!/usr/bin/env python3


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
from hiseq.utils.helper import update_obj, check_file, check_path, file_abspath, combinations, log, Config, listfile, fq_name, file_symlink
from hiseq.utils.seq import Fastx
from .demx_index import DemxIndex
from .demx_barcode import DemxBarcode
from .read_index import IndexTable
from .sample_sheet import SampleSheet


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
        f1_list = listfile(p7_dir, '*.fq.gz')
        f1_list = [i for i in f1_list if fq_name(i, pe_fix=True) in self.samples]
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
            f2_list = listfile(bc_dir, '*.fq.gz')
            f2_list = [i for i in f2_list if fq_name(i, pe_fix=True) in self.samples]
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
        f1 = listfile(self.datadir, '*.fq.gz', recursive=True)
        f2 = listfile(self.datadir, '*.fastq.gz', recursive=True)
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
            prefix = fq_name(fq, pe_fix=True)
            prefix = re.sub('_raw$', '', prefix) 
            s_name = self.d_smp.get(prefix, None) # i7_index -> sample_name
            if s_name and prefix in self.bc_ids:
                bc_fq[prefix] = bc_fq.get(prefix, [])
                bc_fq[prefix].append(fq)
        return bc_fq


#     def demx_barcode(self):
#         cmd_list = []
#         bc_list = self.sheet.to_barcode_table() # barcode.csv
#         for k in bc_list:
#             k = os.path.abspath(k)
#             k_id = os.path.basename(k)
#             k_id = re.sub('^i7_index.|.csv$', '', k_id)
#             k_fq = self.bc_fq.get(k_id, []) # [r1, r2]
#             k_fq = sorted(k_fq)
#             k_outdir = os.path.join(self.outdir, k_id)
#             k_log = os.path.join(k_outdir, 'demx.log')
#             s_name = self.d_smp.get(k_id, None)
#             s_file = os.path.join(self.outdir, s_name+'_1.fq.gz')
#             if os.path.isfile(s_file):
#                 print('barcode done: {}'.format(k_id))
#             else:
#                 check_path(k_outdir)
#                 if len(k_fq) == 2:
#                     cmd = 'hiseq demx -o {} -l 0 -r 1 -1 {} -2 {} -s {} 1> {}'.format(
#                         k_outdir, k_fq[0], k_fq[1], k, k_log)
#                     cmd_list.append(cmd)
#                 else:
#                     print('fq files not found: {}'.format(k))
#         with open(self.demx_bc_shell, 'wt') as w:
#             w.write('\n'.join(cmd_list)+'\n')
#         # run command
#         threads = 8 if len(cmd_list) > 8 else len(cmd_list)
#         if len(cmd_list) > 0:
#             with Pool(processes=threads) as pool:
#                 pool.map(run_shell_cmd, cmd_list)


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
            prefix = fq_name(fq, pe_fix=True)
            ###################################################################
            # Fix sample names:
            # fix for MGI, Next_Ad2.1 <- Next_Ad2-1_raw
            # fix for MGI, Next-Ad2.1 <- Next-Ad2_1_raw
            prefix = re.sub('_raw$', '', prefix)
            prefix = re.sub('Next.Ad', 'Next_Ad', prefix)
            prefix = re.sub('Ad2[^0-9]', 'Ad2.', prefix)
            s_name = self.d_smp.get(prefix, None) # i7_index -> sample_name
            ###################################################################
            s_file = os.path.join(self.outdir, s_name+suffix) # target file
            if s_name and prefix in self.i7_ids:
                if os.path.isfile(s_file):
                    log.warning('renaming skipped, file exists: {}'.format(s_file))
                else:
                    # os.rename(fq, s_file)
                    file_symlink(fq, s_file)
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
                bc_files = listfile(bc_dir, '*.fq.gz')
                for q in bc_files:
                    q_name = fq_name(q)
                    q_new = os.path.join(self.outdir, os.path.basename(q))
                    if q_name.startswith('undemx_'):
                        continue
                    else:
                        if os.path.isfile(q_new):
                            log.info('renaming skipped, file exists: {}'.format(q_new))
                        else:
                            # os.rename(q, q_new)
                            file_symlink(q, q_new)
                    # read count
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

