#!/usr/bin/env python

"""Processing samples_sheet

Example (xlsx)

sheet_name: sample_sheet
header:
    1.SampleID
    2.Lib_number
    3.Lib_user
    4.Sample_name
    5.RBP
    6.Cell_line
    7.Species
    8.Spike-in
    9.P7_index_id
    10.Barcode_id
    11.Seq_type
    12.Lib_type
    13.Reads,M
    14.FCID
    15.Lane
    16.Reference
    17.Spikein_ref
    18.P7_index
    19.Barcode_seq

Mission:
1. Sample_name, sanitize, unique
2. P7_index_id, Barcode_id, unique
3. index_table: p7,p5,bc
4. index_table: p7,bc (sub_table)
5. p7 to name
6. reads(M)

to_MGI:
p7_id,p7_seq,reads(Gb)

"""


import os
import sys
import re
import hiseq
import pathlib
import pandas as pd
from hiseq.utils.helper import update_obj

# mission-1
# 1. convert to MGI table
# 2. split by i7 index


class HiSeqIndex(object):
    """
    Example:

    >>> from sample_sheet import HiSeqIndex
    ## input index_seq
    >>> i = HiSeqIndex('ATCACG')
    >>> i.name
    'TruSeq_Index1'
    >>> i.index
    'ATCACG'

    ## input index_name
    >>> i = HiSeqIndex('TruSeq_Index1')
    >>> i.name
    'TruSeq_Index1'
    >>> i.index
    'ATCACG'

    ## input invalid
    """
    def __init__(self, x):
        x = x
        self.df = self.load_hiseq_index() # name:seq
        self.df2 = {v:k for k,v in self.df.items()}
        if isinstance(x, str):
            self.name = self.get_name(x)
            self.index = self.get_index(x)
            self.is_valid = x in self.df or x in self.df2
        elif isinstance(x, list):
            self.name = [self.get_name(i) for i in x]
            self.index = [self.get_index(i) for i in x]
            self.is_valid = [i in self.df or x in self.df2 for i in x]
        else:
            raise ValueError('illegal x={}, expect str,list got {}'.format(
                x, type(x).__name__))


    def load_hiseq_index(self):
        """
        The index table save in hiseq module
        Read the index_id and index of Illumina: TruSeq, Nextera
        and barcode
        """
        d = os.path.dirname(hiseq.__file__)
        f = os.path.join(d, 'bin', 'illumina_index.csv')
        idx = {}
        with open(f) as r:
            for line in r:
                if line.startswith('#'):
                    continue
                name,seq = line.strip().split(',')
                idx.update({name:seq})
        return idx


    def is_name(self, x):
        return x in self.df


    def is_index(self, x):
        return x in self.df2


    def get_name(self, x):
        if self.is_name(x):
            out = x
        elif self.is_index(x):
            out = self.df2.get(x, None)
        else:
            out = None
        return out


    def get_index(self, x):
        if self.is_index(x):
            out = x
        elif self.is_name(x):
            out = self.df.get(x, None)
        else:
            out = None
        return out


class SampleSheet(object):
    """
    Processing sample table file.xlsx

    ## Example
    >>> s = SampleSheet(x='YY00.xlsx', outdir='data')
    >>> s.to_MGI_table()
    >>> s.to_barcode_table()
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.df = self.load_xlsx()


    def init_args(self):
        args_init = {
            'x': None,
            'outdir': None,
        }
        self = update_obj(self, args_init, force=False)
        if not os.path.isfile(self.x):
            raise ValueError('file not exists: {}'.format(self.x))
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        self.xlsx_prefix = os.path.splitext(os.path.basename(self.x))[0]


    def load_xlsx(self):
        """
        Read the table, save as DataFrame
        sanitize the filename: by "_"
        """
        c = ['Sample_name*', 'P7_index_id*', 'Barcode_id*', 'Reads, M']
        try:
            df = pd.read_excel(self.x, sheet_name='sample_sheet')
            df = df.dropna(thresh=2)
            df = df[c]
            df.columns = ['name', 'i7', 'bc', 'reads']
            df['name'] = self.sanitize(df['name'].to_list()) # sanitize
        except:
            df = pd.DataFrame(columns=['name', 'i7', 'bc', 'reads'])
        return df


    def sanitize(self, x):
        """
        The filename rules:
        1. [A-Za-z0-9_-.]
        2. convert to '_'
        3. unique
        """
        if isinstance(x, str):
            out = re.sub('[^A-Za-z0-9_.\-]', '_', x)
            out = re.sub('_+', '_', out)
        elif isinstance(x, list):
            out = [self.sanitize(i) for i in x]
        else:
            log.warning('illegal x')
            out = x
        return out


    def to_MGI_table(self):
        """
        Convert to MGI table
        1. Assign the names by i7_index_id
        2. concatenate files with same P7 index
        3. Sum the reads
        """
        df2 = self.df.groupby('i7').sum()
        df2 = df2.assign(i7_seq = HiSeqIndex(df2.index.to_list()).index)
        df2 = df2.loc[:, ['i7_seq', 'reads']]
        mgi_csv = os.path.join(self.outdir, self.xlsx_prefix + '.MGI.csv')
        df2.to_csv(mgi_csv)
        return df2


    def to_barcode_table(self):
        """
        Generate index table for barcode:
        format:
        filename,i7,i5,bc
        sample1,NULL,NULL,CCTATA

        filename:
        i7_index.csv
        """
        bc_tables = []
        s = self.df.groupby('i7').size()
        s = s[s>1]
        i7_list = s.index.to_list()
        if len(i7_list) > 0:
            for i7 in i7_list:
                i7_csv = os.path.join(self.outdir, 'i7_index.{}.csv'.format(i7))
                df1 = self.df.loc[self.df['i7'] == i7]
                # convert name to index
                df2 = df1.assign(
                    bc_seq = HiSeqIndex(df1['bc'].to_list()).index,
                    i7_seq = 'NULL', # HiSeqIndex(df1['i7'].to_list()).index,
                    i5_seq = 'NULL',
                )
                df2 = df2.loc[:, ['name', 'i7_seq', 'i5_seq', 'bc_seq']]
                df2.to_csv(i7_csv, index=False, header=False)
                bc_tables.append(i7_csv)
        return bc_tables



def main():
    if len(sys.argv) == 3:
        print('mode=1, convert to MGI table')
        x_table, outdir = sys.argv[1:3]
        p = SampleSheet(x=x_table, outdir=outdir)
        p.to_MGI_table()
        p.to_barcode_table()
    else:
        msg = '\n'.join([
            'Usage: ',
            'sample_sheet.py <x.xlsx> <outdir>',           
        ])
        print(msg)
        sys.exit(1)


if __name__ == '__main__':
    main()

#



# def demx_barcode(x, datadir, outdir):
#     """
#     Arguments
#     ---------
#     table:   xlsx
#     datadir: download data
#     outdir:  target

#     Mission:
#     1. demx barcode:
#     """
#     p = SampleSheet(x=x, outdir=outdir) # load table
#     df = p.df # 'name', 'i7', 'bc', 'reads'
#     dfx = df.loc[:, ['name', 'i7']].set_index('i7').to_dict('dict')
#     d_smp = dfx['name'] # dict, {index:sample_name}
#     s = df.groupby('i7').size()
#     bc_ids = s[s>1]  # i7_index, with barcode
#     bc_fq = {} # fastq files for barcode demultiplexing
#     f1 = listfile(datadir, '*.fq.gz', recursive=True)
#     f2 = listfile(datadir, '*.fastq.gz', recursive=True)
#     fq_list = f1 + f2
#     for fq in fq_list:
#         prefix = fq_name(fq, pe_fix=True)
#         s_name = d_smp.get(prefix, None)
#         if s_name and prefix in bc_ids:
#             bc_fq[prefix] = bc_fq.get(prefix, [])
#             bc_fq[prefix].append(fq)
#     # command for demx
#     # hiseq demx -1 {} -2 {} -o {} -s {} -l 0 -r 1 -p 2
#     cmd_list = []
#     bc_list = p.to_barcode_table() # barcode.csv
#     for k in bc_list:
#         k = os.path.abspath(k)
#         k_id = os.path.basename(k)
#         k_id = re.sub('^i7_index.|.csv$', '', k_id)
#         k_fq = bc_fq.get(k_id, []) # [r1, r2]
#         k_fq = sorted(k_fq)
#         k_outdir = os.path.join(outdir, k_id)
#         k_log = os.path.join(k_outdir, 'demx.log')
#         s_name = d_smp.get(k_id, None)
#         s_file = os.path.join(outdir, s_name+'_1.fq.gz')
#         if os.path.isfile(s_file):
#             print('barcode done: {}'.format(k_id))
#         else:
#             check_path(k_outdir)
#             # os.makedirs(k_outdir)
#             if len(k_fq) == 2:
#                 cmd = 'hiseq demx -o {} -l 0 -r 1 -1 {} -2 {} -s {} 1> {}'.format(
#                     k_outdir, k_fq[0], k_fq[1], k, k_log)
#                 cmd_list.append(cmd)
#             else:
#                 print('fq files not found: {}'.format(k))
#     cmd_file = os.path.join(outdir, 'demx_bc.sh')
#     with open(cmd_file, 'wt') as w:
#         w.write('\n'.join(cmd_list)+'\n')
#     # demx
#     threads = 8 if len(cmd_list) > 8 else len(cmd_list)
#     if len(cmd_list) > 0:
#         with Pool(processes=threads) as pool:
#             pool.map(run_shell_cmd, cmd_list)


# def wrap_mgi_dir(x, datadir, outdir, threads=4):
#     """
#     1. count reads for i7-only
#     2. move barcode files to outdir
#     3. read count
#     """
#     ############
#     # i7 files #
#     ############
#     datadir = os.path.abspath(datadir)
#     outdir = os.path.abspath(outdir)
#     p = SampleSheet(x=x, outdir=outdir) # load table
#     df = p.df # 'name', 'i7', 'bc', 'reads'
#     dfx = df.loc[:, ['name', 'i7']].set_index('i7').to_dict('dict')
#     d_smp = dfx['name'] # dict, {index:sample_name}
#     s = df.groupby('i7').size()
#     i7_ids = s[s==1] # i7 only
#     bc_ids = s[s>1]  # i7_index, with barcode
#     bc_fq = {} # fastq files for barcode demultiplexing
#     f1 = listfile(datadir, '*.fq.gz', recursive=True)
#     f2 = listfile(datadir, '*.fastq.gz', recursive=True)
#     fq_list = f1 + f2
#     for fq in fq_list:
#         is_r1 = re.search('_1.f(ast)?q+.gz', fq, re.IGNORECASE)
#         suffix = '_1.fq.gz' if is_r1 else '_2.fq.gz'
#         prefix = fq_name(fq, pe_fix=True)
#         sname = d_smp.get(prefix, None)
#         if sname and prefix in i7_ids:
#             sfile = os.path.join(outdir, sname+suffix) # target file
#             if os.path.isfile(sfile):
#                 print('renaming skipped, file exists: {}'.format(sfile))
#             else:
#                 os.rename(fq, sfile)
#         # barcode, demultiplexing
#         if sname and prefix in bc_ids:
#             bc_fq[prefix] = bc_fq.get(prefix, [])
#             bc_fq[prefix].append(fq)
#     ##############
#     # i7 count   #
#     ##############
#     q_size = {}
#     for i in listfile(outdir, '*_1.fq.gz', recursive=False):
#         iname = fq_name(i, pe_fix=True)
#         icount = Fastx(i).number_of_seq()
#         q_size.update({iname:icount})
#     ##############
#     # bc count   #
#     ##############
#     n_undemx = 0
#     if len(bc_ids) > 0:
#         for i in bc_ids.index.to_list():
#             bc_dir = os.path.join(outdir, i)
#             bc_files = listfile(bc_dir, '*.fq.gz')
#             for q in bc_files:
#                 q_name = fq_name(q)
#                 q_new = os.path.join(outdir, os.path.basename(q))
#                 if q_name.startswith('undemx_'):
#                     continue
#                 else:
#                     if os.path.isfile(q_new):
#                         print('renaming skipped, file exists: {}'.format(q_new))
#                     else:
#                         os.rename(q, q_new)
#             # read count
#             t = os.path.join(bc_dir, 'read_count.toml')
#             try:
#                 d = Config().load(t)
#                 n_undemx += d.get('undemx', 0)
#                 d.pop('undemx')
#                 q_size.update(d)
#             except:
#                 print('file not exists: {}'.format(t))
#     ## update undemx
#     q_size.update({'undemx': n_undemx})
#     ##############
#     # all count  #
#     ##############
#     stat_toml = os.path.join(outdir, 'read_count.toml')
#     Config().dump(q_size, stat_toml)
#     # log, msg
#     stat_report = os.path.join(outdir, 'report.txt')
#     total = sum(q_size.values())
#     f_stat = []
#     i = 0
#     for k,v in q_size.items():
#         i += 1
#         s = '{:>3} {:<40s} {:>10,} {:6.2f}% {:8.1f}'.format(
#             i, k, v, v/total*100, v/1e6)
#         f_stat.append(s)
#     # message
#     msg = '\n'.join([
#         '-'*80,
#         'Demultiplex report:',
#         '{} : {:>10} {:6.1f}M'.format('Input reads', total, total/1e6),
#         '{:>3} {:<40s} {:>10} {:6} {:8s}'.format(
#             'order', 'filename', 'count', 'percent', 'million'),
#         '\n'.join(f_stat),
#         '-'*80,
#     ])
#     with open(stat_report, 'wt') as w:
#         w.write(msg+'\n')
#     print(msg)
