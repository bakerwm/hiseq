#!/usr/bin/env python3


"""
Description 

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
import argparse
import Levenshtein as lev # distance
from xopen import xopen
import logging
from multiprocessing import Pool
from hiseq.utils.helper import *


class IndexList(object):
    """
    Input the index for demultiplexing
    #file_name, #index1, #index2, #barcode

    criteria:
    1. file_name: unique
    2. index1 + index2 + barcode unique
    3. mismatch required
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.init_args()
        self.idx = self.load()
        self.sample = self.idx.get('sample', {})
        self.count = self.idx.get('count', {})
        self.main = self.idx.get('main', {})
        # sequences for each index
        m = [i.split(',') for i in self.sample.keys()]
        self.index1 = [i[0] for i in m]
        self.index2 = [i[1] for i in m]
        self.barcode = [i[2] for i in m]
        self.status = self.check()


    def init_args(self):
        """
        required arguments: csv/txt/excel
        """
        args_init = {
            'input': None,
            'outdir': os.getcwd(),
            'mismatch': 0
        }
        for k, v in args_init.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        if self.input is None:
            log.error('input, required')


    def load(self):
        """
        Read index info file
        csv file
        index1-index2-barcode
        """
        da = {}
        db = {}
        dx = {}
        with open(self.input) as r:
            for line in r:
                if line.startswith('#'):
                    continue # skip
                s = line.strip().split(',')
                if not len(s) == 4:
                    log.error('input: format unknown: fname,index1,index2,barcode')
                    flag = False
                    break
                fname, ia, ib, bc = s
                ia = ia.upper()
                ib = ib.upper()
                ic = bc.upper()
                n = ','.join([ia, ib, ic])
                da[n] = fname # file name
                db[n] = db.get(n, 0) + 1
                # save all
                dx[ia] = dx.get(ia, {})
                dx[ia][ib] = dx[ia].get(ib, {})
                if ic in dx[ia][ib]:
                    dx[ia][ib][ic] = [dx[ia][ib][ic], fname]
                else:
                    dx[ia][ib][ic] = fname
        # output
        return {
            'sample': da,
            'count': db,
            'main': dx
        }


    def str_distance(self, a, b, partial=True):
        """
        Check index
        """
        if partial:
            a = a[:len(b)]
            b = b[:len(a)]
        return lev.distance(a, b)


    def is_compatible(self, x):
        """
        x list of index
        Check index/barcode, mismatch
        """
        if isinstance(x, str):
            flag = True
        elif isinstance(x, list):        
            # remove null from x
            x = [i.upper() for i in x if not i.upper() == 'NULL']
            flag = True
            for i, k in combinations(x, 2):
                if self.str_distance(i, k) <= self.mismatch:
                    flag = False
            return flag
        elif isinstance(x, dict):
            return self.is_compatible(list(x.keys()))


    def is_valid_index(self, x):
        """
        x, str, list, dict;
        Check the index sequence is valid: ACGTN
        skip: NULL
        """
        if isinstance(x, str):
            x = x.upper()
            return True if x == 'NULL' or re.match('^[ACGTN]*$', x) else False
        elif isinstance(x, list):
            return all([self.is_valid_index(i) for i in x])
        elif isinstance(x, dict):
            return all([self.is_valid_index(i) for i in list(x.keys())])
        else:
            return False


    def check(self):
        """
        Check the index: index1, index2, barcode
        ACGTN only
        barcode, same width
        """
        ## unique all
        m = [i for i in self.count.values() if i > 1] # index unique
        mn = len(self.sample.values()) == len(set(self.sample.values())) # filename unique

        # width
        self.barcode_width = list(set([len(i) for i in self.barcode]))

        # status
        ss = 'ok' if len(m) == 0 else 'failed'
        sn = 'ok' if mn else 'failed'
        s1 = 'ok' if self.is_valid_index(self.index1) else 'failed'
        s2 = 'ok' if self.is_valid_index(self.index2) else 'failed'
        s3 = 'ok' if self.is_valid_index(self.barcode) else 'failed'
        s4 = 'ok' if self.is_compatible(self.index1) else 'failed'
        s5 = 'ok' if self.is_compatible(self.index2) else 'failed'
        s6 = 'ok' if self.is_compatible(self.barcode) else 'failed'
        s7 = 'ok' if len(self.barcode_width) == 1 else 'failed'

        # report
        msg = '\n'.join([
            'Check index status',
            '{0:.<40}: {1:<10}'.format('number of mismatch', self.mismatch),
            '{0:.<40}: {1:<10}'.format('all index unique', ss),
            '{0:.<40}: {1:<10}'.format('filename unique', sn),
            '{0:.<40}: {1:<10}'.format('index-1 [ACGTN]', s1),
            '{0:.<40}: {1:<10}'.format('index-2 [ACGTN]', s2),
            '{0:.<40}: {1:<10}'.format('barcode [ACGTN]', s3),
            '{0:.<40}: {1:<10}'.format('index-1 mismatch', s4),
            '{0:.<40}: {1:<10}'.format('index-2 mismatch', s5),
            '{0:.<40}: {1:<10}'.format('barcode mismatch', s6),
            '{0:.<40}: {1:<10}'.format('barcode width consistent', s7)
        ])

        log.info(msg)
        
        return all([i == 'ok' for i in [ss, sn, s1, s2, s3, s4, s5, s6, s7]])


class Demx(object):
    """
    P7+P5+inline
    SE or PE
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.init_args()


    def init_args(self):
        """
        Deault arguments
        """
        args_init = {
            'fq1': None,
            'fq2': None,
            'outdir': os.getcwd(),
            'index_csv': None,
            'barcode_in_read': 2,
            'barcode_n_left': 0,
            'barcode_n_right': 0,
            'mismatch': 0,
            'threads': 1,
            'parallel_jobs': 4,
            'overwrite': False
        }
        for k, v in args_init.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        # fastq file
        if self.fq1 is None:
            log.error('fq1, file not exists')
            raise ValueError('fq1, fastq file required')

        # sample: csv file
        if self.index_csv is None:
            log.error('index_csv,  file not exists')
            raise ValueError('index_csv, csv file required, filename,index1,index2,barcode')
        else:
            self.idx = IndexList(input=self.index_csv, mismatch=self.mismatch)
            self.sample = self.idx.sample # index1,index2,barcode=sample
            self.index1 = self.idx.index1
            self.index2 = self.idx.index2
            self.barcode = self.idx.barcode
            self.barcode_width = self.idx.barcode_width[0] #
            self.idx_main = self.idx.main
            if not self.idx.status:
                raise ValueError('Check out the index list file')

        # mismatch
        if not self.mismatch in range(3):
            log.error('mismatch, [0,1,2,3] expected, {} got'.format(self.mismatch))
            self.mismatch = 0


    def mission(self):
        """
        Determine the mission of demx
        [pe/se][p7][p5][barcode]
        0101
        0100
        0001
        0000
        1101
        1100
        1001
        1000
        """
        is_pe = 0 if self.fq2 is None else 1 # PE
        is_p7 = 0
        is_p5 = 0
        is_bc = 0
        if not(all([i == 'NULL' for i in self.index1])):
            is_p7 = 1
        if not(all([i == 'NULL' for i in self.index2])):
            is_p5 = 1
        if not(all([i == 'NULL' for i in self.barcode])):
            is_bc = 1

        code = [is_pe, is_p7, is_p5, is_bc]

        if code == [1, 1, 0, 1]:
            self.run_pe_p7_bc()
        elif code == [1, 1, 0, 0]:
            self.run_pe_p7()
        elif code == [1, 0, 0, 1]:
            self.run_pe_bc()
        elif code == [0, 1, 0, 1]:
            self.run_se_p7_bc()
        elif code == [0, 1, 0, 0]:
            self.run_se_p7()
        elif code == [0, 0, 0, 1]:
            self.run_se_bc()
        else:
            log.error('unknown mission')


    def run_pe_p7_bc(self):
        # step1. index1
        self.dir_d1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1) as r1, xopen(self.fq2) as r2:
            self.index_pe(r1, r2, self.dir_d1)
        # save files in dirs/
        self.wrap_dir(self.dir_d1, mode='index1')
        
        # step2. index2
        # step3. barcode
        bc_fqs = []
        for k1, v1 in self.idx_main.items():
            s_dir1 = os.path.join(self.dir_d1, k1) # 
            for k2, v2 in v1.items():
                for k3, v3 in v2.items():
                    if k3 == 'NULL':
                        continue
                    bc_fq = listfile(s_dir1, "*_1.fq")
                    bc_fqs.extend(bc_fq)

        # multiple jobs
        # with Pool(processes=self.parallel_jobs) as pool:
        #     pool.map(self.run_se_barcode_single, bc_fqs)
        for fq in bc_fqs:
            self.run_pe_barcode_single(fq)

        # wrap dirs
        bc_dirs = list(set([os.path.join(os.path.dirname(i), '_tmp') for i in bc_fqs]))
        for i in bc_dirs:
            self.wrap_dir(i, mode='barcode')

        ## step4. rename files
        self.wrap_pe_file()


    def run_pe_p7(self):
        # step1. index1
        self.dir_d1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1) as r1, xopen(self.fq2) as r2:
            self.index_pe(r1, r2, self.dir_d1)
        # save files in dirs/
        self.wrap_dir(self.dir_d1, mode='index1')

        self.wrap_pe_file()
        
        
    def run_pe_bc(self):
        pass


    def run_se_p7_bc(self):
        """
        p7 index
        barcode
        """
        # step1. index1
        self.dir_d1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1) as r:
            self.index_se(r, self.dir_d1)
        # save files in dirs/
        self.wrap_dir(self.dir_d1, mode='index1')

        # step2. index2
        # step3. barcode
        bc_fqs = []
        for k1, v1 in self.idx_main.items():
            s_dir1 = os.path.join(self.dir_d1, k1) # 
            for k2, v2 in v1.items():
                for k3, v3 in v2.items():
                    if k3 == 'NULL':
                        continue
                    bc_fq = listfile(s_dir1, "*.fq")
                    bc_fqs.extend(bc_fq)


        # multiple jobs
        # with Pool(processes=self.parallel_jobs) as pool:
        #     pool.map(self.run_se_barcode_single, bc_fqs)
        for fq in bc_fqs:
            self.run_se_barcode_single(fq)

        # wrap dirs
        bc_dirs = list(set([os.path.join(os.path.dirname(i), '_tmp') for i in bc_fqs]))
        for i in bc_dirs:
            self.wrap_dir(i, mode='barcode')

        # step4. rename files
        self.wrap_se_file()


    def run_se_p7(self):
        """
        structure: 
        |--_tmp
        |--index1
        |    |--barcode
        |    |    |--file.fq
        """
        # step1. index1
        self.dir_d1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1) as r:
            self.index_se(r, self.dir_d1)
        # save files in dirs/
        self.wrap_dir(self.dir_d1, mode='index1')

        # output
        self.wrap_se_file()
     

    def run_se_bc(self):
        pass


    def wrap_dir(self, path, mode='index1'):
        """
        wrap files in dir, according to dict.keys()
        support PE
        """
        f_all = listfile(path, "*.fq")
        for f in listfile(path, "*.fq"):
            fname = fq_name(f, pe_fix=True)
            if '_' in fname:
                fname = fname.partition('_')[2] #second 
            fn = self.search(fname, mode) # multiples index
            if len(fn) == 1:
                fsub = fn.pop()
            elif len(fn) > 1:
                fsub = 'NULL'
            else:
                fsub = 'undemx'
            # move files to fdir
            fdir = os.path.join(path, fsub)
            check_path(fdir)
            fnew = os.path.join(fdir, os.path.basename(f))
            if os.path.exists(fnew):
                log.info('file exists: {}'.format(fnew))
            else:
                os.rename(f, fnew) # target exists


    def wrap_se_file(self):
        """
        rename sub files, in outdir
        """
        print('!DDDD-1', self.sample)
        undemx_fqs = []
        # level-1: index-1
        dir1 = os.path.join(self.outdir, '_tmp')
        for d1 in listdir(dir1, include_dir=True):
            d1_name = os.path.basename(d1)
            # skip exception
            if not os.path.isdir(d1):
                continue
            print('!AAAA-1', d1)
            # loop over dirs
            d2 = os.path.join(d1, '_tmp')
            if os.path.exists(d2):
                # with barcode
                dir2 = os.path.join(d1, '_tmp')
                for d2 in listdir(dir2, include_dir=True):
                    # skip exception
                    if not os.path.isdir(d2):
                        continue
                    # loop over dirs
                    print('!AAAA-2', d2)
                    d2_name = os.path.basename(d2)
                    if '_' in d2_name:
                        d2_name = d2_name.partition('_')[2] #second 
                    idx = ','.join([d1_name, 'NULL', d2_name])
                    fname = self.sample.get(idx, 'undemx')
                    print('!AAAA-3', idx, fname)
                    fout = os.path.join(self.outdir, fname + '.fq.gz')
                    qlist = listfile(d2, '*.fq')
                    if fname == 'undemx':
                        undemx_fqs.extend(qlist)
                    else:
                        self.fq_merge(fout, qlist)
            else:
                # without barcode
                idx = ','.join(os.path.basename(d1), 'NULL', 'NULL')
                fname = self.sample.get(idx, 'undemx')
                fout = os.path.join(self.out, fname + '.fq.gz')
                qlist = listfile(d1, '*.fq')
                self.fq_merge(fout, qlist)

        # undemx
        fout = os.path.join(self.outdir, 'undemx.fq.gz')
        self.fq_merge(fout, undemx_fqs)


    def wrap_pe_file(self):
        """
        rename sub files, in outdir
        support: PE
        """
        print('!DDDD-1', self.sample)
        undemx_r1_fqs = []
        undemx_r2_fqs = []
        # level-1: index-1
        dir1 = os.path.join(self.outdir, '_tmp')
        for d1 in listdir(dir1, include_dir=True):
            d1_name = fq_name(d1, pe_fix=True)
            # skip exception
            if not os.path.isdir(d1):
                continue
            print('!AAAA-1', d1)
            # loop over dirs
            d2 = os.path.join(d1, '_tmp')
            if os.path.exists(d2):
                # with barcode
                dir2 = os.path.join(d1, '_tmp')
                for d2 in listdir(dir2, include_dir=True):
                    # skip exception
                    if not os.path.isdir(d2):
                        continue
                    # loop over dirs
                    print('!AAAA-2', d2)
                    # d2_name = os.path.basename(d2)
                    d2_name = fq_name(d2, pe_fix=True)
                    if '_' in d2_name:
                        d2_name = d2_name.partition('_')[2] #second 
                    idx = ','.join([d1_name, 'NULL', d2_name])
                    fname = self.sample.get(idx, 'undemx')
                    print('!AAAA-3', idx, fname)
                    r1_fout = os.path.join(self.outdir, fname + '_1.fq.gz')
                    r2_fout = os.path.join(self.outdir, fname + '_2.fq.gz')
                    r1_qlist = listfile(d2, '*_1.fq')
                    r2_qlist = listfile(d2, '*_2.fq')
                    if fname == 'undemx':
                        undemx_r1_fqs.extend(r1_qlist)
                        undemx_r2_fqs.extend(r2_qlist)
                    else:
                        self.fq_merge(r1_fout, r1_qlist)
                        self.fq_merge(r2_fout, r2_qlist)
            else:
                # without barcode
                idx = ','.join([os.path.basename(d1), 'NULL', 'NULL'])
                fname = self.sample.get(idx, 'undemx')
                r1_fout = os.path.join(self.outdir, fname + '_1.fq.gz')
                r2_fout = os.path.join(self.outdir, fname + '_2.fq.gz')
                r1_qlist = listfile(d1, '*_1.fq')
                r2_qlist = listfile(d1, '*_2.fq')
                self.fq_merge(r1_fout, r1_qlist)
                self.fq_merge(r2_fout, r2_qlist)

        # undemx
        r1_fout = os.path.join(self.outdir, 'undemx_1.fq.gz')
        r2_fout = os.path.join(self.outdir, 'undemx_2.fq.gz')
        self.fq_merge(r1_fout, undemx_r1_fqs)
        self.fq_merge(r2_fout, undemx_r2_fqs)


    def fq_merge(self, fout, qlist):
        """
        Combine multiple fastq files into single file
        """
        with gzip.open(fout, 'wb') as w:
            for q in qlist:
                with open(q, 'rb') as r:
                    shutil.copyfileobj(r, w)


    # multiple threads
    def run_se_barcode_single(self, fq):
        prefix = fq_name(fq) #
        with open(fq) as r:
            fname = fq_name(fq, pe_fix=False)
            outdir = os.path.join(os.path.dirname(fq), '_tmp')
            self.barcode_se(r, outdir, prefix)


    # multiple threads
    def run_pe_barcode_single(self, fq1):
        fq2 = fq1.replace('1.fq', '2.fq')
        fname = fq_name(fq1, pe_fix=True)
        outdir = os.path.join(os.path.dirname(fq1), '_tmp')
        with open(fq1) as r1, open(fq2) as r2:
            self.barcode_pe(r1, r2, outdir, fname)


    def search(self, x, mode='index1'):
        """
        Search by index
        support for: wrap dir
        return the sample name
        """
        # df = {
        #     'index1': self.index1,
        #     'index2': self.index2,
        #     'barcode': self.barcode
        # }
        # d = df.get(mode, None)
        # h = [i for i in d if self.str_distance(i, x) <= self.mismatch]
        # return h

        if mode == 'index1':
            h = [i for i in self.index1 if self.str_distance(i, x) <= self.mismatch]
        elif mode == 'barcode':
            h = [i for i in self.barcode if self.str_distance(i, x) <= self.mismatch]
        else:
            h = []

        return h


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


    def str_distance(self, a, b, partial=True):
        """
        a, string 
        b, string
        partial, bool, whether compare the same length
        The distance between a and b
        """
        if partial:
            a = a[:len(b)]
            b = b[:len(a)]
        return lev.distance(a, b)


    # p7 index
    def index_se(self, fh, outdir):
        """
        Demultiplex Illumina Hiseq fastq file
        single index or dual index
        1. [ACGTN]{8}
        2. [ACGTN]{8}+[ACGTN]{8} 
        """
        # read fastq files
        # name, seq, qual, comment        
        check_path(outdir)

        # file handler
        dfh = {}

        # go over fastq file
        for name, seq, qual, comment in self.readfq(fh):
            index = comment.split(':')[-1] # '2:N:0:ATCACGAT+AGATCTCG'
            index = index.partition('+')[0] # for single index only
            f_name = index if isinstance(index, str) else 'not'
            f_file = os.path.join(outdir, f_name + '.fq')

            if not f_name in dfh:
                dfh[f_name] = open(f_file, 'wt')
            dfh[f_name].write('\n'.join([
                '@' + name + ' ' + comment,
                seq, '+', qual]) + '\n')

        # close fh
        for k, v in dfh.items():
            v.close()

        return 


    # p7 index: 
    def index_pe(self, fh1, fh2, outdir):
        """
        Demultiplex Illumina Hiseq fastq file
        single index or dual index
        1. [ACGTN]{8}
        2. [ACGTN]{8}+[ACGTN]{8} 
        """
        # read fastq files
        # name, seq, qual, comment    
        # file-handle for files
        check_path(outdir)
        dfh = {}

        for r1, r2 in zip(self.readfq(fh1), self.readfq(fh2)):
        # for name, seq, qual, comment in readfq(fh):
        # determined by read1
            r1_name, r1_seq, r1_qual, r1_comment = r1
            r2_name, r2_seq, r2_qual, r2_comment = r2
            # check output
            index = r1_comment.split(':')[-1] # '2:N:0:ATCACGAT+AGATCTCG'
            index = index.partition('+')[0] # for single index only
            f_name = index if isinstance(index, str) else 'not'
            f_file1 = os.path.join(outdir, f_name + '_1.fq')
            f_file2 = os.path.join(outdir, f_name + '_2.fq')

            if not f_name in dfh:
                dfh[f_name] = [open(f_file1, 'wt'), open(f_file2, 'wt')]

            # file handler
            fw1, fw2 = dfh.get(f_name, [None, None])
            fw1.write('\n'.join([
                '@' + r1_name + ' ' + r1_comment,
                r1_seq, '+', r1_qual]) + '\n')
            fw2.write('\n'.join([
                '@' + r2_name + ' ' + r2_comment,
                r2_seq, '+', r2_qual]) + '\n')

        # close fh
        for k, v in dfh.items():
            v[0].close() # read1
            v[1].close() # read2

        
    # barcode
    def barcode_se(self, fh, outdir, prefix=None):
        """
        check barcode
        save as files
        """
        # read fastq files
        # name, seq, qual, comment        
        check_path(outdir)

        # file handler
        dfh = {}

        # go over fastq file
        for name, seq, qual, comment in self.readfq(fh):
            s = self.barcode_n_left
            w = self.barcode_width
            index = seq[s:(s+w)]
            if prefix is None:
                f_file = os.path.join(outdir, index + '.fq')
            else:
                f_file = os.path.join(outdir, prefix + '_' + index + '.fq')
            if not index in dfh:
                dfh[index] = open(f_file, 'wt')
            dfh[index].write('\n'.join([
                '@' + name + ' ' + comment,
                seq, '+', qual]) + '\n')

        # close fh
        for k, v in dfh.items():
            v.close()


    # barcode
    def barcode_pe(self, fh1, fh2, outdir, prefix=None):
        """
        check barcode
        save as files
        """
        # read fastq files
        # name, seq, qual, comment        
        check_path(outdir)

        # file handler
        dfh = {}

        # go over fastq file        
        for r1, r2 in zip(self.readfq(fh1), self.readfq(fh2)):
            r1_name, r1_seq, r1_qual, r1_comment = r1
            r2_name, r2_seq, r2_qual, r2_comment = r2
            # barcode in index
            seq = r1_seq if self.barcode_in_read == 1 else r2_seq
            s = self.barcode_n_left
            w = self.barcode_width
            f_name = index = seq[s:(s+w)]
            # output
            if prefix is None:
                f_file1 = os.path.join(outdir, f_name + '_1.fq')
                f_file2 = os.path.join(outdir, f_name + '_2.fq')
            else:
                f_file1 = os.path.join(outdir, prefix + '_' + f_name + '_1.fq')
                f_file2 = os.path.join(outdir, prefix + '_' + f_name + '_2.fq')
            if not f_name in dfh:
                dfh[f_name] = [open(f_file1, 'wt'), open(f_file2, 'wt')]
            # write to file
            fw1, fw2 = dfh.get(f_name, [None, None])
            fw1.write('\n'.join([
                '@' + r1_name + ' ' + r1_comment,
                r1_seq, '+', r1_qual]) + '\n')
            fw2.write('\n'.join([
                '@' + r2_name + ' ' + r2_comment,
                r2_seq, '+', r2_qual]) + '\n')
 
        # close fh
        for k, v in dfh.items():
            v[0].close()
            v[1].close()


    def run(self):
        self.mission()



