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
import shutil
import argparse
from contextlib import ExitStack
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
        self.index1 = list(set([i[0] for i in m if not i == 'NULL'])) #list(self.main.keys()) #
        self.index2 = list(set([i[1] for i in m if not i == 'NULL']))
        self.barcode = list(set([i[2] for i in m if not i == 'NULL']))
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
        self.barcode_width = list(set([len(i) for i in self.barcode if not i == 'NULL']))

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

        # clean outdir
        tmp = listdir(self.outdir, include_dir=True)
        if len(tmp) > 0:
            raise Exception('outdir: empty dir required, files detected in: {}'.format(self.outdir))


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
        log.info('Mission: pe:{} p7:{} p5:{} barcode:{}'.format(is_pe, is_p7, is_p5, is_bc))

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


    def run_se_p7_bc(self):
        """
        p7 index
        barcode
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with gzip.open(self.fq1, 'rt') as r:
            self.index_se(r, outdir1)
        flist1 = self.wrap_dir(outdir1, 'index1') # filenames, not-filenames

        ## parallel-jobs
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.run_se_barcode_single, flist1)
        # for fq in flist1:
        #     flist2 = self.run_se_barcode_single(fq)

        # rename files
        self.wrap_read_count()
        self.wrap_file()


    def run_se_p7(self):
        """
        structure: 
        |--_tmp
        |--index1
        |    |--barcode
        |    |    |--file.fq
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with gzip.open(self.fq1, 'rt') as r:
            self.index_se(r, outdir1)
        flist1 = self.wrap_dir(outdir1, 'index1') # filenames, not-filenames

        # save file
        self.wrap_read_count()
        self.wrap_file()


    def run_se_bc(self):
        pass


    def run_pe_p7_bc(self):
        """
        structure: 
        |--_tmp
        |--index1
        |    |--barcode
        |    |    |--file.fq
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with gzip.open(self.fq1, 'rt') as r1, gzip.open(self.fq2, 'rt') as r2:
            self.index_pe(r1, r2, outdir1)
        # rename files in dirs/
        flist1 = self.wrap_dir(outdir1, mode='index1')
        flist1 = sorted(flist1) # sort, read1/read2

        # step2. index2
        # step3. barcode
        # multiple jobs
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.run_pe_barcode_single, flist1[0::2]) # read1
        # for fq1 in flist1[0::2]: # read1
        #     self.run_pe_barcode_single(fq1)

        ## step4. rename files
        self.wrap_read_count()
        self.wrap_file()


    def run_pe_p7(self):
        """
        structure: 
        |--_tmp
        |--index1
        |    |--barcode
        |    |    |--file.fq
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with gzip.open(self.fq1, 'rt') as r1, gzip.open(self.fq2, 'rt') as r2:
            self.index_pe(r1, r2, outdir1)
        # rename files in dirs/
        flist1 = self.wrap_dir(outdir1, mode='index1')

        # save files
        self.wrap_read_count()
        self.wrap_file()
        

    def run_pe_bc(self):
        pass


    def fq_merge(self, fout, qlist):
        """
        Compress, multiple fastq files into single file
        """
        with gzip.open(fout, 'wb') as w:
            for q in qlist:
                with open(q, 'rb') as r:
                    shutil.copyfileobj(r, w)


    def wrap_dir(self, x, mode):
        """
        rename fq files in dir
        """
        dirs = []
        if mode == 'index1':
            for f in listfile(x, "*.fq"):
                fname = fq_name(f, pe_fix=True)
                ix = fname + ',NULL,NULL'
                f_new = os.path.join(x, self.sample.get(ix, fname))
                # check suffix
                f_new += '_1.fq' if f.endswith('_1.fq') else ('_2.fq' if f.endswith('_2.fq') else '.fq')
                if not f == f_new:
                    if os.path.exists(f_new):
                        log.warning('file exists: {}'.format(f_new))
                    else:
                        os.rename(f, f_new)
                else:
                    dirs.append(f)
        elif mode == 'barcode':
            for f in listfile(x, "*.fq"):
                fname = fq_name(f, pe_fix=True)
                ix =  fq_name(x) + ',NULL,' + fname
                f_new = os.path.join(x, self.sample.get(ix, fname))
                f_new += '_1.fq' if f.endswith('_1.fq') else ('_2.fq' if f.endswith('_2.fq') else '.fq')
                if not f == f_new:
                    if os.path.exists(f_new):
                        log.warning('file exists: {}'.format(f_new))
                    else:
                        os.rename(f, f_new)
                else:
                    dirs.append(f)
        else:
            pass

        return dirs


    def wrap_file(self):
        """
        move "_tmp/" files to "outdir", 
        combine undemx file
        """
        f_undemx = []
        f_hits = []
        # level-1: index-1
        dv1 = os.path.join(self.outdir, '_tmp')
        for d1 in listdir(dv1, include_dir=True):
            d1_name = fq_name(d1, pe_fix=True)
            # for exists files
            if d1_name in self.sample.values():
                f_hits.append(d1)
            elif d1_name == 'undemx':
                f_undemx.append(d1)
            elif os.path.isdir(d1): # barcode level
                # level-2: barcode
                for d2 in listfile(d1, "*.fq"):
                    d2_name = fq_name(d2, pe_fix=True)
                    if d2_name in self.sample.values():
                        f_hits.append(d2)
                    elif d2_name == "undemx":
                        f_undemx.append(d2)
                    else:
                        pass
            else:
                pass
        # save files
        for f in f_hits:
            f_out = os.path.join(self.outdir, os.path.basename(f) + '.gz')
            log.info('Saving file: {}'.format(f_out))
            gzip_cmd(f, f_out, decompress=False, rm=True)
        
        # save undemx
        f_undemx = sorted(f_undemx)
        if f_undemx[0].endswith('_1.fq'):
            # PE mode
            f_r1_undemx_out = os.path.join(self.outdir, 'undemx_1.fq.gz')
            f_r2_undemx_out = os.path.join(self.outdir, 'undemx_2.fq.gz')
            self.fq_merge(f_r1_undemx_out, f_undemx[0::2]) # read1
            self.fq_merge(f_r2_undemx_out, f_undemx[1::2]) # read2
        else:
            # SE mode
            f_undemx_out = os.path.join(self.outdir, 'undemx.fq.gz')
            self.fq_merge(f_undemx_out, f_undemx)

        # remove temp files "outdir/_tmp"


    def wrap_read_count(self):
        """
        parse the "read_number.json" file
        save read counts
        """        
        ##----------------------------------##
        # for index1
        # main dir
        f1 = os.path.join(self.outdir, '_tmp', 'read_number.json')
        if os.path.exists(f1):
            da = Json(f1).reader()
        else:
            log.error('file not exists: {}'.format(f1))

        # for barcode
        db = {} # subdir
        for fv1 in listdir(os.path.join(self.outdir, '_tmp'), include_dir=True):
            fv1_name = fq_name(fv1)
            if not os.path.isdir(fv1):
                continue
            # subdir
            fv2 = os.path.join(fv1, 'read_number.json')
            if os.path.exists(fv2):
                db[fv1_name] = Json(fv2).reader()
            else:
                log.error('file not exists: {}'.format(fv2))

        ##----------------------------------##
        dd = {} # all report
        # assign read number for each file
        n_undemx = 0
        dx = {}
        for k, v in da.items():
            ix = k + ',NULL,NULL'
            if ix in self.sample:
                dx[self.sample.get(ix, k)] = v
            elif k == 'undemx':
                n_undemx += v
            else:
                pass

        # for subdir
        for k1, v1 in db.items():
            dd[k1] = {}
            for k2, v2 in v1.items():
                ix = k1 + ',NULL,' + k2
                if ix in self.sample:
                    smp = self.sample.get(ix, k2)
                    dx[smp] = v2
                    dd[k1][smp] = v2
                elif k2 == 'undemx':
                    n_undemx += v2
                    dd[k1]['undemx'] = v2
                else:
                    pass
        # for undemx
        dx['undemx'] = n_undemx

        # for all
        dd['main'] = {}
        for k, v in sorted(dx.items(), key=lambda x: x[0], reverse=False):
            dd['main'][k] = v

        ##----------------------------------##
        # save to final report
        self.demx_report = dd
        f_out = os.path.join(self.outdir, 'demx_report.json')
        Json(dd).writer(f_out)


    def run_se_barcode_single(self, fq):
        """
        Demultiplex barcode: read again
        """
        fname = fq_name(fq, pe_fix=True)
        if fname == 'undemx':
            return None # skip
        outdir = os.path.join(os.path.dirname(fq), fname)
        with xopen(fq) as r:
            self.barcode_se(r, outdir, index1=fname)
        return self.wrap_dir(outdir, 'barcode')


    def run_pe_barcode_single(self, fq1):
        """
        Demultiplex barcode: read again
        """
        fname = fq_name(fq1, pe_fix=True)
        fq2 = fq1.replace('_1.fq', '_2.fq')
        if fname == 'undemx':
            return None # skip
        outdir = os.path.join(os.path.dirname(fq1), fname)
        with open(fq1, 'rt') as r1, open(fq2, 'rt') as r2:
            self.barcode_pe(r1, r2, outdir, index1=fname)
        return self.wrap_dir(outdir, 'barcode')


    def search(self, x, mode='index1', index1=None):
        """
        Search by index
        support for: wrap dir
        return the sample name
        """
        if mode == 'index1':
            h = [i for i in self.index1 if self.str_distance(i, x) <= self.mismatch]
        elif mode == 'barcode':
            h = [i for i in self.get_barcode(index1) if self.str_distance(i, x) <= self.mismatch]
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


    def get_barcode(self, index1=None):
        """
        Return the barcode list
        according to the index1
        """
        if index1 in self.idx_main:
            return list(self.idx_main[index1]['NULL'].keys())
        else:
            return self.barcode


    # p7 index
    def index_se(self, fh, outdir):
        """
        Demultiplex Illumina Hiseq fastq file
        single index or dual index
        1. [ACGTN]{8}
        2. [ACGTN]{8}+[ACGTN]{8} 
        """
        check_path(outdir)

        # all availabe output files
        fn = {} # fq number
        fnames = self.index1
        fnames.append('undemx')
        fouts = [os.path.join(outdir, i + '.fq') for i in fnames]
        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in fouts]
            n = 0 # counter
            for name, seq, qual, comment in self.readfq(fh):
                n += 1
                if n%1000000 == 0 :
                    log.info('Processed reads: {}'.format(n))

                fq = '\n'.join(['@' + name + ' ' + comment, seq, '+', qual])
                index = comment.split(':')[-1] # '2:N:0:ATCACGAT+AGATCTCG'
                # index = index.partition('+')[0] # for single index only
                index = index.split('+')[0]
                fhit = self.search(index, 'index1')
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw = fws[fnames.index(fname)]
                fw.write(fq + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        fn_json = os.path.join(outdir, 'read_number.json')
        Json(fn).writer(fn_json)        

        return fouts


    # p7 index: 
    def index_pe(self, fh1, fh2, outdir):
        """
        Demultiplex Illumina Hiseq fastq file
        single index or dual index
        1. [ACGTN]{8}
        2. [ACGTN]{8}+[ACGTN]{8} 
        """
        check_path(outdir)

        # all availabe output files
        fn = {} # fq read number
        fnames = self.index1
        fnames.append('undemx')
        r1_fouts = [os.path.join(outdir, i + '_1.fq') for i in fnames]
        r2_fouts = [os.path.join(outdir, i + '_2.fq') for i in fnames]
        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in r1_fouts + r2_fouts]
            n = 0 # counter
            for r1, r2 in zip(self.readfq(fh1), self.readfq(fh2)):
                n += 1
                if n%1000000 == 0:
                    log.info('Processed reads: {}'.format(n))

                r1_name, r1_seq, r1_qual, r1_comment = r1
                r2_name, r2_seq, r2_qual, r2_comment = r2
                fq1 = '\n'.join(['@' + r1_name + ' ' + r1_comment, r1_seq, '+', r1_qual])
                fq2 = '\n'.join(['@' + r2_name + ' ' + r2_comment, r2_seq, '+', r2_qual])
                # check output
                index = r1_comment.split(':')[-1] # '2:N:0:ATCACGAT+AGATCTCG'
                # index = index.partition('+')[0] # for single index only
                index = index.split('+')[0]
                fhit = self.search(index, 'index1')
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw1 = fws[fnames.index(fname)]
                fw2 = fws[fnames.index(fname) + len(fnames)]
                fw1.write(fq1 + '\n')
                fw2.write(fq2 + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        fn_json = os.path.join(outdir, 'read_number.json')
        Json(fn).writer(fn_json)


    # barcode
    def barcode_se(self, fh, outdir, index1=None):
        """
        check barcode
        save as files
        """   
        check_path(outdir)

        # add undemx
        fn = {} # fq read number
        bc = self.get_barcode(index1)
        bc.append('undemx')
        fouts = [os.path.join(outdir, i + '.fq') for i in bc]
        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in fouts]
            n = 0 # counter
            for name, seq, qual, comment in self.readfq(fh):
                n += 1
                if n%1000000 == 0 :
                    log.info('Processed reads: {}'.format(n))
                fq = '\n'.join(['@' + name + ' ' + comment, seq, '+', qual])
                # fq = '@' + name + ' ' + comment + '\n' + seq + '\n+\n' + qual
                s = self.barcode_n_left
                w = self.barcode_width
                index = seq[s:(s+w)]
                fhit = self.search(index, 'barcode', index1) # search names
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw = fws[bc.index(fname)]
                fw.write(fq + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        fn_json = os.path.join(outdir, 'read_number.json')
        Json(fn).writer(fn_json)


    # barcode
    def barcode_pe(self, fh1, fh2, outdir, index1=None):
        """
        check barcode
        save as files
        """ 
        check_path(outdir)

        # add undemx
        fn = {} # fq read number 
        bc = self.get_barcode(index1)
        bc.append('undemx')
        r1_fouts = [os.path.join(outdir, i + '_1.fq') for i in bc]
        r2_fouts = [os.path.join(outdir, i + '_2.fq') for i in bc]
        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in r1_fouts + r2_fouts]
            n = 0 # counter
            for r1, r2 in zip(self.readfq(fh1), self.readfq(fh2)):
                n += 1
                if n%1000000 == 0:
                    log.info('Processed reads: {}'.format(n))

                r1_name, r1_seq, r1_qual, r1_comment = r1
                r2_name, r2_seq, r2_qual, r2_comment = r2
                fq1 = '\n'.join(['@' + r1_name + ' ' + r1_comment, r1_seq, '+', r1_qual])
                fq2 = '\n'.join(['@' + r2_name + ' ' + r2_comment, r2_seq, '+', r2_qual])

                # check barcode
                s = self.barcode_n_left
                w = self.barcode_width
                seq = r1_seq if self.barcode_in_read == 1 else r2_seq # which read1/2
                index = seq[s:(s+w)]
                fhit = self.search(index, 'barcode', index1) # search names
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw1 = fws[bc.index(fname)]
                fw2 = fws[bc.index(fname) + len(bc)]
                fw1.write(fq1 + '\n')
                fw2.write(fq2 + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        fn_json = os.path.join(outdir, 'read_number.json')
        Json(fn).writer(fn_json)


    def run(self):
        self.mission()
        # report
        print('RT: {:-^64}'.format('Demx Report'))
        print('RT: {0:>50} {1:>10}  {2:7}'.format('filename', 'count', 'percent'))
        total = sum(self.demx_report['main'].values())
        for k, v in self.demx_report['main'].items():
            print('RT: {0:>50} {1:>10} {2:>7.2f}%'.format(k, v, v / total * 100))

        print('RT: {0:>50} {1:>10}  {2:7}'.format('sum', total, '100.00%'))







    # def wrap_se_dir(self, x, mode):
    #     """
    #     rename fq files in dir
    #     """
    #     dirs = []
    #     if mode == 'index1':
    #         for f in listfile(x, "*.fq"):
    #             fname = fq_name(f, pe_fix=True)
    #             ix = fname + ',NULL,NULL'
    #             f_new = os.path.join(x, self.sample.get(ix, fname) + '.fq')
    #             if not f == f_new:
    #                 if os.path.exists(f_new):
    #                     log.warning('file exists: {}'.format(f_new))
    #                 else:
    #                     os.rename(f, f_new)
    #             else:
    #                 dirs.append(f)
    #     elif mode == 'barcode':
    #         for f in listfile(x, "*.fq"):
    #             fname = fq_name(f, pe_fix=True)
    #             ix =  fq_name(x) + ',NULL,' + fname
    #             f_new = os.path.join(x, self.sample.get(ix, fname) + '.fq')
    #             if not f == f_new:
    #                 if os.path.exists(f_new):
    #                     log.warning('file exists: {}'.format(f_new))
    #                 else:
    #                     os.rename(f, f_new)
    #             else:
    #                 dirs.append(f)
    #     else:
    #         pass

    #     return dirs


    # def wrap_se_file(self):
    #     """
    #     move "_tmp/" files to "outdir", 
    #     combine undemx file
    #     """
    #     f_undemx = []
    #     f_hits = []
    #     # level-1: index-1
    #     dv1 = os.path.join(self.outdir, '_tmp')
    #     for d1 in listdir(dv1, include_dir=True):
    #         d1_name = fq_name(d1)
    #         # for exists files
    #         if d1_name in self.sample.values():
    #             f_hits.append(d1)
    #         elif d1_name == 'undemx':
    #             f_undemx.append(d1)
    #         elif os.path.isdir(d1): # barcode level
    #             # level-2: barcode
    #             for d2 in listfile(d1, "*.fq"):
    #                 d2_name = fq_name(d2)
    #                 if d2_name in self.sample.values():
    #                     f_hits.append(d2)
    #                 elif d2_name == "undemx":
    #                     f_undemx.append(d2)
    #                 else:
    #                     pass
    #         else:
    #             pass
    #     # save files
    #     for f in f_hits:
    #         f_out = os.path.join(self.outdir, os.path.basename(f) + '.gz')
    #         log.info('Saving file: {}'.format(f_out))
    #         gzip_cmd(f, f_out, decompress=False, rm=True)
        
    #     # save undemx
    #     f_undemx_out = os.path.join(self.outdir, 'undemx.fq.gz')
    #     self.fq_merge(f_undemx_out, f_undemx)

    #     # remove temp files "outdir/_tmp"

