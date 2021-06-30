
# -*- coding: utf-8 -*-

"""
Functions for sequence, fastx

"""

import os
import sys
import re
import logging
import shutil
import numpy as np
import pandas as pd
from xopen import xopen
import pyfastx
import collections # Fastx().collapse()

# hiseq.utils.utils
# circular import
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')



# hiseq.utils.file.check_path()
# circular import
def check_path(x, **kwargs):
    """Check if x is path
    Parameters
    ----------
    x : str
        Path to a file
    
    Keyword Parameters
    ------------------
    show_error : bool
        Show the error messages
        
    show_log : bool
        Show the log messages 
    
    create_dirs : bool
        Create the dirs
    """
    show_error = kwargs.get('show_error', False)
    show_log = kwargs.get('show_log', False)
    create_dirs = kwargs.get('create_dirs', True)
    if isinstance(x, str):
        out = False
        if os.path.isdir(x):
            out = True
        elif os.path.isfile(x):
            if show_error:
                log.error('file exists, not a directory: {}'.format(x))
        else:
            if create_dirs:
                try:
                    os.makedirs(x)
                    out = True
                except:
                    if show_error:
                        log.error('`os.makedirs` failed: {}'.format(x))
        # show log
        flag = 'ok' if out else 'failed'
        if show_log is True:
            log.info('{:<6s} : {}'.format(flag, x))
    elif isinstance(x, list):
        out = all([check_path(i, **kwargs) for i in x])
    else:
        if show_error:
            log.error('x expect str or list, got {}'.format(type(x).__name__))
    return out



# hiseq.utils.utils
# circular import
def update_obj(obj, d, force=True, remove=False):
    """Update the object, by dict
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)
    return obj



class Fastx(object):
    """
    Collection of tools to manipulate fastx file
    1. trimmer: cutadapt [trimmomatic, ...]
    2. collapse: fastx_collapse [fastx_toolkit]
    3. collapse: seqkit rmdup [faster, discard read numbers]
    4. fq2fa: fastq_to_fasta [fastx_toolkit]
    5. fa2fq: [fasta_to_fastq.pl]
    6. revcomp: [seqkit seq]
    7. sample: [seqkit sample -n]
    ...    
    """
    def __init__(self, input, **kwargs):
        """
        read fastq/a file
        """
        self = update_obj(self, kwargs, force=True)
        self.input = input
        self.format = self.fx_type(input)
        self.cat = 'zcat' if input.endswith('.gz') else 'cat'


    def is_not_empty(self, x):
        """Check the file is empty or not
        plain text: 0
        gzip.txt: 20
        """
        return os.stat(x).st_size > 20


    def is_fastq(self):
        return self.fx_type(self.input) == 'fastq'


    def is_fasta(self):
        return self.fx_type(self.input) == 'fasta'


    def file_type(self, fn, top_n=1000):
        """
        Check the file type by top 10000 rows:
        identify @ for fastq, > for fasta, * unknown
        """
        assert isinstance(fn, str)

        if self.is_not_empty(fn):
            d = {}
            counter = 0
            with xopen(fn) as fh:
                for line in fh:
                    counter += 1
                    if counter > top_n:
                        break
                    elif counter % 4 == 1: # for 1st line of fastq; and fa
                        base = line[0] # the first character
                        if base.lower() in 'acgtn': # sequence line in fasta
                            continue
                        d[base] = d.get(base, 0) + 1
                    else:
                        continue

            ## percentage
            x = sorted(d.items(), key=lambda kv:kv[1], reverse=True)
            ## the top1 character
            x_top1 = x[0][0]
            x_top1_pct = x[0][1] / sum(d.values())

            ## check
            if x_top1 == '@':
                fx_type = 'fastq'
            elif x_top1 == '>':
                fx_type = 'fasta'
            else:
                fx_type = None

            ## if top1_pct < 90%
            if x_top1_pct < 0.9:
                fx_type = None
        else:
            fx_type = None

        return fx_type


    def file_ext(self, fn):
        """
        Check the file type by extension: fa/fq/fasta/fastq
        gzip supported
        """
        fname = os.path.basename(fn)
        fname = fname.lower()
        if fn.endswith('gz'):
            fname = os.path.splitext(fname)[0]

        fext = os.path.splitext(fname)[1]

        if fext.lower() in ['.fa', '.fasta']:
            fx_type = 'fasta'
        elif fext.lower() in ['.fq', '.fastq']:
            fx_type = 'fastq'
        else:
            fx_type = None

        return fx_type


    def fx_type(self, fn):
        fx1 = self.file_ext(fn) # extension
        if os.path.exists(fn):
            fx2 = self.file_type(fn) # content
            if fx1 is None or fx2 is None:
                fx_out = fx2 if fx1 is None else fx1
            else:
                if fx1 == fx2:
                    fx_out = fx1
                else:
                    raise Exception('error, filename {} and content {} not \
                        match: {}'.format(fx1, fx2, fn))
        else:
            fx_out = fx1

        return fx_out


    def revcomp(self, out):
        """
        Rev comp the fastx file
        fx_reader()
        """
        with xopen(self.input) as r, xopen(out, 'wt') as w:
            for name, seq, qual, comment in self.readfq(r):
                base_from = 'ACGTNacgtn'
                base_to = 'TGCANtgcan'
                tab = str.maketrans(base_from, base_to)
                seq = seq.translate(tab)[::-1]
                name = name + ' ' + comment
                if qual is None:
                    w.write('\n'.join('>'+name, seq))
                else:
                    qual = qual[::-1]
                    w.write('\n'.join('@'+name, seq, '+', qual))
        

    def sample_random(self, out, n=1000, p=0.01):
        """
        Extract subset of fastx file using seqkit
        fx_reader()
        """
        assert isinstance(n, int)
        assert isinstance(p, float)
        outdir = os.path.dirname(out)
        check_path(outdir, create_dirs=True)

        ## warning
        if n > 1000000:
            log.warning('too-big n={}, change to 1000000'.format(n))
            n = 1000000
        if p > 0.1:
            p = 0.1
            
        ## cmd
        cmd = ' '.join([
            '{}'.format(shutil.which('seqkit')),
            'sample -n {}'.format(n),
            '-o {}'.format(out),
            '{}'.format(self.input)
            ])

        try:
            run_shell_cmd(cmd)
        except:
            log.error('sample_random() faied, {}'.format(self.input))


    def sample(self, out, n=1000):
        """
        Create a subsample of input fastq files, default: 1M reads
        Run the whole process for demostration
        """
        outdir = os.path.dirname(out)
        check_path(outdir, create_dirs=True)
        
        if os.path.exists(out):
            log.info('file eixsts, {}'.format(out))
        else:
            i = 0
            with xopen(self.input, 'rt') as r, xopen(out, 'wt') as w:
                for name, seq, qual, comment in self.readfq(r):
                    i += 1
                    if i > n:
                        break
                    if isinstance(comment, str):
                        if len(comment) > 0:
                            name += ' '+comment
                    if self.format == 'fastq':
                        fx = '@{}\n{}\n+\n{}'.format(name, seq, qual)
                    elif self.format == 'fasta':
                        fx = '>{}\n{}'.format(name, seq)
                    else:
                        continue
                    w.write(fx+'\n')


    def fa2fq(self, out):
        """
        Convert fasta to fastq
        quality='J' Phred = 33
        """
        if not self.format == 'fasta':
            raise Exception('fasta file expected, input: {}'.format(self.input))

        with xopen(self.input) as r, xopen(out, 'wt') as w:
            for name, seq, qual, comment in self.readfq(r):
                qual = 'J' * len(seq) # Phred33, 41
                if isinstance(comment, str):
                    if len(comment) > 0:
                        name += ' '+comment
                w.write('\n'.join('@'+name, seq, '+', qual) + '\n')


    def readfa(self, fh):
        """
        Read fasta file,
        return [name, seq]
        """
        name = seq = ''
        for line in fh:
            if line.startswith('>'):
                head = line.strip().split()[0] # the first item
                head = head.replace('>', '')
                if len(seq) > 0:
                    yield [name, seq]
                name = head
                seq = ''
                continue
            seq += line.strip()

        # last one
        if len(seq) > 0:
            yield [name, seq]


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


    def count(self, x):
        """
        count the number of lines
        source: by Michael Bacon on StackOverflow forum: 
        url: https://stackoverflow.com/a/27518377/2530783.
        """
        def _make_gen(fh):
            block = fh(1024*1024)
            while block:
                yield block
                block = fh(1024*1024)

        with xopen(x, 'rb') as fh:
            return sum(buf.count(b'\n') for buf in _make_gen(fh.read))


    def fq_counter(self, x):
        """
        Count fastq records
        N = (total lines) / 4
        """
        return int(self.count(x) / 4)


    def fa_counter(self, x):
        """
        Count fasta records
        N = sum('>')
        """
        def _make_gen(fh):
            block = fh(1024*1024)
            while block:
                yield block
                block = fh(1024*1024)

        with xopen(x, 'rb') as fh:
            return sum(buf.count(b'\n>') for buf in _make_gen(fh.read))


    def number_of_seq(self):
        """
        Number of sequences
        fa, fq
        """
        return self.fq_counter(self.input) if self.format == 'fastq' else self.fa_counter(self.input)


    def region_to_pos(self, region=None):
        """
        Extract substring by region: 
        Convert 1-indexed region to 0-indexed python style
        The region cloud be: 1-indexed
        1:20,   the first 20 bases
        -20:-1, the last 20 bases
        1:-1,   the full length
        """
        if region is None:
            region = '1:-1' # the full length
        # convert "7,-7" to "7:-7"
        region = region.replace(',', ":")
        p = re.compile('^(-?\d+):(-?\d+)$')
        m = p.search(region)
        if m:
            pass
        else:
            log.error('unknown region format, expect: 1:-1, got: {}'.format(region))
        # 0-indexed
        start = int(m.group(1))
        end = int(m.group(2))
        if start > 0:
            start = start - 1
        if end < 0:
            end = end + 1
        return (start, end)


    def sub_string(self, s, start, end):
        """
        Extract substring, by start, end (0-index)
        convert to python style
        """
        if isinstance(s, str):
            if start == 0 and end == 0:
                pass
            elif start == 0:
                s = s[:end]
            elif end == 0:
                s = s[start:]
            else:
                # exception: start:-end > length
                if start - end > len(s):
                    s = None
                else:
                    s = s[start:end]
        else:
            s = None
        return s


    def subseq(self, out, region=None):
        """
        Get subseq by region

        The region cloud be: 1-indexed
        1:20,   the first 20 bases
        -20:-1, the last 20 bases
        1:-1,   the full length
        """
        start, end = self.region_to_pos(region)
        with xopen(self.input) as r, xopen(out, 'wt') as w:
            for name, seq, qual, comment in self.readfq(r):
                name = name.strip()
                seq = self.sub_string(seq, start, end)
                qual = self.sub_string(qual, start, end)
                # update name
                # !!!! error !!!!
                # STAR 2.5.2a, 
                # does not allow 'white spaces' at the end of name-line of fastq
                # throw error:
                # ReadAlignChunk_processChunks.cpp:115:processChunks EXITING because of FATAL ERROR in input reads: unknown file format: the read ID should start with @ or >
                if isinstance(comment, str):
                    if len(comment) > 0:
                        name += ' '+comment #fix comment=None
                # filter length
                if isinstance(seq, str):
                    if len(seq) < self.len_min:
                        continue
                    if isinstance(qual, str):
                        out = '@{}\n{}\n+\n{}'.format(name, seq, qual)
                    else:
                        out = '>{}\n{}'.format(name, seq)
                    w.write(out+'\n')
                else:
                    continue


    def subseq_pe(self, input2, out1, out2, region=None):
        """
        Get subseq by region, for paired end reads; filter by length

        The region cloud be: 1-indexed
        1:20,   0:20 , the first 20 bases
        -20:-1, -21: , the last 20 bases
        1:-1,   0:   , the full length
        """
        start, end = self.region_to_pos(region)

        with xopen(self.input) as r1, xopen(input2) as r2, \
            xopen(out1, 'wt') as w1, xopen(out2, 'wt') as w2:
            for read1, read2 in zip(self.readfq(r1), self.readfq(r2)):
                name1, seq1, qual1, comment1 = read1
                name2, seq2, qual2, comment2 = read2
                seq1 = self.sub_string(seq1, start, end)
                seq2 = self.sub_string(seq2, start, end)
                if isinstance(comment1, str):
                    if len(comment1) > 0:
                        name1 = name1 + ' ' + comment1 #fix comment=None
                if isinstance(comment2, str):
                    if len(comment2) > 0:
                        name2 = name2 + ' ' + comment2 #fix comment=None
#                 name1 = name1 + ' ' + comment1
#                 name2 = name2 + ' ' + comment2
                # specific length
                if len(seq1) < self.len_min or len(seq2) < self.len_min:
                    continue
                # write
                if qual1 is None: # fa
                    w1.write('\n'.join(['>' + name1, seq1]) + '\n')
                    w2.write('\n'.join(['>' + name2, seq2]) + '\n')
                else:
                    qual1 = self.sub_string(qual1, start, end)
                    qual2 = self.sub_string(qual2, start, end)
                    w1.write('\n'.join(['@' + name1, seq1, '+', qual1]) + '\n')
                    w2.write('\n'.join(['@' + name2, seq2, '+', qual2]) + '\n')


    # Deprecated: (see: subseq)
    def cut(self, out, len_min=15, **kwargs):
        """
        Cut bases from either ends of fasta/q

        7, cut 7-nt from right of sequence (3')
        -5, cut 5-nt from left of sequence (5')
        
        Cut to specific length, from right/left
        

        len_min:
        cut: [7, -5, '7,-5']
        cut_to_length: [30, -25]
        discard_tooshort: [True, False]

        """
        # arguments
        args = kwargs

        cut = args.get('cut', 0)
        cut_to_length = args.get('cut_to_length', 0) # defualt: skip
        discard_tooshort = args.get('discard_tooshort', True)

        # subseq = seq[start:end]
        def cut_sub(x):
            if isinstance(cut, int):
                return x[cut:] if cut > 0 else x[:cut]
            elif isinstance(cut, str):
                if re.match('^\d+,-\d+$', cut):
                    s, e = cut.split(',', 1)
                    s = eval(s)
                    e = eval(e)
                    return x[s:e]
                else:
                    raise Exception('unknown format for cut={}'.format(cut))

        # cut to length
        def cut_sub2(x):
            if isinstance(cut_to_length, int):
                if abs(cut_to_length) < len(x):
                    cut_n = len(x) - abs(cut_to_length)
                    return x[cut_n:] if cut_to_length > 0 else x[:-cut_n]
                else:
                    return x
            else:
                raise Exception('unknown format, cut_to_length={}'.format(cut_to_length))
            # if isinstance(cut_to_length, int):
            #     if cut_to_length > 0: # cut from 3' end
            #         x2 = x[:cut_to_length]
            #     else: # cut from 5' end
            #         n = len(x) - abs(cut_to_length)
            #         if n < 0:
            #             n = 0
            #         x2 = x[n:]
            #     return x2
            # else:
            #     log.error('unknown x, expect int, got {}'.format(cut_to_length))
            #     return x


        # merge two funcs
        def cut_cut(x):
            # cut
            x_cut = cut_sub(x)

            # cut to length
            if not cut_to_length == 0:
                x_cut = cut_sub2(x_cut)

            return x_cut

        with xopen(self.input) as r, xopen(out, 'wt') as w:
            for name, seq, qual, comment in self.readfq(r):
                seq = cut_cut(seq)
                if len(seq) < len_min and discard_tooshort:
                    continue # skip
                # write
                name = name + ' ' + comment
                if qual is None: # fasta
                    w.write('\n'.join(['>'+name, seq]) + '\n')
                else:
                    qual = cut_cut(qual)
                    w.write('\n'.join(['@'+name, seq, '+', qual]) + '\n')


    # Deprecated: (see: subseq)
    def cut_pe(self, input2, out1, out2, len_min=15, **kwargs):
        """
        Cut bases from either ends of fasta/q

        7, cut 7-nt from right of sequence (3')
        -5, cut 5-nt from left of sequence (5')
        
        Cut to specific length, from right/left
        

        len_min:
        cut: [7, -5, '7,-5']
        cut_to_length: [30, -25]
        discard_tooshort: [True, False]

        """
        # arguments
        args = kwargs

        cut = args.get('cut', 0)
        cut_to_length = args.get('cut_to_length', 0) # defualt: skip
        discard_tooshort = args.get('discard_tooshort', True)

        # subseq = seq[start:end]
        def cut_sub(x):
            if isinstance(cut, int):
                return x[cut:] if cut > 0 else x[:cut]
            elif isinstance(cut, str):
                if re.match('^\d+,-\d+$', cut):
                    s, e = cut.split(',', 1)
                    s = eval(s)
                    e = eval(e)
                    return x[s:e]
                else:
                    raise Exception('unknown format for cut={}'.format(cut))

        # cut to length
        def cut_sub2(x):
            if isinstance(cut_to_length, int):
                if abs(cut_to_length) < len(x):
                    cut_n = len(x) - abs(cut_to_length)
                    return x[cut_n:] if cut_to_length > 0 else x[:-cut_n]
                else:
                    return x
            else:
                raise Exception('unknown format, cut_to_length={}'.format(cut_to_length))

        # merge two funcs
        def cut_cut(x):
            # cut
            x_cut = cut_sub(x)

            # cut to length
            if not cut_to_length == 0:
                x_cut = cut_sub2(x_cut)

            return x_cut

        with xopen(self.input) as r1, xopen(input2) as r2, \
            xopen(out1, 'wt') as w1, xopen(out2, 'wt') as w2:
            for read1, read2 in zip(self.readfq(r1), self.readfq(r2)):
                name1, seq1, qual1 = read1
                name2, seq2, qual2 = read2
                seq1_cut = cut_cut(seq1)
                seq2_cut = cut_cut(seq2)
                if discard_tooshort:
                    if len(seq1_cut) < len_min or len(seq2_cut) < len_min:
                        continue # skip pair reads
                # write
                if qual1 is None: # fa
                    w1.write('\n'.join(['>' + name1, seq1]) + '\n')
                    w2.write('\n'.join(['>' + name2, seq2]) + '\n')
                else:
                    qual1_cut = cut_cut(qual1)
                    qual2_cut = cut_cut(qual2)
                    w1.write('\n'.join(['@' + name1, seq1_cut, '+', qual1_cut]) + '\n')
                    w2.write('\n'.join(['@' + name2, seq2_cut, '+', qual2_cut]) + '\n')


    def collapse(self, out, fq_out=False):
        """
        Collapse fastx file, remove PCR duplicates
        sort by counts
        """

        d = {}
        with xopen(self.input) as r:
            for _, seq, _ in self.readfq(r):
                d[seq] = d.get(seq, 0) + 1

        # sort by value
        tmp = sorted(d.items(), key=lambda kv: kv[1], reverse=True)
        dd = collections.OrderedDict(tmp)

        # save to file
        n = 0
        with xopen(out, 'wt') as w:

            for key, value in dd.items():
                n += 1
                name = str(n) + '-' + str(value)
                if fq_out:
                    qual = 'J' * len(key)
                    w.write('\n'.join(['@'+name, key, '+', qual]) + '\n')
                else:
                    w.write('\n'.join(['>'+name, key]) + '\n')


    def detect_adapter(self):
        """
        Guess adapters, sampling the first 1000000 records
        TruSeq    AGATCGGAAGAGC
        Nextera   CTGTCTCTTATACACATCT
        smallRNA  TGGAATTCTCGG

        to-do
        specific type of adapters
        """
        ad = {
            'truseq': 'AGATCGGAAGAGC',
            'nextera': 'CTGTCTCTTATA',
            'smallrna': 'TGGAATTCTCGG'
        }

        # count 
        d = {}
        n_max = 1000000
        n = 0
        with xopen(self.input) as r:
            for _, seq, _, _ in self.readfq(r):
                n += 1
                if n > n_max:
                    break
                # check
                if ad['truseq'] in seq:
                    d['truseq'] = d.get('truseq', 0) + 1
                elif ad['nextera'] in seq:
                    d['nextera'] = d.get('nextera', 0) + 1
                elif ad['smallrna'] in seq:
                    d['smallrna'] = d.get('smallrna', 0) + 1
                else:
                    continue

        # summary
        msg = '\n'.join([
            '{}\t{}\t{}\t{}\t{}'.format('Type', 'sequence', 'total', 'count', 'percent'),
            '{}\t{}\t{}\t{}\t{:.2f}%'.format('TruSeq', ad['truseq'], n_max, d.get('truseq', 0), d.get('truseq', 0)/n_max * 100),
            '{}\t{}\t{}\t{}\t{:.2f}%'.format('Nextera', ad['nextera'], n_max, d.get('nextera', 0), d.get('nextera', 0)/n_max * 100),
            '{}\t{}\t{}\t{}\t{:.2f}%'.format('smallRNA', ad['smallrna'], n_max, d.get('smallrna', 0), d.get('smallrna', 0)/n_max * 100)
            ])

        print(msg)

        # sort
        ds = sorted(d.items(), key=lambda kv: kv[1])
        dd = collections.OrderedDict(ds)
        return dd


    def calFreq(self, x):
        """Calculate the frequency of list
        return dataframe

        index count
        """
        if isinstance(x, list):
            var, freq = np.unique(x, return_counts=True)
            df = pd.DataFrame(data=freq, index=var, columns=['count']).reset_index()
            df.columns = ['length', 'count']
        else:
            df = pd.DataFrame(columns=['length', 'count'])
        return df


    def len_dist(self, n_max=0, csv_file=None):
        """Calculate the length distribution of the fx
        fq
        fa

        >fragsize.csv
        length count
        """
        chunk = 1000000
        counter = 0
        fragSizes = []
        frames = []
        with xopen(self.input) as r:
            for name, seq, qual, comment in self.readfq(r):
                counter += 1
                fragSizes.append(len(seq))
                # last record
                if n_max > 0 and counter  >= n_max:
                    log.info('Stopped at limit: {}'.format(n_max))
                    break # stop
                # chunk
                if counter > 0 and counter%chunk == 0:
                    frames.append(self.calFreq(fragSizes))
                    fragSizes = [] # empty
                    log.info('{} : {}'.format('Processed', counter))
            # last chunk
            if len(fragSizes) > 0:
                frames.append(self.calFreq(fragSizes))
                fragSizes = [] # empty
                log.info('{} : {}'.format('Processed', counter))
        # overall
        df = pd.concat(frames, axis=0).groupby(['length']).sum().reset_index()
        df['id'] = fq_name(self.input)
        # save to file
        if isinstance(csv_file, str):
            if os.path.exists(os.path.dirname(csv_file)):
                df.to_csv(csv_file, index=False)
        # output
        return df
        


