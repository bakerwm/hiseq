
# -*- coding: utf-8 -*-

"""
Functions for sequence, fastx

"""

import os
import sys
import re
from xopen import xopen
import collections # Fastx().collapse()
from .helper import *


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
        self.input = input
        self.format = self.fx_type(input)
        self.count = self.fq_counter(input) if self.format == 'fastq' else self.fa_counter(input)
        self.cat = 'zcat' if input.endswith('.gz') else 'cat'


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
            for name, seq, qual in self.readfq(r):
                base_from = 'ACGTNacgtn'
                base_to = 'TGCANtgcan'
                tab = str.maketrans(base_from, base_to)
                seq = seq.translate(tab)[::-1]
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

        seqkit_exe = which('seqkit')
        if seqkit_exe is None:
            raise Exception('command not found in $PATH: seqkit')

        # log.info('extracting sub-sample: n=%d' % n)

        ## warning
        if n > 1000000 or p > 0.1:
            log.warning('choose a smaller number.')

        cmd = '{} | {} sample -n {} > {}'.format(
            self.cat,
            self.input,
            seqkit_exe,
            n,
            out)

        run_shell_cmd(cmd)

        return self.output


    def sample(self, outdir, sample_size=1000000, gzipped=True):
        """
        Create a subsample of input fastq files, default: 1M reads
        Run the whole process for demostration
        """
        s = 4 if self.format == 'fastq' else 2
        nsize = s * sample_size

        # update args
        check_path(outdir)
        
        # subsample
        fx_out = os.path.join(outdir, os.path.basename(self.input))
        
        # gzip output
        if not self.input.endswith('.gz') and gzipped:
            fx_out += '.gz'

        # run
        if os.path.exists(fx_out):
            log.info('file eixsts, {}'.format(fx_out))
        else:
            with xopen(self.input, 'rt') as r, xopen(fx_out, 'wt') as w:
                i = 0 # counter
                for line in r:
                    i += 1
                    if i > nsize: # fastq file: 4 lines per read
                        break
                    w.write(line)


    def fa2fq(self, out):
        """
        Convert fasta to fastq
        quality='J' Phred = 33
        """
        if not self.format == 'fasta':
            raise Exception('fasta file expected, input: {}'.format(self.input))

        with xopen(self.input) as r, xopen(out, 'wt') as w:
            for name, seq, qual in self.readfq(r):
                qual = 'J' * len(seq) # Phred33, 41
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
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fh: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fh: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
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
        return self.count(x) / 4


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

        # merge two funcs
        def cut_cut(x):
            # cut
            x_cut = cut_sub(x)

            # cut to length
            if not cut_to_length == 0:
                x_cut = cut_sub2(x_cut)

            return x_cut

        with xopen(self.input) as r, xopen(out, 'wt') as w:
            for name, seq, qual in self.readfq(r):
                seq = cut_cut(seq)
                if len(seq) < len_min and discard_tooshort:
                    continue # skip
                # write
                if qual is None: # fasta
                    w.write('\n'.join(['>'+name, seq]) + '\n')
                else:
                    qual = cut_cut(qual)
                    w.write('\n'.join(['@'+name, seq, '+', qual]) + '\n')


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


