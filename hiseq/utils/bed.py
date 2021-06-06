#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function, processing BED files

functions/classes:

Bed
BedOverlap
bed2gtf
bed_to_saf
TableReader

1. convert formats between: GTF, GFF, 
"""

import os
import pybedtools
from hiseq.utils.helper import *
import sys
import tempfile
import pathlib
import pybedtools
from shutil import which
from itertools import combinations
from itertools import count
from hiseq.utils.utils import convert_image, log, update_obj, Config, \
    run_shell_cmd
from hiseq.utils.file import check_path, file_exists, file_nrows, file_prefix, \
    remove_file, Genome, read_lines
from hiseq.utils.bam import Bam
from hiseq.utils.featurecounts import FeatureCounts


def bed_to_saf(file_in, file_out):
    """Convert BED to SAF format, for featureCounts
    GeneID Chr Start End Strand

    see: https://www.biostars.org/p/228636/#319624
    """
    if file_exists(file_out) and not self.overwrite:
        log.info('bed_to_saf() skipped, file exists: {}'.format(file_out))
    else:
        try:
            with open(file_in, 'rt') as r, open(file_out, 'wt') as w:
                for line in r:
                    tabs = line.strip().split('\t')
                    chr, start, end = tabs[:3]
                    # strand
                    strand = '.'
                    if len(tabs) > 5:
                        s = tabs[5]
                        if s in ['+', '-', '.']:
                            strand = s
                    # id
                    if len(tabs) > 3:
                        name = os.path.basename(tabs[3])
                    else:
                        name = '_'.join([chr, start, end, strand])
                    # output
                    w.write('\t'.join([name, chr, start, end, strand])+'\n')
        except:
            log.error('bed_to_saf() failed, see: {}'.format(file_out))
    return file_out



###############################################################################
## code from bx-python package
## 
## To-Do, not finish 
###############################################################################
class TableReader(object):
    """Reader for tab data
    from bx-python package
    
    Parameters
    ----------
    input : str 
    
    return_header : bool
    
    return_comments : bool
    
    force_header : None
    
    comment_lines_startswith : str
    """
    def __init__(self, input, return_header=True, return_comments=True,
        force_header=None, comment_lines_startswith=['#']):
        self.input = input
        # options
        self.return_comments = return_comments
        self.return_header = return_header
        self.linenum = 0 # init
        self.header = force_header
        self.comment_lines_startswith = comment_lines_startswith
        self.input_iter = iter(input)


    def __iter__(self):
        return self

    def next(self):
        line = self.input_iter.next()
        self.linenum += 1
        line = line.rstrip('\r\n')
        # for blank lines
        # add '#'
        if line.strip() == '':
            if self.return_comments:
                return Comment(line)
            else:
                return self.next()
        # comment lines
        for comment_line_start in self.comment_lines_startswith:
            if line.startswith(comment_line_start):
                # ithe first line is header
                if self.header is None and self.linenum == 1:
                    self.header = self.parse_header(line)
                    if self.return_header:
                        return self.header
                    else:
                        return self.next()
                else:
                    if self.return_comments:
                        return self.parse_comment(line)
                    else:
                        return self.next()
        # not a comment
        try:
            return self.parse_row(line)
        except: # ParseError, err:
            err.linenum = self.linenum
            raise err


    def parse_header(self, line):
        if line.startswith('#'):
            fields = line.lstrip('#').split('\t')
        else:
            fields = line.split('\t')

        return Header(fields)


    def parse_comment(self, line):
        return Comment(line)


    def parse_row(self, line):
        return TableRow(self, line.split('\t'))


class TableRow(object):
    """A row of a table
    """
    def __init__(self, reader, fields):
        self.reader = reader
        self.fields = fields

    def __getitem__(self, key):
        if type(key) == int:
            return self.fields[key]
        elif type(key) == str:
            if self.reader.header:
                if self.reader.header.field_to_column[key] > len(self.fileds):
                    return None
                return self.fileds[self.reader.header.field_to_column[key]]
            else:
                raise TypeError('column names only supported for files with \
                    headers')
        else:
            raise TypeError('field indices must be integers or strings')

    @property
    def fieldnames(self):
        return self.reader.header.fields

    def __str__(self):
        return '\t'.join(self.fields)


class Header(object):
    """
    Header of a table, contains column names and a mapping from them
    to column indexes
    """
    def __init__(self, fields):
        self.set_fields(fields)

    def set_fields(self, fields):
        self.fields = fields
        self.field_to_column = dict(zip(fields, count()))

    def __getitem__(self, key):
        if type(key) == int:
            return self.fields[key]
        elif type(key) == str:
            if key in self.fields_to_column:
                return key
        else:
            raise TypeError('field indices must be integers or strings')

    def __str__(self):
        return '#' + '\t'.join(self.fields)


class Comment(object):
    """
    Comment lines in table input
    start with '#'
    """
    def __init__(self, line):
        self.line = line

    def __str__(self):
        if self.line.startswith('#'):
            return self.line
        return '#' + self.line


class ParseError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args)
        self.linnum = kwargs.get('linenum', None)

    def __str__(self):
        if self.linenum:
            return Exception.__str__(self) + ' on line ' + str(self.linenum)
        else:
            return Exception.__str__(self)



class Bed(object):
    def __init__(self, fields, reader=None):
        """Convertions for Bed file
        """
        if not reader is None:
            TableRow.__init__(self, reader, fields)
        else:
            self.fields = fields

        self.nfields = len(fields)
        if self.nfields < 3:
            raise TypeError('Not enough fields, at least 3')

    def __str__(self):
        return '\t'.join(map(str, self.fields))

    def __cmp__(self, other):
        return cmp(self.fields, other.fields)

    def get_nfields(self):
        """
        number of fields
        """
        return self.nfields

    def get_chrom(self):
        """
        chrom
        """
        return self.__getitem__('chrom')

    def get_start(self):
        """
        chromStart position
        """
        return int(self.__getitem__('chromStart'))

    def get_end(self):
        """
        chromEnd position
        """
        return int(self.__getitem__('chromEnd'))

    def get_name(self):
        return self.__getitem__('name')

    def get_score(self):
        return self.__getitem__('score')

    def get_strand(self):
        return self.__getitem__('strand')

    def get_thickStart(self):
        return int(self.__getitem__('thickStart'))    

    def get_thickEnd(self):
        return int(self.__getitem__('thickEnd'))

    def get_blockCount(self):
        return int(self.__getitem__('blockCount'))

    def get_blockSies(self):
        return int(self.__getitem__('blockSizes'))    

    def copy(self):
        """
        Create new Bed record 
        """
        return Bed(list(self.fields), self.reader)

    def get_length(self):
        return self.get_end() - self.get_start()

    def find_middle(self):
        length = self.get_length()
        return self.get_start() + length/2 # int

    def find_midpoint_interval(self):
        mid = self.find_middle()
        if self.get_strand() == '+':
            s, e = (mid, s + 1)
        else:
            s, e = (mid, mid - 1)
        return (s, e)

    def extend_twosides(self, extent=100, by=True, chrom2size=None):
        return self.extend('twosizes', extend, by, chrom2size)

    def extend_3end(self, extend=100, by=True, chrom2size=None):
        if self.get_strand() == '+':
            direction = 'right'
        else:
            direction = 'left'
        return self.extend(direction, extend, by, chrom2size)

    def extend_5end(self, extend=100, by=True, chrom2size=None):
        if self.get_strand() == '+':
            direction = 'left'
        else:
            direction = 'right'
        return self.extend(direction, extend, by, chrom2size)

    def extend(self, direction, extent=100, by=True, chrom2size=None):
        """
        Return new Bed record with extension
        direction: twosides|left|right
        extlen: extended length
        by: True, add extlen according to direction;
            False, only extend when current length < extlen;
        chrom2size: if True, chop the extended record
        """
        assert direction in ['twosides', 'left', 'right']
        bed_ext = self.copy()
        if by:
            if direction == 'twosides':
                bed_ext.fields[1] = str(self.get_start() - extlen)
                bed_ext.fields[2] = str(self.get_end() + extlen)
            elif direction == 'left':
                bed_ext.fields[1] = str(self.get_start() - extlen)
            else:
                bed_ext.fields[2] = str(self.get_end() + extlen)
        else:
            if self.get_length() >= extlen:
                return bed_ext
            else:
                if direction == 'twosides':
                    middle = self.find_middle()
                    leftlen = int(extlen / 2)
                    bed_ext.fields[1] = str(middle - leftlen)
                    bed_ext.fields[2] = str(middle + extlen - leftlen)
                elif direction == 'left':
                    bed_ext.fields[1] = str(self.get_end() - extlen)
                else:
                    bed_ext.fields[2] = str(self.get_start() + extlen)

        if chrom2size:
            bed_ext.chop(chrom2size)
        return bed_ext

    def chop(self, chrom2size):
        """Chop start, end so as to not pass chromStart, chromEnd"""
        chrom = self.get_chrom()
        assert chrom in chrom2size
        self.fields[1] = str(max(0, self.get_start()))
        self.fields[2] = str(min(chrom2size[chrom], self.get_end()))



class BedReader(TableReader):
    """
    Reader for BED format
    Reader for BED file
    BED3, BED6, BED12,
    UCSC: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    """
    fieldNames = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
        'blockStarts']

    def __init__(self, fhd, discard_first_column=False, return_header=False,
                 return_comments=False, force_header=None,
                 comment_lines_startswith = ["#", "track ", "browser"]):

        TableReader.__init__(self, fhd, return_header, return_comments,
                             force_header, comment_lines_startswith)
        self.discard_first_column = discard_first_column
        if not self.header:
            self.header = Header(BedReader.fieldNames)

    def parse_row(self, line):
        """return Bed object"""
        fields = line.split("\t")
        if self.discard_first_column:
            fields.pop(0)
        return Bed(fields, self) #self as argument reader
    
    
class BedOverlap(object):
    """
    pybedtools API, to calculate the overlaps between bed files
    pybedtools.contrib.venn_maker.venn_maker(beds, names, figure_filename, script_filename, run=True)
    make venn plot
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'peak_list': None,
            'outdir': None,
            'flag': False,
            'prefix': 'peak_overlap',
            'overwrite': False
        }
        self = update_obj(self, args_init, force=False)
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        # update
        self.init_peaks()
        self.init_files()


    def init_peaks(self):
        p = []
        if isinstance(self.peak_list, str):
            p_ext = os.path.splitext(self.peak_list)[1]
            if p_ext == '.txt':
                p = read_lines(self.peak_list, comment='#')
            elif p_ext in ['.bed', '.narrowPeak']:
                p = [self.peak_list]
            else:
                log.error('unknown format, {}'.format(p_ext))
        elif isinstance(self.peak_list, list):
            p = self.peak_list
        else:
            log.error('illegal peak, expect str,list, got {}'.format(
                type(self.peak).__name__))
        # check all
        p = [i for i in p if file_nrows(i) > 0] # filter by rows
        if not all([
            isinstance(p, list),
            all(file_exists(p)),
            len(p) > 1,
        ]):
            log.error('illegal peak, [>=2 files; >0 peaks]')
        # overlap between 2-4 bed files
        if len(p) > 4:
            p = p[:4] # at most 4 files
            log.error('supported for <=4 files, choose the first 4')
        self.peak_list = p # assign


    def init_files(self):
        self.peak_names = file_prefix(self.peak_list)
        if not isinstance(self.prefix, str):
            self.prefix = 'peak_overlap'
        # files
        self.config_yaml = os.path.join(self.outdir, 'config.yaml')
        self.tiff = os.path.join(self.outdir, self.prefix + '.tiff')
        self.png = os.path.join(self.outdir, self.prefix + '.png')
        self.venn_R = os.path.join(self.outdir, self.prefix + '.venn.R')


    def run_overlap(self):
        plt = pybedtools.contrib.venn_maker.venn_maker
        if os.path.exists(self.tiff) and self.overwrite is False:
            log.info('BedOverlap() skipped, file exists')
        if len(self.peak_list) > 1:
            log.info('Calculating overlaps between BED files')
            try:
                plt(self.peak_list, self.peak_names, figure_filename=self.tiff,
                    script_filename=self.venn_R, run=True)
            except:
                log.error('run_overlap() failed')
        # tiff -> png
        if os.path.exists(self.png) and self.overwrite is False:
            pass
        else:
            log.info('Coverting Tiff to png')
            if os.path.exists(self.tiff):
                convert_image(self.tiff, 'PNG')

 
    def run(self):
        self.run_overlap()
        

class PeakFRiP(object):
    """
    Calculate the FRiP 
    see ENCODE: https://www.encodeproject.org/data-standards/terms/#enrichment
    1. bedtools intersect
    2. featureCounts
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'peak': None,
            'bam': None,
            'outdir': None,
            'method': 'bedtools', #option: featureCounts
            'genome': None,
            'gsize': None,
        }
        self = update_obj(self, args_init, force=False)
        if not all([
            isinstance(self.peak, str),
            isinstance(self.bam, str),
            file_exists(self.peak),
            file_exists(self.bam),
            self.method in ['bedtools', 'featureCounts'],
        ]):
            raise ValueError('PeakFRiP() failed, check arguments')
        if not isinstance(self.outdir, str):
            self.outdir = self._tmp(dir=True)
        check_path(self.outdir, create_dirs=True)
        # get gsize file
        if isinstance(self.genome, str):
            self.gsize = Genome(self.genome).get_fasize()
        

    def run_bedtools(self):
        """
        see: https://www.biostars.org/p/337872/#338646
        sort -k1 -k2,2n -o peak.bed peak.bed
        bedtools intersect -c -a peak.bed -b file.bam | awk '{i+=$n}END{print i}'
        """
        total = Bam(self.bam).count()
#         rip_txt = self._tmp(delete=False)
        rip_txt = os.path.join(self.outdir, 'frip_matrix.txt')
        if file_exists(self.gsize):
            args_sorted = '--sorted -g {}'.format(self.gsize)
        else:
            args_sorted = ''
        cmd = ' '.join([
            'sort -k1,1 -k2,2n -o {} {}'.format(self.peak, self.peak),
            '&& bedtools intersect -c',
            '-a {} -b {}'.format(self.peak, self.bam),
            args_sorted,
            r"| awk '{i+=$NF}END{print i}'",
            '> {}'.format(rip_txt)
            ])
        run_shell_cmd(cmd)
        with open(rip_txt, 'rt') as r:
            rip = int(r.read().strip())
        if total > 0:
            frip = round(rip/total, 4) #  '{:.2f}'.format(rip/total*100)
        else:
            frip = 0
        # fragments (PE)
        if Bam(self.bam).is_paired():
            total = int(total / 2.0)
            rip = int(rip / 2.0)
        remove_file(rip_txt, ask=False)
        # output
        out = {
            'index': os.path.basename(self.bam),
            'total': total,
            'map': rip,
            'pct': frip,
        }
        return out


    def run_featureCounts(self):
        """
        see: https://www.biostars.org/p/337872/#337890
        -p , fragments
        """
        fc = FeatureCounts(gtf=self.peak, bam_list=self.bam, outdir=self.outdir)
        fc.run()
        df = Config().load(fc.summary_json) # multiple bam
        return list(df.values())[0] # index,total,map,pct #first one; dict

    
    def _tmp(self, delete=True, dir=False, suffix='.txt'):
        """
        Create a tmp file/dir
        """
        if dir:
            tmp = tempfile.TemporaryDirectory()
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix, delete=delete)
        return tmp.name



    def run(self):
        if self.method == 'bedtools':
            out = self.run_bedtools()
        elif self.method == 'featureCounts':
            out = self.run_featureCounts()
        else:
            out = None
            log.error('method, unknown, [bedtools|featureCounts], got {}'.format(
                self.method))
        return out


class PeakIDR(object):
    """
    Check IDR
    Irreproducibility Discovery Rate
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'idr_cmd': which('idr'),
            'peak_list': None,
            'outdir': None,
            'prefix': None,
            'input_type': 'narrowPeak', # broadPeak, bed, gff
            'cor_method': 'pearson', # spearman
            'overwrite': False,
            'flag': True # run all
        }
        self = update_obj(self, args_init, force=False)
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.init_peaks()
        check_path(self.outdir, create_dirs=True)
        self.config_yaml = os.path.join(self.outdir, self.prefix+'.config.yaml')
        Config().dump(self.__dict__, self.config_yaml)

    
    def init_peaks(self):
        p = []
        if isinstance(self.peak_list, str):
            p_ext = os.path.splitext(self.peak_list)[1]
            if p_ext == '.txt':
                p = read_lines(self.peak_list, comment='#')
            elif p_ext in ['.bed', '.narrowPeak']:
                p = [self.peak_list]
            else:
                log.error('unknown format, {}'.format(p_ext))
        elif isinstance(self.peak_list, list):
            p = self.peak_list
        else:
            log.error('illegal peak, expect str,list, got {}'.format(
                type(self.peak).__name__))
        # check all
        p = [i for i in p if file_nrows(i) > 100] # filter by rows
        if not all([
            isinstance(p, list),
            all(file_exists(p)),
            len(p) > 1,
        ]):
            log.error('illegal peak, [>=2 files; >100peaks]')
        self.peak_list = p # assign


    def run_pair_idr(self, peakA=None, peakB=None):
        """
        Compute IDR for two group of peaks
        """
        if not all([
            isinstance(peakA, str),
            isinstance(peakB, str),
            file_exists(peakA),
            file_exists(peakB),
        ]):
            log.error('run_pair_idr() skipped, illegal peakA/peakB')
            return None
        # generate prefix/names
        pA, pB = file_prefix([peakA, peakB])
        prefix = '{}.vs.{}.idr'.format(pA, pB)
        if isinstance(self.prefix, str):
            prefix = self.prefix + prefix
        # files
        cmd_txt = os.path.join(self.outdir, prefix + '.cmd.sh')
        idr_txt = os.path.join(self.outdir, prefix + '.txt')
        idr_png = idr_txt + '.png'
        idr_log = os.path.join(self.outdir, prefix + '.log')
        # command
        cmd = ' '.join([
            '{} --samples {} {}'.format(self.idr_cmd, peakA, peakB),
            '--input-file-type {}'.format(self.input_type),
            '--plot', #.format(idr_png),
            '--output-file {}'.format(idr_txt),
            '--log-output-file {}'.format(idr_log)
            ])
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(idr_png) and not self.overwrite:
            log.warning('run_idr() skipped, file exists: {}'.format(idr_png))
        else:
            try:
                _, stdout, stderr = run_shell_cmd(cmd)
                with open(idr_log, 'wt') as w:
                    w.write(stdout + '\n' + stderr + '\n')
            except:
                log.error('idr command failed')
        if not os.path.exists(idr_png):
            log.error('Peak().idr() failed: {}'.format(idr_log))


    def run(self):
        """
        Run A, B
        """
        # prepare config
        msg = '\n'.join([
            '-'*80,
            '{:>14} : {}'.format('program', 'hiseq.utils.bed.PeakIDR()'),
            '{:>14} : {}'.format('peaks', self.peak_list),
            '{:>14} : {}'.format('outdir', self.outdir),
            '{:>14} : {}'.format('prefix', self.prefix),
            '{:>14} : {}'.format('peak_type', self.input_type),
            '{:>14} : {}'.format('cor_method', self.cor_method),
            '{:>14} : {}'.format('n_peaks', len(self.peak_list)),
            '-'*80,
        ])
        print(msg)
        if len(self.peak_list) > 1:
            for peakA, peakB in combinations(self.peak_list, 2):
                self.run_pair_idr(peakA, peakB)


#