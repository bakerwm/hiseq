#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function, processing BED files

1. convert formats between: GTF, GFF, 
"""

import os
import pybedtools
from hiseq.utils.helper import *
import sys
from itertools import count






def bed2gtf(file_in, file_out):
    """Convert BED to GTF
    chrom chromStart chromEnd name score strand
    """
    with open(file_in) as r, open(file_out, 'wt') as w:
        for line in r:
            fields = line.strip().split('\t')
            start = int(fields[1]) + 1
            w.write('\t'.join([
                fields[0],
                'BED_file',
                'gene',
                str(start),
                fields[2],
                '.',
                fields[5],
                '.',
                'gene_id "{}"; gene_name "{}"'.format(fields[3], fields[3])
                ]) + '\n')
    return file_out



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

