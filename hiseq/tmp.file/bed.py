#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function, processing BED files

1. convert formats between: GTF, GFF, 
"""

import os
import pybedtools
from hiseq.utils.helper import *






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