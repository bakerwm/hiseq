#!/usr/bin/env python3


"""
Function:
1. Read index list 
2. mission type: bc/p7/bc+p7

index_table:
name,p7,p5,bc

Criteria:
1. p7+p6+bc no duplication
"""


import os
import sys
import Levenshtein as lev # distance
import logging
from hiseq.utils.helper import *
from hiseq.utils.helper import update_obj




class IndexTable(object):
    """
    Load index, csv format
    #file_name, p7, p5, barcode
    
    Example:
    >>> args = {
    'index_table': 'yy58.csv',
    'mismatch': 1
    }

    >>> IndexTable(**args).run()
    --------------------------------------------------------------------------------
    Check the index table:
             No. samples | 33      
            max mismatch | 1       
         filename unique | ok      
            index unique | ok      
                p7_index | ok      
                p5_index | skipped 
                 barcode | ok      
    --------------------------------------------------------------------------------
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.read_index()


    def init_args(self):
        args_init = {
            'index_table': None,
            'mismatch': 0
        }
        self = update_obj(self, args_init, force=False)
        if not os.path.exists(self.index_table):
            log.error('index file not exists')
        if not self.mismatch in range(4):
            log.warning('mismatch illegal, use mm==0 instead')
            self.mismatch = 0
        # available index
        self.idx = illumina_index(seq_type='all', by_name=False)


    def is_valid_index(self, x):
        """
        x, str
        """
        # return x == 'NULL' or re.match('^[ACGTN]+$', x) is not None
        return x == 'NULL' or x in self.idx

        
    def str_distance(self, x, y, partial=True):
        """
        Check distance between a, b
        """
        if partial:
            x = x[:len(y)]
            y = y[:len(x)]
        return lev.distance(x, y)


    def is_compatiable(self, x, y, mm=0):
        """
        Check if the index x and y, compatiable
        """
        if self.is_valid_index(x) and self.is_valid_index(y):
            out = self.str_distance(x, y) > mm
        elif x == 'NULL' or y == 'NULL':
            out = True
        else:
            out = False
        return out
    
    
    def is_compatiable2(self, x):
        """
        Check if the index in x, compatiable between each other
        """
        if isinstance(x, list):
            d = []
            for i, j in combinations(x, 2):
                d.append(self.is_compatiable(i, j, self.mismatch))
            out = all(d)
        else:
            out = False
        return out

    
    def read_index(self):
        """Check: p7, p5, bc"""
        p7_list = []
        p5_list = []
        bc_list = []        
        d = {}
        with open(self.index_table) as r:
            for line in r:
                if line.startswith('#'):
                    continue
                s = re.split('[,\s\t]', line.strip())
                if not len(s) == 4:
                    log.error('illegal format, expect:name,p7,p5,bc')
                    break
                name,p7,p5,bc = s
                p7,p5,bc = [i.upper() for i in [p7, p5, bc]] # upper case
                if name in d:
                    log.error('duplicate filenames: {}'.name)
                    break
                # save each index
                if not p7 == 'NULL':
                    p7_list.append(p7)
                if not p5 == 'NULL':
                    p5_list.append(p5)
                if not bc == 'NULL':
                    bc_list.append(bc)                
                d[name] = (p7,p5,bc)
        # update:
        args = {
            'samples': d,
            'p7_list': list(set(p7_list)),
            'p5_list': list(set(p5_list)),
            'bc_list': list(set(bc_list)),
        }
        self = update_obj(self, args, force=True)
        # return [d, p7_list, p5_list, bc_list]
    
    
    def check_index(self):
        # d, p7_list, p5_list, bc_list = self.read_index()
        # for p7
        p7_flag = 'skipped' if len(self.p7_list) == 0 \
            else 'ok' if self.is_compatiable2(self.p7_list) else 'failed'
        p5_flag = 'skipped' if len(self.p5_list) == 0 \
            else 'ok' if self.is_compatiable2(self.p5_list) else 'failed'
        bc_flag = 'skipped' if len(self.bc_list) == 0 \
            else 'ok' if self.is_compatiable2(self.bc_list) else 'failed'
        # to-do: bit compute
        # name, nodup
        name_list = list(self.samples.keys())
        name_flag = 'ok' if len(name_list) == len(set(name_list)) else 'failed'
        # for all index
        idx_list = list(self.samples.values())
        idx_flag = 'ok' if len(idx_list) == len(set(idx_list)) else 'failed'
        # output
        return (name_flag, idx_flag, p7_flag, p5_flag, bc_flag)


    def run(self):
        name_flag, idx_flag, p7_flag, p5_flag, bc_flag = self.check_index()
        # message
        msg = '\n'.join([
            '-'*80,
            'Check index table:',
            '{:>20} | {:<8}'.format('No. samples', len(self.samples)),
            '{:>20} | {:<8}'.format('max mismatch', self.mismatch),
            '{:>20} | {:<8}'.format('filename unique', name_flag),
            '{:>20} | {:<8}'.format('index unique', idx_flag),
            '{:>20} | {:<8}'.format('p7_index', p7_flag),
            '{:>20} | {:<8}'.format('p5_index', p5_flag),
            '{:>20} | {:<8}'.format('barcode', bc_flag),
            '-'*80
        ])
        print(msg)
        out = 'failed' in [name_flag, idx_flag, p7_flag, p5_flag, bc_flag]
        if out:
            log.error('index failed, try mismatch={}, run again.'.format(
                self.mismatch - 1))
        return not out

    

def illumina_index(seq_type="truseq", by_name=True):
    """
    Parameters
    ----------
    seq_type: str
        The group of seq_type, options: [all, truseq, nextera, bc]
    
    by_name: bool
        Save the name in key
    """
    index_list = """
        TruSeq_Index1,ATCACG
        TruSeq_Index2,CGATGT
        TruSeq_Index3,TTAGGC
        TruSeq_Index4,TGACCA
        TruSeq_Index5,ACAGTG
        TruSeq_Index6,GCCAAT
        TruSeq_Index7,CAGATC
        TruSeq_Index8,ACTTGA
        TruSeq_Index9,GATCAG
        TruSeq_Index10,TAGCTT
        TruSeq_Index11,GGCTAC
        TruSeq_Index12,CTTGTA
        TruSeq_Index13,AGTCAA
        TruSeq_Index14,AGTTCC
        TruSeq_Index15,ATGTCA
        TruSeq_Index16,CCGTCC
        TruSeq_Index17,GTAGAG
        TruSeq_Index18,GTCCGC
        TruSeq_Index19,GTGAAA
        TruSeq_Index20,GTGGCC
        TruSeq_Index21,GTTTCG
        TruSeq_Index22,CGTACG
        TruSeq_Index23,GAGTGG
        TruSeq_Index24,GGTAGC
        TruSeq_Index25,ACTGAT
        TruSeq_Index26,ATGAGC
        TruSeq_Index27,ATTCCT
        TruSeq_Index28,CAAAAG
        TruSeq_Index29,CAACTA
        TruSeq_Index30,CACCGG
        TruSeq_Index31,CACGAT
        TruSeq_Index32,CACTCA
        TruSeq_Index33,CAGGCG
        TruSeq_Index34,CATGGC
        TruSeq_Index35,CATTTT
        TruSeq_Index36,CCAACA
        TruSeq_Index37,CGGAAT
        TruSeq_Index38,CTAGCT
        TruSeq_Index39,CTATAC
        TruSeq_Index40,CTCAGA
        TruSeq_Index41,GACGAC
        TruSeq_Index42,TAATCG
        TruSeq_Index43,TACAGC
        TruSeq_Index44,TATAAT
        TruSeq_Index45,TCATTC
        TruSeq_Index46,TCCCGA
        TruSeq_Index47,TCGAAG
        TruSeq_Index48,TCGGCA
        Next_Ad2.1,TAAGGCGA
        Next_Ad2.2,CGTACTAG
        Next_Ad2.3,AGGCAGAA
        Next_Ad2.4,TCCTGAGC
        Next_Ad2.5,GGACTCCT
        Next_Ad2.6,TAGGCATG
        Next_Ad2.7,CTCTCTAC
        Next_Ad2.8,CAGAGAGG
        Next_Ad2.9,GCTACGCT
        Next_Ad2.10,CGAGGCTG
        Next_Ad2.11,AAGAGGCA
        Next_Ad2.12,GTAGAGGA
        Next_Ad2.13,GTCGTGAT
        Next_Ad2.14,ACCACTGT
        Next_Ad2.15,TGGATCTG
        Next_Ad2.16,CCGTTTGT
        Next_Ad2.17,TGCTGGGT
        Next_Ad2.18,GAGGGGTT
        Next_Ad2.19,AGGTTGGG
        Next_Ad2.20,GTGTGGTG
        Next_Ad2.21,TGGGTTTC
        Next_Ad2.22,TGGTCACA
        Next_Ad2.23,TTGACCCT
        Next_Ad2.24,CCACTCCT
        P7_1A,CCTATA
        P7_1B,TGCTAT
        P7_2A,TATACT
        P7_2B,ATCTTC
        P7_3A,CGTGAT
        P7_3B,ACATCG
        P7_4A,GCCTAA
        P7_4B,ATTGGC
        P7_5A,TGGTCA
        P7_5B,CACAGT
        P7_6A,GATCTG
        P7_6B,TCAAGT
        P7_7A,CTGTAC
        P7_7B,AGACTA
        P7_8A,GTAGCC
        P7_8B,TACAAG
        """
    # determine the seq_type
    d = {'truseq': 'TruSeq', 'nextera': 'Next', 'bc': 'P7', 'all': 'all'}
    s = d.get(seq_type.lower(), None)
    # extract index
    idx = index_list.strip().split('\n')
    f = [i.strip().split(',') for i in idx ]
    d1 = {k:v for k,v in f if k.startswith(s) or s == 'all'}
    d2 = {v:k for k,v in f if k.startswith(s) or s == 'all'}
    return d1 if by_name else d2    
    
    