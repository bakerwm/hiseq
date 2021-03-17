#!/usr/bin/env python3 

"""For aligner index
1. search/fetch index
2. validate index
"""

import os
from hiseq.utils.helper import *
from utils import *

class AlignIndex(object):
    """Validate the index for aligner
    This is a general function/class for aligner index
    functions:
    1. is_valid()
    2. guess_aligner()
    3. index_name()
    4. index_size()
    
    Parameters
    ----------
    index : str
        Path to the index for aligner
    
    aligner : Nonr or str
        The aligner (supported: ...)
        
        
    >>> p = AlignIndex(x).is_valid()
    >>> p = AlignIndex(x, 'bowtie').is_valid()
    >>> AlignIndex(x).index_name()
    >>> AlignIndex(x).index_size(out_file='chrom.sizes.txt')
    """
    def __init__(self, index=None, aligner=None, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.index = index
        self.aligner = aligner
        self.supported_aligner = [
            'bowtie', 'bowtie2', 'hisat2', 'bwa', 'star',
            ]


    def is_star_index(self, index=None):
        """The index is for STAR
        !!!to-do!!!
        make sure the STAR version is consistent

        require files:
        SAindex, Genome, SA, ...
        """
        if index is None:
            index = self.index
        if index is None:
#             log.error('index not valid, expect str, got NoneType')
            out = False
        f = [
            'SAindex', 'Genome', 'SA', 'genomeParameters.txt',
            'chrLength.txt', 'chrNameLength.txt', 'chrName.txt',
            ]
        if isinstance(index, str):
            f_list = [os.path.join(index, i) for i in f]
            out = check_file(f_list, check_empty=True)
        else:
            out = False
        return out


    def is_bowtie_index(self, index=None):
        """The index is for bowtie
        check the files
        .[1234].ebwt, .rev.[12].ebwt
        """
        if index is None:
            index = self.index
        if index is None:
#             log.error('index not valid, expect str, got NoneType')
            out = False
        f = ['.' + i + '.ebwt' for i in [
            '1', '2', '3', '4', 'rev.1', 'rev.2'
            ]]
        if isinstance(index, str):
            f_list = [index+i for i in f]
            out = check_file(f_list, check_empty=True)
        else:
            out = False
        return out


    def is_bowtie2_index(self, index=None):
        """The index is for bowtie2
        check the files
        .[1234].bt2, .rev.[12].bt2
        """
        if index is None:
            index = self.index
        if index is None:
#             log.error('index not valid, expect str, got NoneType')
            out = False
        f = ['.' + i + '.bt2' for i in [
            '1', '2', '3', '4', 'rev.1', 'rev.2'
            ]]
        if isinstance(index, str):
            f_list = [index+i for i in f]
            out = check_file(f_list, check_empty=True)
        else:
            out = False
        return out


    def is_hisat2_index(self, index=None):
        """The index is for hisat2
        check the files
        .[1234].ht2, .rev.[12].ht2
        """
        if index is None:
            index = self.index
        if index is None:
#             log.error('index not valid, expect str, got NoneType')
            out = False
        f = ['.' + i + '.ht2' for i in [
            '1', '2', '3', '4', 'rev.1', 'rev.2'
            ]]
        if isinstance(index, str):
            f_list = [index+i for i in f]
            out = check_file(f_list, check_empty=True)
        else:
            out = False
        return out


    def is_bwa_index(self, index=None):
        """The index is for bowtie
        check the files
        .[1234].ebwt, .rev.[12].ebwt
        """
        if index is None:
            index = self.index
        if index is None:
#             log.error('index not valid, expect str, got NoneType')
            out = False
        f = ['.sa', '.amb', '.ann', '.pac', '.bwt']
        if isinstance(index, str):
            f_list = [index+i for i in f]
            out = check_file(f_list, check_empty=True)
        else:
            out = False
        return out


    def guess_aligner(self, index=None):
        """Guess the aligner of the index
        Rules:
        see: is_{}_index()
        """
        fn = {
            'bowtie': self.is_bowtie_index,
            'bowtie2': self.is_bowtie2_index,
            'hisat2': self.is_hisat2_index,
            'bwa': self.is_bwa_index,
            'star': self.is_star_index
        }
        for a in self.supported_aligner:
            is_index = fn.get(a, None)
            if is_index(index):
                aligner = a
                break
            else:
                aligner = None
        return aligner


    def is_valid(self, index=None):
        """The input index is valid"""
        # match the aligner
        if isinstance(self.aligner, str):
            out = self.guess_aligner(index) == self.aligner.lower()
        else:
            out = self.guess_aligner(index) is not None
        return out


    def index_name(self, index=None):
        """Fetch the name of the index
        basename : bowtie, bowtie2, hisat2, bwa
        is_dir, STAR
        """
        if index is None:
            index = self.index
        if index is None:
#             log.error('index not valid, expect str, got NoneType')
            out = None
        if self.is_valid(index):
            index = index.rstrip('/')
            out = os.path.basename(index)
        else:
            log.error('not valid index: {}'.format(index))
            out = None
        return out


    def _tmp(self, is_dir=False, suffix='.txt'):
        """
        Create a tmp file to save json object
        """
        if is_dir:
            tmp = tempfile.TemporaryDirectory(prefix='tmp')
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
                delete=False)
        return tmp.name


    def inspect_index(self, index=None, out_file=None):
        """Inspect the index for {bowtie|bowtie2|hisat2}
        bowtie-inspect -s index
        """        
        if index is None:
            index = self.index
#         if index is None:
#             log.error('index not valid, expect str, got NoneType')
        out = None
        if self.is_valid(index):
            tmp_f1 = self._tmp()
            tmp_f2 = out_file if isinstance(out_file, str) else self._tmp()
            aligner = self.guess_aligner(index)
            if aligner in ['bowtie', 'bowtie2', 'hisat2']:
                cmd = ' '.join([
                    '{}-inspect'.format(aligner),
                    '-s {}'.format(index),
                    '> {}'.format(tmp_f1)
                    ])
                try:
                    os.system(cmd)
                except:
                    log.error('failed to run: {}-inspect'.foramt(aligner))
                # processing the file
                # gsize file
                gsizes = []
                with open(tmp_f1) as r:
                    for line in r:
                        if line.startswith('Sequence'):
                            _,name,size = line.strip().split('\t')
                            gsizes.append((name, size))
                if len(gsizes) > 0:
                    with open(tmp_f2, 'wt') as w:
                        w.write('\n'.join(['{}\t{}'.format(a,b) for a,b in 
                            gsizes])+'\n')
                    out = tmp_f2
                    file_remove(tmp_f1, ask=False)
                else:
                    out = None
        return out


    def index_size(self, index=None, out_file=None):
        """Extract the chr size from index
        bowtie:  bowtie-inspect -s <index>
        bowtie2: bowtie-inspect -s <index>
        STAR:  chrNameLength.txt
        """
        if index is None:
            index = self.index
        if index is None:
#             log.error('index not valid, expect str, got NoneType')
            out = None
        if self.is_valid(index):
            if self.is_star_index(index):
                gsize = os.path.join(index, 'chrNameLength.txt')
            elif self.is_bowtie_index(index) or \
                self.is_bowtie2_index(index) or \
                self.is_hisat2_index(index):
                gsize = self.inspect_index(index, out_file)
            else:
                log.error('index_size() not support the index: {}'.format(
                    index))
                gsize = None
        # output
        s = 0
        if file_exists(gsize):
            with open(gsize) as r:
                for line in r:
                    s += eval(line.strip().split('\t')[1])
        # return s if s > 0 else None
        return s
        


def fetch_index(genome, group=None, aligner='bowtie', genome_path=None):
    """fetch index
    Genome, tags, aligner
    
    Parameters
    ----------
    genome : str
        The UCSC name of the genome, [dm6, mm10, hg38, GRCh38]

    group : str
        The group of the index, [genome, rRNA, rRNA, chrM, MT, MT_trRNA, ...]
        if group=`None`, return the name of the genome

    aligner : str
        The aligner, [bowtie, bowtie2, star, bwa, ...] 

    genome_path : str
        The root of the genome data

    Structure of genome_path:
    default: {HOME}/data/genome/{genome_version}/{aligner}/

    ## bowtie/bowtie2/hisat2/...
    path-to-genome/
        |- Bowtie_index /
            |- genome
            |- rRNA
            |- MT_trRNA
            |- transposon
            |- piRNA_cluster

    ## STAR
    path-to-genome/
        |- Bowtie_index /
            |- genome/
            |- rRNA/
            |- MT_trRNA/
            |- transposon/
            |- piRNA_cluster/
    """
    supported_genomes = ['dm3', 'dm6', 'mm9', 'mm10', 'hg19', 'hg38', 'GRCh38']
    supported_aligner = ['bowtie', 'bowtie2', 'bwa', 'star', 'hisat2']
    # check arguments
    tag_err = False
    if not genome in supported_genomes:
        log.error('genome={}, unknown, supported: {}'.format(
            genome, supported_genome))
        tag_err = True
    if not aligner in supported_aligner:
        log.error('aligner={}, unknown, supported: {}'.format(
            aligner, supported_aligner))
        tag_err = True
    if genome_path is None:
        genome_path =  os.path.join(str(pathlib.Path.home()), 'data', 'genome')
    if not os.path.exists(genome_path):
        log.error('genome_path={} not exists'.format(genome_path))
        tag_err = True
    if group is None:
        group = genome
    if not isinstance(group, str):
        log.error('group={} expect str, got {}'.format(
            type(group).__name__, group))
        tag_err = True
    if tag_err:
        log.error('invalid arguments: genome={}, group={}, aligner={}\
            genome_path={}'.format(genome, group, aligner, genome_path))
        return None
    # index
    index = os.path.join(
        genome_path, genome, '{}_index'.format(aligner), group)
    # output
    if AlignIndex(index, aligner).is_valid():
        out = index
    else:
        log.error('index not found: genome={}, group={}, aligner={}'.format(
            genome, group, aligner))
        out = None
    return out




