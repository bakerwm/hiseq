#!/usr/bin/env python3 

"""For aligner index
1. search/fetch index
2. validate index
3. check arguments

Including functions/classes:

AlignIndex
fetch_index
check_index_args

Utilities:
is_supported: check genome, aligner, ...
"""

import os
import pathlib
import tempfile
from hiseq.utils.utils import Config, update_obj, log, is_supported
from hiseq.utils.file import check_file, file_exists, remove_file


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
            out = False
        f = ['.sa', '.amb', '.ann', '.pac', '.bwt']
        if isinstance(index, str):
            f_list = [index+i for i in f]
            out = check_file(f_list, check_empty=True)
        else:
            out = False
        return out


    def is_salmon_index(self, index=None):
        """The index is for salmon
        check the files: info.json
        "SeqHash":
        "NameHash":
        """
        if index is None:
            index = self.index
        out = False
        if isinstance(index, str):
            f = os.path.join(index, 'info.json')
            if file_exists(f):
                df = Config().load(f)
                out = 'SeqHash' in df
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
            'STAR': self.is_star_index,
            'salmon': self.is_salmon_index,
        }
        for a in is_supported(key='supported_aligner', return_values=True):
            is_index = fn.get(a, None)
            if is_index is None:
                continue
            if is_index(index):
                aligner = a
                break
            else:
                aligner = None
        return aligner


    def is_valid(self, index=None):
        """The input index is valid"""
        if index is None:
            index = self.index
        index_aligner = self.guess_aligner(index)
        # match the aligner
        if isinstance(self.aligner, str):
            out = index_aligner == self.aligner
        else:
            out = index_aligner is not None
        return out


    def index_name(self, index=None):
        """Fetch the name of the index
        basename : bowtie, bowtie2, hisat2, bwa
        is_dir, STAR
        """
        if index is None:
            index = self.index
        if index is None:
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
                    log.error('failed to run: {}-inspect'.format(aligner))
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
                    remove_file(tmp_f1, ask=False)
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
            out = None
        # gsize: file, chrom_sizes
        gsize = None
        if self.is_valid(index):
            if self.is_star_index(index):
                gsize = os.path.join(index, 'chrNameLength.txt')
            elif self.is_bowtie_index(index) or \
                self.is_bowtie2_index(index) or \
                self.is_hisat2_index(index):
                gsize = self.inspect_index(index, out_file)
            elif self.is_salmon_index(index):
                f = os.path.join(index, 'info.json')
                if file_exists(f):
                    df = Config().load(f)
                    s = df.get('seq_length', None)
                else:
                    s = 0
            else:
                log.error('index_size() not support the index: {}'.format(
                    index))
        # output
        s = 0
        if file_exists(gsize):
            with open(gsize) as r:
                for line in r:
                    s += eval(line.strip().split('\t')[1])
        # return s if s > 0 else None
        return gsize if out_file else s


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
        The aligner, [bowtie, bowtie2, star, bwa, salmon, ...] 

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
    # check arguments
    tag_err = False
    if not is_supported(genome):
        log.error('genome={}, unknown, supported: {}'.format(
            genome, is_supported(key='supported_genome', return_values=True)))
        tag_err = True
    # if not aligner in supported_aligner:
    if not is_supported(aligner):
        log.error('aligner={}, unknown, supported: {}'.format(
            aligner, is_supported(key='supported_aligner', return_values=True)))
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


def check_index_args(**kwargs):
    """Check the index for aligner
    exists
    valid
    aligner
    name
    ...
    Parameters
    ----------
    index_list (ignore all)
    extra_index (append)
    genome
    genome_index
    spikein
    spikein_index
    to_rRNA
    to_MT
    to_chrM
    to_MT_trRNA
    """
    # default arguments
    aligner = kwargs.get('aligner', None)
    genome = kwargs.get('genome', None)
    genome_index = kwargs.get('genome_index', None)
    spikein = kwargs.get('spikein', None)
    spikein_index = kwargs.get('spikein_index', None)
    index_list = kwargs.get('index_list', None)
    extra_index = kwargs.get('extra_index', None)
    to_rRNA = kwargs.get('to_rRNA', False)
    to_MT = kwargs.get('to_MT', False)
    to_chrM = kwargs.get('to_chrM', False)
    to_MT_trRNA = kwargs.get('to_MT_trRNA', False)
    verbose = kwargs.get('verbose', False)
    # for message
    msg = '\n'.join([
        '-'*80,
        'The arguments for index:',
        '{:>14s} : {}'.format('aligner', aligner),
        '{:>14s} : {}'.format('genome', genome),
        '{:>14s} : {}'.format('genome_index', genome_index),
        '{:>14s} : {}'.format('spikein', spikein),
        '{:>14s} : {}'.format('spikein_index', spikein_index),
        '{:>14s} : {}'.format('index_list', index_list),
        '{:>14s} : {}'.format('extra_index', index_list),
        '{:>14s} : {}'.format('to_rRNA', to_rRNA),
        '{:>14s} : {}'.format('to_chrM', to_chrM),
        '{:>14s} : {}'.format('to_MT_trRNA', to_MT_trRNA),
        '-'*80,
    ])
    if verbose:
        print(msg)
    # index group:
    if to_rRNA:
        group = 'rRNA'
    elif to_MT or to_chrM:
        group = 'chrM'
    elif to_MT_trRNA:
        group = 'MT_trRNA'
    else:
        group = None
    group_index_g = None
    group_index_sp = None
    # level-1
    if isinstance(index_list, list):
        pass
    else:
        index_list = [] # init
        # level-1.1 : extra
        if isinstance(extra_index, str):
            if AlignIndex(extra_index, aligner).is_valid():
                index_list = [extra_index]
            else:
                raise ValueError('extra_index not valid, {}'.format(extra_index))
        else:
            # level-2.1 : spikein
            if spikein_index:
                if AlignIndex(spikein_index, aligner).is_valid():
                    index_list.append(spikein_index)
                else:
                    log.error('not valid spikein_index: {}'.format(spikein_index))
            elif is_supported(spikein):
                spikein_index = fetch_index(spikein, aligner=aligner)
                if group:
                    group_index_sp = fetch_index(
                        spikein, group=group, aligner=aligner
                    )
                    if group_index_sp:
                        index_list.append(group_index_sp)
                if spikein_index:
                    index_list.append(spikein_index)
                else:
                    log.error('AlignIndex({}, {}) not available'.format(
                        spikein, aligner))
            # level-2.2 : genome
            if genome_index:
                if AlignIndex(genome_index, aligner).is_valid():
                    index_list.append(genome_index)
                else:
                    log.error('not valid index: {}'.format(genome_index))
            elif is_supported(genome):
                genome_index = fetch_index(genome, aligner=aligner)
                if group:
                    group_index_g = fetch_index(
                        genome, group=group, aligner=aligner)
                    if group_index_g:
                        index_list.append(group_index_g)
                if genome_index:
                    index_list.append(genome_index)
                else:
                    log.error('AlignIndex({}, {}) not available'.format(
                        genome, aligner))
    # check all
    index_list = [i for i in index_list if 
        AlignIndex(i, aligner).is_valid()]
    index_list_msg = '\n'.join(index_list) if len(index_list) > 0 else '-'
    # msg
    msg_out = '\n'.join([
        '-'*80,
        'The index output:',
        index_list_msg,
        '-'*80,
        ])
    print(msg_out)
    return index_list

