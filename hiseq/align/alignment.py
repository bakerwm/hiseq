# -*- coding:utf-8 -*-

"""
###############
1. unmap files: 
  - bowtie  : _1.fastq, _2.fastq
  - bowtie2 : .1.fastq, .2.fastq




###############


Align fasta/q files to reference sequence

bowtie
bowtie2
STAR
bwa
hisat2
...


arguments:

fq
index
outdir

# optional

unique
mismatch
insertion (bowtie2 -X)

max-multi

# output
sam

unmap

# custom args
arg_str


## standard

input: fq(s), index, args
output: sam, unmap, log
"""

import os
import sys
import re
import shutil
import logging

from hiseq.utils.args import Adapter, args_init
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions


def sam2bam(sam, bam, sort=True, extra_para=''):
    """
    Convert sam to bam, sort by position
    """
    samtools_exe = shutil.which('samtools')

    cmd = '{} view -Subh {} {} | samtools sort -o {} -'.format(
        samtools_exe,
        extra_para,
        sam,
        bam)

    run_shell_cmd(cmd)

    return bam


## top level
@Logger('INFO')
class Alignment(object):

    def __init__(self, **kwargs):
        """
        The top level of alignment function
        A port for actual usage
        
        args:
          - aligner
          - fq
          - outdir
          - fq2
          - index_list
        """
        args = args_init(kwargs, align=True)
        
        if args['fq'] is None:
            raise Exception('{:>10} : --fq no input detected'.format('error'))

        self.kwargs = args


    def run(self):
        bam_list = AlignNFastx(**self.kwargs).run()
        print(bam_list)


## not tested
class bwa(object):

    def __init__(self, fq, index, outdir, fq2=None, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...

        """

        assert os.path.exists(fq)
        assert AlignIndex(index).is_index('bwa')
        is_path(outdir)

        ## init
        args = args_init(kwargs, align=True) # init

        ## output
        fqname = file_prefix(fq)[0]
        if not fq2 is None:
            fqname = re.sub('_[rR]?1$', '', fqname)
        align_prefix = os.path.join(outdir, fqname)
        
        self.fq_type = Fastx(fq).format()
        self.fq = fq
        self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.fqname = fqname
        self.kwargs = args

        self.sam = align_prefix + '.sam'
        self.bam = align_prefix + '.bam'
        self.log = align_prefix + '.log'
        self.unmap = align_prefix + '.unmap.fq' # fa/fq
        if Fastx(fq).format() == 'fasta':
            self.unmap = align_prefix + '.unmap.fa' # fa/fq

        self.aligner_exe = shutil.which('bowtie2')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. bowtie [options] -U <fq> 
        2. bowtie [options] -1 fq1 -2 fq2
        """
        args = self.kwargs.copy()

        n_map = args.get('n_map', 0)

        ## multi map
        if n_map < 1:
            n_map = 1

        ## unique
        if args['unique_only']:
            arg_unique = '-m 1'
        else:
            arg_unique = '-v 2 -k {}'.format(n_map)

        arg_fx = '-f' if self.fq_type == 'fasta' else '-q'

        cmd = '        '

        return cmd


    def align(self):
        """
        Run bowtie command line
        """
        args = self.kwargs.copy()

        cmd = self.get_cmd()

        if os.path.exists(self.bam) and args['overwrite'] is False:
            log.info('{:>20} : file exists, alignment skipped'.format(self.fqname))
        else:
            run_shell_cmd(cmd)

        ## convert sam to bam
        sam2bam(self.sam, self.bam)

        return self.bam

## not tested
class Hiset2(object):

    def __init__(self, fq, index, outdir, fq2=None, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...

        """

        assert os.path.exists(fq)
        assert AlignIndex(index).is_index('bowtie2')
        is_path(outdir)

        ## init
        args = args_init(kwargs, align=True) # init

        ## output
        fqname = file_prefix(fq)[0]
        if not fq2 is None:
            fqname = re.sub('_[rR]?1$', '', fqname)
        align_prefix = os.path.join(outdir, fqname)
        
        self.fq_type = Fastx(fq).format()
        self.fq = fq
        self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.fqname = fqname
        self.kwargs = args

        self.sam = align_prefix + '.sam'
        self.bam = align_prefix + '.bam'
        self.log = align_prefix + '.log'
        self.unmap = align_prefix + '.unmap.fq' # fa/fq
        if Fastx(fq).format() == 'fasta':
            self.unmap = align_prefix + '.unmap.fa' # fa/fq

        self.aligner_exe = shutil.which('bowtie2')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. bowtie [options] -U <fq> 
        2. bowtie [options] -1 fq1 -2 fq2
        """
        args = self.kwargs.copy()

        n_map = args.get('n_map', 0)

        ## multi map
        if n_map < 1:
            n_map = 1

        ## unique
        if args['unique_only']:
            arg_unique = '-m 1'
        else:
            arg_unique = '-v 2 -k {}'.format(n_map)

        arg_fx = '-f' if self.fq_type == 'fasta' else '-q'

        ## cmd
        cmd = '{} -x {} {} {} -p {} --mm --no-unal --un {}'.format(
            self.aligner_exe,
            self.index,
            arg_fx,
            arg_unique,
            args.get('threads', 4),
            self.unmap)

        ## se or pe
        if self.fq2 is None: 
            cmd += ' {} 1>{} 2>{}'.format(
                self.fq,
                self.sam,
                self.log)
        else:
            cmd += ' -1 {} -2 {} 2>{}'.format(
                self.fq1,
                self.fq2,
                self.log)

        return cmd


    def align(self):
        """
        Run bowtie command line
        """
        args = self.kwargs.copy()

        cmd = self.get_cmd()

        if os.path.exists(self.bam) and args['overwrite'] is False:
            log.info('{:>20} : file exists, alignment skipped'.format(self.fqname))
        else:
            run_shell_cmd(cmd)

        ## convert sam to bam
        sam2bam(self.sam, self.bam)

        return self.bam


## for one index
class AlignIndex(object):

    def __init__(self, aligner='bowtie', index=None, **kwargs):

        args = args_init(kwargs, align=True) # init

        self.aligner = aligner

        if not index is None:
            self.index = index.rstrip('/') # 
            self.name = self.get_name()
            self.aligner_supported = self.get_aligner() # all
            self.check = index if self.is_index() else None
        
        self.kwargs = args        


    def get_aligner(self, index=None):
        """
        Search the available index for aligner:
        bowtie, [*.[1234].ebwt,  *.rev.[12].ebwt]
        bowtie2, [*.[1234].bt2, *.rev.[12].bt2]  
        STAR,
        bwa, 
        hisat2, 
        """
        if index is None:
            index = self.index

        bowtie_files = [index + i for i in [
            '.1.ebwt',
            '.2.ebwt',
            '.3.ebwt',
            '.4.ebwt',
            '.rev.1.ebwt',
            '.rev.2.ebwt']]

        bowtie2_files = [index + i for i in [
            '.1.bt2',
            '.2.bt2',
            '.3.bt2',
            '.4.bt2',
            '.rev.1.bt2',
            '.rev.2.bt2']]

        hisat2_files = ['{}.{}.ht2'.format(index, i) for i in range(1, 9)]


        bwa_files = [index + i for i in [
            '.sa',
            '.amb',
            '.ann',
            '.pac',
            '.bwt']]

        STAR_files = [os.path.join(index, i) for i in [
            'SAindex',
            'Genome',
            'SA',
            'chrLength.txt',
            'chrNameLength.txt',
            'chrName.txt',
            'chrStart.txt',
            'genomeParameters.txt']]

        ## check exists
        bowtie_chk = [os.path.exists(i) for i in bowtie_files]
        bowtie2_chk = [os.path.exists(i) for i in bowtie2_files]
        hisat2_chk = [os.path.exists(i) for i in hisat2_files]
        bwa_chk = [os.path.exists(i) for i in bwa_files]
        STAR_chk = [os.path.exists(i) for i in STAR_files]

        ## check file exists
        aligner = []

        if all(bowtie_chk):
            aligner.append('bowtie')
        elif all(bowtie2_chk):
            aligner.append('bowtie2')
        elif all(hisat2_chk):
            aligner.append('hisat2')
        elif all(bwa_chk):
            aligner.append('bwa')
        elif all(STAR_chk):
            aligner.append('STAR')
        else:
            pass

        return aligner


    def is_index(self, index=None):
        """
        Check if index support for aligner
        """
        if index is None:
            index = self.index
        return self.aligner in self.get_aligner(index)

        # return self.aligner in self.aligner_supported if 
        #    index is None else self.aligner in self.get_aligner(index)


    def get_name(self, index=None):
        """
        Get the name of index
        basename: bowtie, bowtie2, hisqt2, bwa
        folder: STAR
        """
        if index is None:
            index = self.index
        if os.path.isdir(index):
            # STAR
            iname = os.path.basename(index)
        else:
            # bowtie, bowtie2, bwa, hisat2
            iname = os.path.basename(index)

        return iname


    def search(self, genome=None, group='genome'):
        """
        Search the index for aligner: STAR, bowtie, bowtie2, bwa, hisat2
        para:

        *genome*    The ucsc name of the genome, dm3, dm6, mm9, mm10, hg19, hg38, ...
        *type*      Choose from: genome, rRNA, transposon, piRNA_cluster, ...

        structure of genome_path:
        default: {HOME}/data/genome/{genome_version}/{aligner}/

        path-to-genome/
            |- Bowtie_index /
                |- genome
                |- rRNA
                |- MT_trRNA
            |- transposon  
            |- piRNA cluster

        """
        args = self.kwargs.copy()

        g_path = args.get('genome_path', './')
        i_prefix = os.path.join(g_path, genome, self.aligner + '_index')

        # all index
        d = {
            'genome': [
                os.path.join(i_prefix, 'genome'),
                os.path.join(i_prefix, genome)],
            'genome_rm': [
                os.path.join(i_prefix, 'genome_rm'),
                os.path.join(i_prefix, genome + '_rm')],
            'rRNA': [
                os.path.join(i_prefix, 'rRNA'),
                os.path.join(i_prefix, 'MT_trRNA')],
            'te': [
                os.path.join(i_prefix, 'transposon')],
            'piRNA_cluster': [
                os.path.join(i_prefix, 'piRNA_cluster')],
            'miRNA': [
                os.path.join(i_prefix, 'miRNA')],
            'miRNA_hairpin': [
                os.path.join(i_prefix, 'miRNA_hairpin')],
            'structureRNA': [
                os.path.join(i_prefix, 'structureRNA')]}

        # hit
        i_list = d.get(group, ['temp_temp'])
        i_list = [i for i in i_list if self.is_index(i)]

        return i_list[0] if len(i_list) > 0 else None

## for one fq, all index
class AlignIndexBuilder(object):
    """
    Search index :
        for alignment
        search index, parameters from command/arguments/parameters
        extra_index (priority=1)
        genome
        spikein
        to_rRNA
        to_te
        to_piRNA_cluster
        ...
    : arguments :
    1. genomme (str), spikein (str|None), align_to_rRNA (True|False)
    2. align_to_te (True|False), te_index (str|None) (mapping)
    3. extra_index (str|None)
    """
    def __init__(self, aligner='bowtie', **kwargs):

        ## args init
        args = args_init(kwargs, align=True)

        ## update arguments
        self.genome = args.get('genome', 'dm3')
        self.spikein = args.get('spikein', 'dm3')
        self.align_to_rRNA = args.get('align_to_rRNA', True)
        self.extra_index = args.get('extra_index', None)
        self.aligner = aligner
        self.kwargs = args
        self.index_list = self.builder()


    def builder(self):
        """
        Return the list of index
        extra
        genome (rRNA) 
        spikein (rRNA)
        ...
        """
        args = self.kwargs.copy()

        ## extra
        ##   : list of index, str
        if not self.extra_index is None:
            if isinstance(self.extra_index, str):
                self.extra_index = [self.extra_index]
            elif isinstance(self.extra_index, list):
                pass
            else:
                raise Exception('{:>10} : --extra-index, expect str or list, \
                    get {}'.format('error', type(self.extra_index)))
            index_list = self.extra_index
        elif isinstance(self.genome, str): # Genome(self.genome).check
            index_list = []
            ## genome
            if self.align_to_rRNA:
                index_list.append(
                    AlignIndex(self.aligner).search(self.genome, 'rRNA'))
            index_list.append(AlignIndex(self.aligner).search(self.genome, 'genome'))
            ## spikein
            if isinstance(self.spikein, str):
                if self.align_to_rRNA:
                    index_list.append(
                        AlignIndex(self.aligner).search(self.spikein, 'rRNA'))
                index_list.append(
                    AlignIndex(self.aligner).search(self.spikein, 'genome'))

        else:
            raise Exception('{:>10} : --genome and --extra-index not valid'.format('error'))

        ## check
        return [i for i in index_list if AlignIndex(self.aligner, i).check]


## For one fq, one index
class AlignConfig(object):

    def __init__(self, aligner, fq, index, outdir, fq2=None, outdir_fixed=False, **kwargs):
        """
        Config for alignment, for one fq, one index
        directories
        files
        ...
        """
        assert isinstance(fq, str)
        assert isinstance(index, str)
        assert os.path.exists(fq)
        assert AlignIndex(aligner, index).is_index()

        ## init
        args = args_init(kwargs, align=True) # init

        ## file format        
        self.fx_type = Fastx(fq).format
        self.format = self.fx_type
        self.aligner = aligner
        self.fq = self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.outdir_fixed = outdir_fixed
        self.kwargs = args
        tmp = self.config() # update self
        self.check = self.check_exists()


    def check_exists(self):
        """
        Check arguments, target file exists, 
        """
        args = self.kwargs.copy()

        ## save parameters
        ## update
        args['aligner'] = self.aligner
        args['fq'] = args['fq1'] = self.fq
        args['fq2'] = self.fq2
        args['index'] = self.index
        args['outdir'] = self.outdir
        args['outdir_fixed'] = self.outdir_fixed

        ## chk files
        args_pickle = self.out_prefix + '.arguments.pickle' # out_prefix
        args_file = self.out_prefix + '.arguments.txt'

        ## chk
        chk1 = args_checker(args, args_pickle)
        chk2 = args['overwrite'] is False
        chk3 = os.path.exists(self.bam)

        args_logger(args, args_file, True) # update arguments.txt

        return all([chk2, chk3]) # arguments.pickle, changed 


    def config(self):
        """
        Create directories, filenames
        *fq*     The fastq file (or fq1 of PE)
        *index*  The aligner index file
        *args*   Align_path, output directory of the alignment    

        return files:
        prefix, bam, bed, log, unmap
        """
        args = self.kwargs.copy()

        ## fqname
        fqname = file_prefix(self.fq)[0]

        if not self.fq2 is None:
            fqname = re.sub('[._][rR]?1$', '', fqname)

        ## remove 'unmap' suffix
        fqname = re.sub('.unmap$', '', fqname) #!!!

        if args['simplify_name']:
            fqname = re.sub('.not_\w+|.map_\w+', '', fqname)

        if args['smp_name']:
            fqname = args['smp_name']

        ## index name
        index_name = args.get('index_name', None)
        if index_name is None:
            index_name = AlignIndex(self.aligner, self.index).name

        ## output (fixed?!)
        if self.outdir_fixed:
            sub_outdir = self.outdir
        else:
            sub_outdir = os.path.join(self.outdir, fqname, index_name)

        out_prefix = os.path.join(sub_outdir, fqname)
        self.sam = out_prefix + '.sam'
        self.bam = out_prefix + '.bam'
        self.log = out_prefix + '.log'
        self.unmap = out_prefix + '.unmap.' + self.fx_type # default

        if self.fq2 is None:
            self.unmap1 = self.unmap2 = None
        else:
            if self.aligner == 'bowtie':
                self.unmap1 = out_prefix + '.unmap_1.' + self.fx_type
                self.unmap2 = out_prefix + '.unmap_2.' + self.fx_type
            else:
                self.unmap1 = out_prefix + '.unmap.1.' + self.fx_type
                self.unmap2 = out_prefix + '.unmap.2.' + self.fx_type


        self.fqname = fqname
        self.out_prefix = out_prefix        
        self.index_name = index_name
        self.sub_outdir = sub_outdir

        ## create directory
        is_path(sub_outdir)

        # print(self.fq, fqname, self.bam)

        return fqname, out_prefix, self.sam, self.bam, self.log, self.unmap, self.unmap1, self.unmap2


class AlignNFastx(object):

    def __init__(self, aligner, fq, outdir, fq2=None, index_list=None, **kwargs):
        """
        Top level for Alignment
        for N fastq and N index
        """

        ## fastq files
        self.fq = self.fq1 = fq
        self.fq2 = fq2

        if isinstance(fq, str):
            self.fq = self.fq1 = [fq]
            self.fq2 = [fq2]

        if not self.fq2 is None:
            if not len(self.fq) == len(self.fq2):
                raise Exception('{:>10} : fq1 and fq2 not paired for PE reads'.format('error'))

        ## init
        args = args_init(kwargs, align=True)

        ## update args
        args.pop('aligner', None)

        ## get index_list
        if index_list is None:
            index_list = AlignIndexBuilder(aligner, **args).index_list

        # update args
        args.pop('aligner', None)
        args.pop('outdir', None)
        args.pop('index_list', None)
        args.pop('fq', None)
        self.fq = self.fq1 = fq
        args['fq2'] = fq2

        self.aligner = aligner
        self.outdir = outdir
        self.index_list = index_list
        self.kwargs = args


    def run(self):
        args = self.kwargs.copy()

        bam_list = []
        for fq1, fq2 in zip(self.fq, self.fq2):
            # print('AAA')
            fq_bams = AlignNIndex(self.aligner, self.fq, self.outdir, self.index_list, **args).run()
            bam_list.append(fq_bams)

        return bam_list


class AlignNIndex(object):
    """
    Align fq/fa to multiple reference
      - in_parallel
      - sequential
    """
    def __init__(self, aligner, fq, outdir, index_list=None, fq2=None, **kwargs):
        """
        Get multiple index 
        """
        args = args_init(kwargs, align=True)

        ## prepare list
        if index_list is None:
            index_list = AlignIndexBuilder(aligner, **args).index_list

        if len(index_list) == 0:
            raise Exception('{:>10} : index not detected'.format('error'))

        # remove dup keys
        args.pop('aligner', None)
        args.pop('fq', None)
        args.pop('index', None)
        args.pop('outdir', None)
        args.pop('fq2', None)

        self.aligner = aligner
        self.fq = fq
        self.fq2 = fq2
        self.outdir = outdir
        self.index_list = index_list
        self.kwargs = args
        self.index_parallel = args.get('index_parallel', False) # for multiple index

    def run(self):
        """
        Run Alignment() one by one
        """
        args = self.kwargs.copy()

        bam_list = []
        for i, index in enumerate(self.index_list):
            ## custom para
            args.pop('index_name', None) # pop
            index_name = str(i + 1) + '.' + AlignIndex(self.aligner, index).name
            config = AlignConfig(
                self.aligner, 
                self.fq, 
                index, 
                self.outdir,
                fq2 = self.fq2,
                index_name=index_name,
                **args)

            print(config.unmap, config.unmap1, config.unmap2)

            bam = AlignOneIndex(
                self.aligner, 
                self.fq, 
                index, 
                config.sub_outdir, 
                self.fq2,
                index_name=index_name, 
                **args).run()
            bam_list.append(bam)

            ## update fastq file!? parallel
            if not self.index_parallel:
                if self.fq2 is None:
                    self.fq = config.unmap
                else:
                    self.fq = self.fq1 = config.unmap1
                    self.fq2 = config.unmap2

        return bam_list


class AlignOneIndex(object):
    """
    Align fa/fq to reference, 1 to 1
    """
    def __init__(self, aligner, fq, index, outdir, fq2=None, **kwargs):
        assert isinstance(index, str) # 1 index

        ## args
        args = args_init(kwargs, align=True)

        ## update args
        args['aligner'] = aligner
        args['fq'] = args['fq1'] = fq
        args['fq2'] = fq2
        args['index'] = index
        args['outdir'] = outdir

        self.aligner = aligner
        self.args = args
        self.fq = self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.kwargs = args


    def get_aligner_port(self):
        """
        Determine the class object for aligner
        """
        args = self.kwargs.copy()

        aligner = self.aligner.lower()

        if aligner == 'bowtie':
            port = Bowtie # class
        elif aligner == 'bowtie2':
            port = Bowtie2 # class
        elif aligner == 'star':
            port = Star # class
        elif aligner == 'bwa':
            port = Bwa # class
        elif aligner == 'hisat2':
            port = Hisat2 # class
        else:
            raise Exception('{:>10} : aligner not supported {}'.format(self.aligner))

        return port


    def run(self):
        args = self.kwargs.copy()

        aligner_port = self.get_aligner_port()

        return aligner_port(**args).run()


class Bowtie(object):

    def __init__(self, fq, index, outdir, fq2=None, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...
        """
        assert os.path.exists(fq)
        args = args_init(kwargs, align=True) # init

        ## update
        args.pop('aligner', None)
        args.pop('fq', None)
        args.pop('fq2', None)
        args.pop('index', None)
        args.pop('outdir', None)

        ## check
        ## saving common configs for alignment
        ## in Class, ignore outdir
        self.config = AlignConfig('bowtie', fq, index, outdir, fq2, outdir_fixed=True, **args)
        self.fq = self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.kwargs = args
        self.aligner_exe = shutil.which('bowtie')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. bowtie [options] -U <fq> 
        2. bowtie [options] -1 fq1 -2 fq2
        """
        args = self.kwargs.copy()

        n_map = args.get('n_map', 0)

        ## multi map
        if n_map < 1:
            n_map = 1

        ## unique
        if args['unique_only']:
            arg_unique = '-m 1'
        else:
            arg_unique = '-v 2 -k {}'.format(n_map)

        arg_fx = '-f' if self.config.format == 'fasta' else '-q'

        ## cmd
        cmd = '{} {} {} {} -p {} --mm --best --sam --no-unal --un {}'.format(
            self.aligner_exe,
            self.index,
            arg_fx,
            arg_unique,
            args.get('threads', 4),
            self.config.unmap)

        ## se or pe
        if self.fq2 is None: 
            cmd += ' {} 1>{} 2>{}'.format(
                self.fq,
                self.config.sam,
                self.config.log)
        else:
            cmd += ' -1 {} -2 {} 1> {} 2>{}'.format(
                self.fq1,
                self.fq2,
                self.config.sam,
                self.config.log)

        return cmd


    def run(self):
        """
        Run bowtie command line
        """
        args = self.kwargs.copy()
        cmd = self.get_cmd() # output

        ## para, bam, overwrite
        if self.config.check:
            log.info('{:>20} : file exists, alignment skipped'.format(self.config.fqname))
        else:
            run_shell_cmd(cmd)
            sam2bam(self.config.sam, self.config.bam, sort=True) # convert to bam

        return self.config.bam
        

class Bowtie2(object):

    def __init__(self, fq, index, outdir, fq2=None, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...
        """
        assert os.path.exists(fq)
        args = args_init(kwargs, align=True) # init

        ## update
        args.pop('aligner', None)
        args.pop('fq', None)
        args.pop('fq2', None)
        args.pop('index', None)
        args.pop('outdir', None)

        ## check
        ## saving common configs for alignment
        ## in Class, ignore outdir
        self.config = AlignConfig('bowtie2', fq, index, outdir, fq2, outdir_fixed=True, **args)
        self.fq = self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.kwargs = args
        self.aligner_exe = shutil.which('bowtie2')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. bowtie [options] -U <fq> 
        2. bowtie [options] -1 fq1 -2 fq2
        """
        args = self.kwargs.copy()

        ## unique
        if args['unique_only']:
            arg_unique = '-q 30' # samtools filtering
        else:
            arg_unique = ''

        self.arg_unique = arg_unique # pass to other func

        ## multi map
        n_map = args.get('n_map', 1)
        if n_map < 2:
            # n_map = 1 # default 1, report 1 hit for each read
            arg_multi = '' # default: 1
        else:
            arg_multi = '-k {}'.format(n_map)

        ## fx type
        arg_fx = '-f' if self.config.format == 'fasta' else '-q'

        ## cmd
        cmd = '{} -x {} {} {} -p {} --very-sensitive-local --mm --no-unal'.format(
            self.aligner_exe,
            self.index,
            arg_fx,
            arg_unique,
            args.get('threads', 4))

        ## se or pe
        if self.fq2 is None: 
            cmd += ' --un {} {} 1>{} 2>{}'.format(
                self.config.unmap,
                self.fq,
                self.config.sam,
                self.config.log)
        else:
            cmd += ' --un-conc {} -1 {} -2 {} 1>{} 2>{}'.format(
                self.config.unmap,
                self.fq1,
                self.fq2,
                self.config.sam,
                self.config.log)

        return cmd


    def run(self):
        """
        Run bowtie command line
        """
        args = self.kwargs.copy()
        cmd = self.get_cmd() # output

        ## para, bam, overwrite
        if self.config.check:
            log.info('{:>20} : file exists, alignment skipped'.format(self.config.fqname))
        else:
            run_shell_cmd(cmd)
            sam2bam(self.config.sam, self.config.bam, sort=True, extra_para=self.arg_unique) # convert to bam

        return self.config.bam
        

class Star(object):

    def __init__(self, fq, index, outdir, fq2=None, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...

        use '--outFilterMultimapNmax' to control uniquely mapped reads

        # update 2019-03-21
        since 99% of the reads do not map to Blumeria, STAR takes a lot of time trying to squeeze the reads into the wrong genome. The best solution is to make a combined genome of Barley and Blumeria, which will alloy mapping simultaneously to the two genomes, with the best alignments winning.
        Another option (if you insist on mapping to Blumeria alone) is to reduce --seedPerWindowNmax to 20 or even smaller values. More discussion on it here: https://groups.google.com/d/msg/rna-star/hJL_DUtliCY/HtpiePlMBtYJ .
        see: https://github.com/alexdobin/STAR/issues/329#issuecomment-334879474
        """
        assert os.path.exists(fq)
        args = args_init(kwargs, align=True) # init

        ## update
        args.pop('aligner', None)
        args.pop('fq', None)
        args.pop('fq2', None)
        args.pop('index', None)
        args.pop('outdir', None)

        ## check
        ## saving common configs for alignment
        ## in Class, ignore outdir
        self.config = AlignConfig('STAR', fq, index, outdir, fq2, outdir_fixed=True, **args)
        self.fq = self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.kwargs = args
        self.aligner_exe = shutil.which('STAR')

    def get_cmd(self):
        """
        Create basic alignment command line
        1. STAR
        """
        args = self.kwargs.copy()

        ## for small genome
        ## mapping takes too long,
        ## 99% of the reads not mapped
        ## change --seedPerWindowNmax 
        ## max number of seeds per window
        ## https://github.com/alexdobin/STAR/issues/329
        ## https://groups.google.com/d/msg/rna-star/hJL_DUtliCY/HtpiePlMBtYJ 
        if args['small_genome']:
            seed_max = 5 # even smaller
        else:
            seed_max = 50 # default

        ## for unique map
        ## --outFilterMultimapNmax
        ## maximum number of loci the read is allowed to map to. [default: 10]
        ##
        ## filt unique reads by: samtools view -q 30
        n_map = args.get('n_map', 10)
        if n_map == 0:
            n_map = 10 # default
        if args['unique_only']:
            n_map = 1 # unique # not working
        arg_unique = '--outFilterMultimapNmax {} --seedPerWindowNmax {}'.format(n_map, seed_max)
        
        ## file type
        arg_reader = 'zcat' if self.fq.endswith('.gz') else '-'

        ## input fq
        # self.fq2 = '' if self.fq2 is None else self.fq2
        if self.fq2 is None:
            self.fq2 = ''

        # cmd
        ## output prefix
        out_prefix = self.config.out_prefix
        cmd = '{} --runMode alignReads \
            --genomeDir {} \
            --readFilesIn {} {} \
            --readFilesCommand {} \
            --outFileNamePrefix {} \
            --runThreadN {} \
            --limitOutSAMoneReadBytes 1000000 \
            --genomeLoad NoSharedMemory  \
            --limitBAMsortRAM 10000000000 \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMismatchNoverLmax 0.07 \
            --seedSearchStartLmax 20 \
            --outReadsUnmapped Fastx {} {}'.format(
                self.aligner_exe,
                self.index,
                self.fq1,
                self.fq2,
                arg_reader,
                self.config.out_prefix,
                args.get('threads', 4),
                self.config.unmap,
                arg_unique)

        return cmd


    def update_names(self, keep_old=False):
        """
        Update the filenames of STAR output
        bam: *Aligned.sortedByCoord.out.bam -> *.bam # mv
        log: *Log.final.out -> *.log # copy
        log: *Log.out -> *.out
        unmap: *Unmapped.out.mate1 -> *.unmap.1.fastq
               *Unmapped.out.mate1 -> *.unmap.1.fastq
        """
        args = self.kwargs.copy()

        bam_from = self.config.out_prefix + 'Aligned.sortedByCoord.out.bam'
        log_from = self.config.out_prefix + 'Log.final.out'
        unmap = unmap1 = self.config.out_prefix + 'Unmapped.out.mate1'
        unmap2 = self.config.out_prefix + 'Unmapped.out.mate2'

        if not keep_old:
            ## rename files
            ## move BAM, unmap
            ## copy log
            flag = 0
            if os.path.exists(bam_from) and not os.path.exists(self.config.bam):
                # shutil.copy(bam_from, self.config.bam) # original, contains unique + multiple
                os.symlink(os.path.basename(bam_from), self.config.bam) # symlink
                flag += 1
            if os.path.exists(log_from):
                shutil.copy(log_from, self.config.log)
                flag += 1
            # print('1', unmap, self.config.unmap)
            # print(self.fq2)
            ## why self.fq2 changed, from None to '' !!!
            if not self.fq2:
                if os.path.exists(unmap):
                    # print('2', unmap, self.config.unmap)
                    shutil.move(unmap, self.config.unmap)
                    flag += 1
            else:
                if os.path.exists(unmap1):
                    # print('3', unmap1, self.config.unmap1)
                    shutil.move(unmap1, self.config.unmap1)
                    flag += 1
                if os.path.exists(unmap2):
                    shutil.move(unmap2, self.config.unmap2)

            return self.config.bam
        else:
            return bam_from



    def run(self):
        """
        Run STAR command line
        """
        args = self.kwargs.copy()
        cmd = self.get_cmd() # output

        ## para, bam, overwrite
        if self.config.check:
            log.info('{:>20} : file exists, alignment skipped'.format(self.config.fqname))
        else:
            run_shell_cmd(cmd)
            ## unique
            ## -q 30
            bam_from = self.update_names(False)
            if args['unique_only']:
                sam2bam(bam_from, self.config.bam, sort=True, extra_para='-q 30')
                os.remove(bam_from)
            self.update_names(True)

        return self.config.bam


