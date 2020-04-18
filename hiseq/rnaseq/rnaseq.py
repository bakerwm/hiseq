#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Working mode:

1. build design

Create RNAseq_auto_design.txt for further analysis

2. pickle

Pass all arguments from *.pickle file

3. design

pass values from design.txt



## Design

input:

- control (rep1, rep2, ...)
- treatment (rep1, rep2, ...)
- genome

output:

  - single:
    - raw_data
    - clean_data
    - align
    - count
    - qc
    - report
    ...

  - pair:
    - count
    - deseq
    - qc
    - report


## figures/tables
1. scatter plot
2. MA plot
3. volcano
4. heatplot
...

## subfunctions:

class RNAseqSingle()
class RNAseqPair()
class RNAseqConfig()
...
object to dict, ...


DE analysis:

DEseq2, edgeR, ...



Alignment:

gene,
TE
piRNA cluster
...

"""

import os
import re
import hiseq
import pysam
import shutil
import tempfile
import pandas as pd
from hiseq.utils.helper import *
from hiseq.qc.trimmer import Trimmer
from hiseq.align.alignment import Alignment


def print_df(d):
    if isinstance(d, dict):
        for k, v in d.items():
            print('{:>15} : {}'.format(k, v))
    else:
        print(d)


class RNAseqReader(object):
    """
    Return the config, files for RNAseq direcotory

    RNAseqSingle: config/
    RNAseqMultiple: config/
    DeseqSingle: config/

    # to-do
    DEseqMultiple: ? (RNAseqMultiple)
    """
    def __init__(self, path, **kwargs):
        self.path = path
        self.feature = kwargs.get('feature', 'gene')
        # self.return_config = kwargs.get('return_args', False)
        # self.rnaseq_type, self.args = self.check_rnaseq_type(return_args=True)


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def check(self, return_args=False):
        """
        Check if the directory is of RNAseq 
        [single|multiple|deseq_single|deseq_multiple]
        single - config:
            outdir/feature/config/arguments.pickle
        multiple - config:
            outdir/config/feature/config/arguments.pickle
        deseq - single:
            outdir/feature/config/arguments.pickle
        """
        x1 = os.path.join(self.path, self.feature, 'config', 'arguments.pickle')
        x2 = os.path.join(self.path, 'config', self.feature, 'config', 'arguments.pickle')

        if os.path.exists(x1):
            args = pickle_to_dict(x1)
        elif os.path.exists(x2):
            args = pickle_to_dict(x2)
        else:
            args = {}
            log.error("""
                unknown directory, expect config file:
                RNAseq single: {}
                RNAseq multiple: {}
                RNAseq deseq: {}""".format(x1, x2, x1))
        
        # update self
        self.update(args) # fresh new

        if return_args:
            return (tag, args)
        else:
            return tag


    def is_rnaseq_single(self):
        """
        Check if the directory is of RNAseq single sample
        outdir/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'rnaseq_single'


    def is_rnaseq_multi(self):
        """
        Check if the directory is of RNAseq multiple sample
        outdir/config/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'rnaseq_multiple'


    def is_deseq_single(self):
        """
        Check if the directory if of RNAseq DESeq analysis
        outdir/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'deseq_single'


    def is_deseq_multiple(self):
        """
        Check if the directory if of RNAseq DESeq analysis
        outdir/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'deseq_multiple'


class RNAseqConfig(object):
    """
    UPDATED:
    directory type: single, pair
    config: init_atac, single, merge, multiple
    list all files (objects)

    required args:
    fq1, fq2, genome, outdir, ...

    options:
    """
    def __init__(self, **kwargs):
        self.update(kwargs) # update from args
        self.init_rnaseq() # prepare args
        self.rnaseq_type = self.mission_type()
        self.init_subgroup()


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def init_rnaseq(self):
        """
        check arguments conflicts, defaults...
        """
        # default values
        ## critical args:
        ## PE mode get low pct unique mapped reads, but SE mode not.
        ## so force to SE mode
        args_default = {
            'read1_only': False,
            'build_design': False,
            'pickle': None,
            'design': None,
            'genome': 'mm10',
            'outdir': str(pathlib.Path.cwd()),
            'feature': 'gene',
            'fq1': None, 
            'fq2': None,
            'dirs_ctl': None,
            'dirs_exp': None,
            'smp_path': None,
            'smp_name': None,
            'group': None,
            'gtf': None,
            'overwrite': False
        }
        self.update(args_default, overwrite=False) # update missing attrs
        # 1st level: build design

        # 2nd level: pickle / update all config
        if not self.pickle is None:
            if os.path.exists(self.pickle) and self.pickle.endswith('.pickle'):
                args_pickle = pickle_to_dict(self.pickle)
                ## clear all attributes
                for k in self.__dict__:
                    self.__delattr__(k)
                # assign attrs from pickle 
                self.update(args_pickle) # update all
            else:
                raise Exception('--pickle, failed: {}'.format(self.pickle))

        # 3rd level: design (input)
        # update args from design.txt
        # 1. group, smp_name, feature, genome, outdir, fq1, fq2
        # 2. group, smp_name, feature, genome, outdir, smp_path
        if not self.design is None:
            args_design = DesignReader(design).to_dict()
            self.update(args_design) # update specific args


        # 4th level: read1 only
        if self.read1_only is True:
            self.fq2 = None

        # 5th level: default path, files

        ## smp_name [from fq1, smp_path]
        if self.smp_name is None:
            if self.smp_path:
                self.smp_name = fq_name(self.smp_path)
            elif self.fq1:
                self.smp_name = fq_name(self.fq1)
            else:
                pass # None

        ## group
        if self.group is None:
            if self.smp_name:
                self.group = fq_name_rmrep(self.smp_name)
            else:
                pass # None

        ## outdir
        self.projectdir = os.path.join(self.outdir, self.feature)
        self.configdir = os.path.join(self.projectdir, 'config')

        ## config
        self.auto_design = os.path.join(self.configdir, 'RNAseq_auto_design.txt') # new created
        self.config_txt = os.path.join(self.configdir, 'arguments.txt')
        self.config_pickle = os.path.join(self.configdir, 'arguments.pickle')
        self.config_json = os.path.join(self.configdir, 'arguments.json')


    def mission_type(self):
        """
        Determine the purpose the RNAseq analysis
        1. single
        2. pair
        ...
        """
        # print('!CCCC')
        # print_df(self.__dict__)
        if self.build_design is True:
            flag = 'build_design'
        # elif not self.pickle is None:
        #    flag = 'from_pickle'
        elif isinstance(self.smp_path, list):
            flag = 'deseq_multiple'
        elif isinstance(self.dirs_ctl, list) and isinstance(self.dirs_exp, list):
            flag = 'deseq_single'
        elif isinstance(self.fq1, list):
            flag = 'rnaseq_multiple'
        elif isinstance(self.fq1, str):
            flag = 'rnaseq_single'
        else:
            raise Exception("""unknown RNAseq() arguments:
                RNAseq-single:   fq1, fq2, genome, outdir;
                RNAseq-multiple: [fq1], [fq2], genome, outdir | design.txt
                DESeq-single:    dirs_ctl, dirs_exp;
                DESeq-multiple:  smp_path
                """)

        return flag


    def init_subgroup(self):
        # if self.rnaseq_type == 'from_pickle':
        #     log.info('read args from pickle')
        #    args_pickle = pickle_to_dict(self.pickle)
        #    for k, v in args.items():
        #        setattr(self, k, v) # update all
        #elif self.rnaseq_type == 'build_design':
        # create_dirs = self.args.get('create_dirs', True)
        create_dirs = getattr(self, 'create_dirs', True)


        if self.rnaseq_type == 'build_design':
            self.init_build_design(create_dirs)
        elif self.rnaseq_type == 'rnaseq_single':
            self.init_rnaseq_single(create_dirs)
        elif self.rnaseq_type == 'rnaseq_multiple':
            self.init_rnaseq_multiple(create_dirs)
        elif self.rnaseq_type == 'deseq_single':
            self.init_deseq_single(create_dirs)
        elif self.rnaseq_type == 'deseq_multiple':
            self.init_deseq_multiple(create_dirs)
        else:
            log.error('unknown rnaseq type: {}'.format(self.rnaseq_type))
            pass


    def init_build_design(self, create_dirs=True):
        """
        Create design.txt for project
        required:
            fq1, (fq2)
            genome
            outdir
            feature
            gtf

        rnaseq:
        names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
        rnaseq:
        names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']

        atacseq:
        names=['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']
        """
        if create_dirs is True:
            check_path([self.outdir, self.configdir])


    def init_rnaseq_single(self, create_dirs=True):
        """
        Initiate the config, directories, files for rnaseq single
        update self.
        fq1, fq2, genome, outdir: args,
        or:
        design
        """
        # args = self.args.copy() # global
        # self.fq1 = os.path.abspath(args['fq1'])
        # self.fq2 = os.path.abspath(args['fq2']) if args['fq2'] else args['fq2']
        # if self.read1_only: self.fq2 = None  # read1_only
        # self.genome = args['genome']
        # self.outdir = os.path.abspath(args['outdir'])
        # self.feature = args.get('feature', 'gene')

        chk1 = isinstance(self.fq1, str) # checked in: mission_type

        ## sample name
        if isinstance(self.smp_name , list):
            self.smp_name = self.smp_name
        # if self.smp_name is None:
        #     self.smp_name = fq_name(self.fq1)
        # self.fqname= self.smp_name

        ## genome GTF
        ## optional: ucsc, ensembl
        if self.gtf is None:
            self.gtf = Genome(genome=self.genome).gene_gtf('ensembl')

        ## paths
        self.rawdir = os.path.join(self.projectdir, 'raw_data')
        self.cleandir = os.path.join(self.projectdir, 'clean_data')
        self.aligndir = os.path.join(self.projectdir, 'align')
        self.bamdir = os.path.join(self.projectdir, 'bam_files')
        self.bwdir = os.path.join(self.projectdir, 'bw_files')
        self.countdir = os.path.join(self.projectdir, 'count')
        self.reportdir = os.path.join(self.projectdir, 'report')
        self.out_prefix = os.path.join(self.projectdir, self.smp_name)

        ## files
        ## raw data
        self.raw_fq_list = [os.path.join(self.rawdir, os.path.basename(self.fq1))]
        fq2_raw = None if self.fq2 is None else os.path.join(self.rawdir, os.path.basename(self.fq2))
        self.raw_fq_list.append(fq2_raw)

        ## clean data
        self.clean_fq_list = [os.path.join(self.cleandir, fq_name(self.fq1) + '.fq.gz')]
        fq2_clean = None if self.fq2 is None else os.path.join(self.cleandir, fq_name(self.fq2) + '.fq.gz')
        self.clean_fq_list.append(fq2_clean)

        ## files
        self.trim_stat = os.path.join(self.cleandir, self.smp_name + '.qc.stat')
        self.bam_raw = os.path.join(self.aligndir, self.smp_name, '2.*', self.smp_name + '.bam')
        self.align_stat = os.path.join(self.aligndir, self.smp_name + '.align.txt')
        self.bw_fwd = os.path.join(self.bwdir, self.smp_name + '.fwd.bigWig')
        self.bw_rev = os.path.join(self.bwdir, self.smp_name + '.rev.bigWig')
        self.count_sens = os.path.join(self.countdir, 'count.sens.txt')
        self.count_anti = os.path.join(self.countdir, 'count.anti.txt')
        self.strandness_status = os.path.join(self.countdir, 'strandness_status.out')

        ## create directories
        if create_dirs is True:
            check_path([
                self.configdir,
                self.rawdir,
                self.cleandir,
                self.aligndir,
                self.bamdir,
                self.bwdir,
                self.countdir,
                self.reportdir])


    def init_rnaseq_multiple(self, create_dirs=True):
        """
        Initiate the config, for multiple RNAseq samples
        """
        chk1 = isinstance(self.fq1, list) # make sure

        ## fq files
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = None if self.fq2 is None else file_abspath(self.fq2)

        ## sample name
        chkb0 = isinstance(self.smp_name, list) # list
        chkb1 = len(self.smp_name) == len(self.fq1) # number of samples, names
        chkb2 = len(self.smp_name) == len(set(self.smp_name)) # unique names
        if not all([chkb0, chkb1, chkb2]):
            self.smp_name = fq_name(self.fq1) # auto-genrrated 

        ## get GTF file
        if self.gtf is None:
            self.gtf = Genome(genome=self.genome).gene_gtf('ensembl') # ucsc version

        # paths
        self.project_outdir = os.path.join(self.outdir, 'config', self.feature)
        self.configdir = os.path.join(self.project_outdir, 'config')
        # project files
        self.config_txt = os.path.join(self.configdir, 'arguments.txt')
        self.config_pickle = os.path.join(self.configdir, 'arguments.pickle')
        self.config_json = os.path.join(self.configdir, 'arguments.json')

        # create
        if create_dirs is True:
            check_path([self.configdir])


    def init_deseq_single(self, create_dirs=True):
        """
        1-vs-1:
        DE analysis for two group of samples
        dirs_ctl:
        dirs_exp:
        genome:
        outdir:
        feature:
        """
        # is list
        chk1 = isinstance(self.dirs_ctl, list)
        chk2 = isinstance(self.dirs_exp, list)
        # is RNAseq single
        chk3 = all(RNAseqReader(i).is_rnaseq_single() for i in self.dirs_ctl + self.dirs_exp)

        # smp_name
        chkb0 = isinstance(self.smp_name, list)  # is list
        if not chkb0: self.smp_name = fq_name(self.dirs_ctl, dirs_exp)
        chkb1 = len(self.smp_name) == len(self.dirs_ctl + self.dirs_exp) # length
        chkb2 = len(self.smp_name) == len(set(self.smp_name)) # unique
        if not all([chkb1, chkb2]):
            # auto-generated
            self.smp_name = [RNAseqReader(i).smp_name for i in self.dirs_ctl + self.dirs_exp]

        # group
        chkc0 = isinstance(self.group, list) # is list
        if not chkc0: self.group = fq_name_rmrep(self.smp_name)
        chkc1 = len(self.group) == len(self.dirs_ctl + self.dirs_exp) # length
        chkc2 = len(set(self.group)) == 2 # paired
        chkc3 = len(set(self.group[:len(self.dirs_ctl)])) == 1 # ctl
        chkc4 = len(set(self.group[len(self.dirs_ctl):])) == 1 # exp
        if not all([chkc1, chkc2, chkc3, chkc4]):
            # auto-generated
            self.group = fq_name_rmrep(self.smp_name)

        ## paths
        self.prefix_ctl, self.prefix_exp = list_uniquer(self.group, sorted=False)
        self.project_name = '{}.vs.{}'.format(self.prefix_ctl, self.prefix_exp)
        self.project_outdir = os.path.join(self.outdir, self.project_name, self.feature)
        # project dir
        self.configdir = os.path.join(self.project_outdir, 'config')
        self.countdir = os.path.join(self.project_outdir, 'count')
        self.deseqdir = os.path.join(self.project_outdir, 'deseq')
        self.enrichdir = os.path.join(self.project_outdir, 'enrich')
        self.pdfdir = os.path.join(self.project_outdir, 'pdf')
        self.reportdir = os.path.join(self.project_outdir, 'report')
        # project files
        self.config_txt = os.path.join(self.configdir, 'arguments.txt')
        self.config_pickle = os.path.join(self.configdir, 'arguments.pickle')
        self.config_json = os.path.join(self.configdir, 'arguments.json')
        self.deseq_design = os.path.join(self.configdir, 'deseq_design.txt')
        # list of count
        self.count_ctl = [RNAseqReader(i).count_sens for i in self.dirs_ctl]
        self.count_exp = [RNAseqReader(i).count_sens for i in self.dirs_exp]

        ## create directories
        if create_dirs is True:
            check_path([
                self.configdir,
                self.countdir,
                self.deseqdir,
                self.enrichdir,
                self.reportdir])


    def init_deseq_multiple(self, create_dirs=True):
        """
        N: pairs deseq
        1-vs-1:
        DE analysis for N groups of samples
        smp_path:
        group:
        genome:
        outdir:
        feature:
        """
        pass


class RNAseqSingle(object):
    def __init__(self, **kwargs):
        """
        fq1
        genomt
        outdir
        fq2 (optional)

        align single file to reference genome
        """
        self.update(kwargs) # init self, fresh new
        # self.args = kwargs
        # self.args.update(self.__dict__) # for sub-func
        self.status = self.init_rnaseq() # update all variables: *.config, *.args


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def init_rnaseq(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        self.config = RNAseqConfig(**self.__dict__) # update, global
        # self.args.update(self.config.__dict__) # update, global
        self.update(self.config.__dict__, overwrite=False) # update local attributes
        # assert in_attr(self.config, ['fq1', 'genome', 'outdir'])
        # assert self.config.rnaseq_type == 'rnaseq_single'
        chka1 = isinstance(self.fq1, str)
        chka2 = isinstance(self.genome, str)
        chka3 = isinstance(self.outdir, str)
        if not all([chka1, chka2, chka3]):
            raise Exception('value error: \nfq1: {}\ngenome: {}\noutdir: {}'.format(
                self.fq1, self.genome, self.outdir))

        # save arguments
        chk1 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        args_logger(self.__dict__, self.config_txt)
        chk2 = True

        # status
        return all([chk1, chk2])


    #######################################
    ## main pipeline
    def prep_raw(self, copy=False):
        """
        Copy raw data to dest dir
        if not: create a symlink

        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.raw_fq_list

        # copy
        if copy is True:
            shutil.copy(self.fq1, raw_fq1)
            shutil.copy(self.fq2, raw_fq2)
        else:
            symlink(self.fq1, raw_fq1, absolute_path=True)
            symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=False):
        """
        Trim reads:
        hiseq.qc.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        args_trim = self.__dict__.copy()
        fq1, fq2 = self.raw_fq_list
        clean_fq1, clean_fq2 = self.clean_fq_list

        # update args
        # args = self.args.copy()         
        # args = self.config.args # all
        args_trim['fq1'] = args_trim['fq'] = fq1
        args_trim['fq2'] = fq2
        args_trim['outdir'] = self.cleandir

        if trimmed is True:
            # create symlink from rawdir
            # if raw is not gzipped, do it
            ## fq1
            if is_gz(fq1):
                symlink(fq1, clean_fq1)
            else:
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
            ## fq2
            if not fq2 is None:
                if is_gz(fq2):
                    symlink(fq2, clean_fq2)
                else:
                    gzip_cmd(fq2, clean_fq2, decompress=False, rm=False)

        else:
            if check_file(self.clean_fq_list):
                log.info('trim() skipped, file exists: {}'.format(
                    self.clean_fq_list))
            else:
                Trimmer(**args_trim).run()


    def align(self):
        """
        Alignment PE reads to reference genome, using STAR
        """
        args_align = self.__dict__.copy()
        fq1, fq2 = self.clean_fq_list

        # update arguments
        # args = self.config.args
        args_align['fq1'] = args_align['fq'] = fq1
        args_align['fq2'] = fq2
        args_align['outdir'] = self.aligndir
        args_align['aligner'] = 'STAR' # repeat,

        if check_file(self.bam_raw):
            log.info('align() skipped, file exists: {}'.format(
                self.bam_raw))
        else:
            Alignment(**args_align).run()


    def get_raw_bam(self):
        """
        Get the align bam file
        from bamdir/1., 2., ...
        !!! specific: 2.genome/*.bam
        """
        # bamdir = os.path.join(self.config.aligndir, '2.*')
        # bamdir = self.align_stat.rstrip('.align.txt')
        bamlist = listfile(self.aligndir, '*.bam', recursive=True)
        print('!AAAA1', self.bamdir, bamlist)
        bamlist = [b for b in bamlist if not b.endswith('.raw.bam')]
        print('!AAAA2', self.bamdir, bamlist)
        # [spike-in]? [rRNA, genome]
        return(bamlist[-1]) # last one


    def fc_count(self):
        """
        Run FeatureCounts for the bam file
        """
        args = self.__dict__.copy()

        # determine the strandness
        strand, strand_status = RNAseqLibrary(
            bam=self.get_raw_bam(),
            gtf=self.gtf).run(with_status=True)

        # save to status
        with open(self.strandness_status, 'w') as w:
            w.write(strand_status)

        # check strand
        if strand == 'read1':
            strand_fwd, strand_rev = (1, 2)
        elif strand == 'read2':
            strand_fwd, strand_rev = (2, 1)
        else:
            strand_fwd, strand_rev = (0, 0)

        # run sense
        args_fwd = {
            'gtf': self.config.gtf,
            'bam_list': self.get_raw_bam(),
            'outdir': self.countdir,
            'strandness': strand_fwd,
            'outname': 'count.sens.txt'}
        _, _, assign_fwd = FeatureCounts(**args_fwd).run()

        # run anti
        args_rev = {
            'gtf': self.config.gtf,
            'bam_list': self.get_raw_bam(),
            'outdir': self.countdir,
            'strandness': strand_rev,
            'outname': 'count.anti.txt'}
        _, _, assign_rev = FeatureCounts(**args_rev).run()


    def bam2bw(self):
        """
        Create bigWig files
        """
        pass


    def report(self):
        """
        Create alignment report for RNAseq
        1. trim
        2. align
        3. quant
        ...
        """
        pass


    def run(self):
        """
        Run pipeline for RNAseq,
        process
        """
        # init dir
        # args = self.config.args.copy()

        # copy_raw_fq = # args.get('copy_raw_fq', False)
        # trimmed = args.get('trimmed', False)
        copy_raw_fq = getattr(self, 'copy_raw_fq', False)
        trimmed = getattr(self, 'trimmed', False)        

        # 1. copy raw data
        self.prep_raw(copy_raw_fq)

        # 2. trim
        self.trim(trimmed)

        # 3. align
        self.align()

        # 4. quant
        self.fc_count()

        # 5. bam2bw
        self.bam2bw()

        # 6.report
        self.report()

        return self.outdir


class RNAseqMultiple(object):
    """
    Run RNAseq for each sample, one by one
    use RNAseqSingle
    """
    def __init__(self, **kwargs):
        """
        fq1
        fq2 (optional)
        genome
        outdir
        """
        # self.args = kwargs
        self.update(kwargs) # fresh new
        self.config = RNAseqConfig(**self.__dict__) # update, global
        self.update(self.config.__dict__, overwrite=False) # update unknown args
        self.status = self.init_rnaseq() # update all variables: *.config, *.args


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def init_rnaseq(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        # assert in_attr(self.config, ['fq1', 'genome', 'outdir'])
        chk1 = all([hasattr(self, i) for i in ['fq1', 'genome', 'outdir']])
        assert self.rnaseq_type == 'rnaseq_multiple'

        # check arguments
        chk2 = args_checker(self.__dict__, self.config_pickle)
        args = self.__dict__
        jfile = self.config_json
        print('!AAAA', jfile, args)
        # Json(args).writer(self.config_json)
        # Json(self.__dict__).writer(self.config_json)
        args_logger(self.__dict__, self.config_txt)

        # status
        return all([chk1, chk2])


    def run(self):
        """
        check each sample
        save to outdir/fq_name
        """
        args = self.__dict__.copy()

        ## run each sample
        ## Pool for parallel #!!!! features
        smp_dirs = []
        for i, fq1 in enumerate(args['fq1']):
            ## update required args
            args_i = args.copy()
            args_i['fq1'] = fq1
            args_i['fq2'] = args['fq2'][i]
            args_i['smp_name'] = args_i['smp_name'][i]
            args_i['outdir'] = os.path.join(args['outdir'], args_i['smp_name'])
            args_i['rnaseq_type'] = 'rnaseq_single'
            smp_dirs.append(args_i['outdir'])
            # make sure to run RNAseq single
            args_i['smp_path'] = args_i['dirs_ctl'] = args_i['dirs_exp'] = None
            args_i['pickle'] = None # if pickle input !!!
            RNAseqSingle(**args_i).run()

        return smp_dirs


class RNAseqDeseqSingle(object):
    """
    Run RNAseq for multiple samples

    required args:
        -dirs_ctl
        -dirs_exp
        -genme
        -outdir
        -feature

    config
    count/txt, rpkm, up, down, criteria
    de
    pdf/ma, pca, volcano
    report/html
    enrich/GO, GSEA, ...
    """
    def __init__(self, **kwargs):
        """
        dirs_ctl
        dirs_exp
        genome
        outdir
        feature
        group_ctl (None)
        group_exp (None)
        merge replicates
        """
        # self.args = kwargs
        # self.status = self.init_rnaseq() # update all variables: *.config, *.args
        self.update(kwargs) # fresh new
        self.status = self.init_rnaseq() # update all variables: *.config, *.args


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def init_rnaseq(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        self.config = RNAseqConfig(**self.args) # update, global
        self.update(self.config.__dict__, overwrite=False) # update local attributes
        # self.args.update(self.config.__dict__) # update, global
        assert self.rnaseq_type == 'deseq_single'

        # check arguments
        chk1 = args_checker(self.args, self.config_pickle)
        Json(self.__dict__).writer(self.config_json)
        args_logger(self.__dict__, self.config_txt)
        chk2 = True
        # self.args['overwrite'] = self.args.get('overwrite', False)
        # chk2 = self.args['overwrite'] is False

        # status
        return all([chk1, chk2])


    #######################################
    ## main pipeline
    def copy_count_files(self):
        """
        Copy count files to target dirs
        rename the count files by fqname
        """
        # get count_txt file
        for n, f in zip(self.prefix_ctl, self.count_ctl):
            # copy new file
            f_new = os.path.join(self.countdir, n + '.count_sens.txt')
            shutil.copy(f, f_new)

        for n, f in zip(self.prefix_exp, self.count_exp):
            # copy new file
            f_new = os.path.join(self.countdir, n + '.count_sens.txt')
            shutil.copy(f, f_new)


    ## create design.txt
    def get_design(self):
        """
        Create design.txt for this experiment
        colnames: group name gene count.txt
        """
        dlines = []
        for i, n in enumerate(self.smp_name):
            f_new = os.path.join(self.countdir, n + '.count_sens.txt')
            dlines.append('\t'.join(
                [self.group[i], n, self.feature, f_new]))
        # # control
        # for i, n in enumerate(self.prefix_ctl):
        #     f_new = os.path.join(self.countdir, n + '.count_sens.txt')
        #     dlines.append('\t'.join(
        #         [self.group[i], n, self.feature, f_new]))
        # # treatment
        # for i, n in enumerate(self.prefix_exp):
        #     f_new = os.path.join(self.countdir, n + '.count_sens.txt')
        #     dlines.append('\t'.join(
        #         [self.group[i], n, self.feature, f_new]))

        if os.path.exists(self.deseq_design):
            log.info('file exists - {}'.format(self.deseq_design))
        else:
            with open(self.deseq_design, 'wt') as w:
                w.write('\n'.join(dlines) + '\n')


    def deseq2(self):
        """
        ## quality control
        ## PCA/Cor/...
        """
        pass


    def run(self):
        """
        Run all
        """
        # 1. copy count
        self.copy_count_files()

        # 2. design
        self.get_design()

        # 3. run DESeq2
        self.deseq2()


class RNAseqDeseqMultiple(object):
    """
    args:
    smp_path:

    description:
    Create paires for input `smp_path`
    run RNAseqDeseqSingle() for each pair
    """
    def __init__(self, **kwargs):
        """
        smp_path: (output of RNAseqSingle)
        group: (optional, parse from basename(smp_path))
        """
        # self.args = kwargs
        self.update(kwargs) # fresh new
        self.status = self.init_rnaseq() # update all variables: *.config, *.args

        # self.smp_path = kwargs.get('smp_path', [])


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def init_rnaseq(self):
        if self.group is None:
            self.group = [fq_name_rmrep(i) for i in self.smp_path]
        self.group_pairs = design_combinations(self.group, n=2, return_index=True)

        ## check
        chk0 = len(self.smp_path) > 1
        chk1 = isinstance(self.smp_path, list)
        chk2 = len(self.group) > 1

        return all([chk0, chk1, chk2])


    def run(self):
        ## index for groups
        if self.status is True:
            for (ia, ib) in self.group_pairs:
                args_i = self.__dict__.copy()
                args_i['dirs_ctl'] = [self.smp_path[i] for i in ia]
                args_i['dirs_exp'] = [self.smp_path[i] for i in ib]
                args_i['smp_name'] = [self.smp_name[i] for i in ia + ib]
                args_i['group'] = [self.group[i] for i in ia + ib]
                args_i['smp_path'] = None # clear
                # update args
                config_i = RNAseqConfig(**args_i) # init
                args_i.update(config_i.__dict__) # update args
                RNAseqDeseqSingle(**args_i).run() # run single
        else:
            log.warning('DEseq multiple, skipped')


class RNAseqBuildDesign(object):
    """
    Create design.txt for experiment
    """
    def __init__(self, **kwargs):
        """
        required arguments:
        fq1, (fq2)
        genome
        outdir
        """
        self.update(kwargs) # fresh new
        self.config = RNAseqConfig(**self.__dict__)
        self.update(self.config.__dict__)
        self.status = self.init_rnaseq() # update all variables: *.config, *.args

        # self.args = kwargs
        # self.config = RNAseqConfig(**self.args)
        # self.args.update(self.config.__dict__)
        # self.status = self.init_rnaseq()


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def init_rnaseq(self):
        self.design = self.auto_design

        # check arguments
        chk1 = args_checker(self.__dict__, self.config_pickle)
        Json(self.__dict__).writer(self.config_json)
        args_logger(args, self.config_txt)
        args['overwrite'] = self.args.get('overwrite', False)
        chk2 = args['overwrite'] is False
        return all([chk1, chk2])


    def run(self):
        DesignBuilder(**self.__dict__).to_file(self.auto_design)
        ## log
        lines = '\n' + '#'*104 + '\n'
        lines += '# {:<100s} #\n'.format('Create design for RNAseq anslysis')
        lines += '# {:<100s} #\n'.format('1. design file:')
        lines += '# {:<100s} #\n'.format(self.design)
        lines += '# {:<100s} #\n'.format('2. arguments in pickle: ')
        lines += '# {:<100s} #\n'.format(self.config_pickle)
        lines += '# {:<100s} #\n'.format('')
        lines += '# {:<100s} #\n'.format('Run the following command to finish analysis:')
        lines += '# {:<100s} #\n'.format('')
        lines += '$ hiseq rnaseq -d {} \n'.format(self.design)
        lines += 'or\n'
        lines += '$ hiseq rnaseq --pickle {}\n'.format(self.config_pickle)
        log.info(lines)


class RNAseq(object):
    def __init__(self, **kwargs):
        """
        options:
        global options: feature(gene/te/piRNAcluster/...), gtf,
        1. design: design.txt (single/multiple/deseq)
        2. single/multiple: --fq1, --fq2, --genome, --outdir
        3. deseq: --dirs-ctl, --dirs-exp, --outdir
        ...
        """
        # self.update(kwargs) # fresh new
        # # self.args = kwargs
        # self.config = RNAseqConfig(**self.__dict__) # init
        # self.args.update(self.config.__dict__) # update global
        # self.args['pickle'] = None # terminate `pickle`: top-level, after 1st round RNAseqConfig(), !!!
        # self.group = self.args.get('group', [])

        self.update(kwargs) # fresh new
        self.config = RNAseqConfig(**self.__dict__)
        self.update(self.config.__dict__, overwrite=False)
        self.pickle = None # terminate `pickle` option
        # self.status = self.init_rnaseq() # update all variables: *.config, *.args


    def update(self, d, overwrite=True):
        """
        Update attributes from dict
        overwrite exists attr
        """
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or overwrite:
                    setattr(self, k, v)


    def run(self):
        """
        Run all RNAseq analysis
        """
        args = self.__dict__.copy()

        if self.config.rnaseq_type == 'build_design':
            RNAseqBuildDesign(**args).run()
        elif self.config.rnaseq_type == 'deseq_single':
            RNAseqDeseqSingle(**args).run()
        elif self.config.rnaseq_type == 'deseq_multiple':
            RNAseqDeseqMultiple(**args).run()
        elif self.config.rnaseq_type == 'rnaseq_single':
            RNAseqSingle(**args).run()
        elif self.config.rnaseq_type == 'rnaseq_multiple':
            # env: rnaseq_multiple
            smp_path = RNAseqMultiple(**args).run()
            # env: deseq_multiple
            args_i = args.copy()
            args_i['smp_path'] = smp_path
            self.config = RNAseqConfig(**args_i)
            args_i.update(self.config.__dict__)
            args_i['smp_path'] = smp_path # covered by RNAseqConfig
            RNAseqDeseqMultiple(**args_i).run()
        else:
            pass


class RNAseqLibrary(object):
    """
    Determine the RNAseq library type:
    strandness: 1 ++, 1 --, / 2 +-, 2 -+ :
    dUTP, NSR: 1 +-, 1 -+, / 2 ++, 2 -- :
    ###
    or:
    infer_experiment.py from RSeQC package.
    """
    def __init__(self, bam, gtf, size=200000, **kwargs):
        self.bam = bam
        self.gtf = gtf
        self.size = size
        self.args = kwargs
        self.config()


    def config(self):
        """
        Prepare Bam file (subsset)
        """
        if isinstance(self.bam, str):
            pass
        elif isinstance(self.bam, list):
            self.bam = self.bam[0] # the first one
        else:
            raise Exception('Expect str and list, get {}'.format(type(bam)))

        self.outdir = self.args.get('outdir', None)
        if self.outdir is None:
            self.outdir = self._tmp()

        self.bam_sub = self.subset_bam(self.outdir, self.size) # subset

        ## subset bam list: 20000 reads
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        ## clean
        self.cleanup = self.args.get('cleanup', True)

        ## input a BED file
        if self.gtf.lower().endswith('.bed'):
            self.bed = self.gtf
            self.gtf = os.path.join(self.outdir,
                os.path.basename(os.path.splitext(self.gtf)[0]) + '.gtf')
            # convert
            self.bed2gtf(self.bed, self.gtf)


    def _tmp(self):
        """
        Create a tmp filename
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=True)
        return tmp.name


    def bed2gtf(self, bed_in, gtf_out):
        """
        Convert BED to GTF format
        """
        with open(bed_in) as r, open(gtf_out) as w:
            for line in r:
                fields = line.rstrip().split('\t')
                start = int(fields[1]) + 1 # 0-index to 1-index
                end = fields[2]
                des = 'gene_id "{}"; gene_name "{}";'.format(fields[3], fields[3])
                gtf_line = '\t'.join([
                    fields[0],
                    'BED_file',
                    'exon',
                    str(start),
                    end,
                    '.',
                    fields[5],
                    '.',
                    des])
                w.write(gtf_line + '\n')


    def subset_bam(self, dest, size=20000):
        """
        Extract N reads from bam list
        """
        if not os.path.exists(dest):
            os.mkdir(dest)

        ## create bam index
        Bam(self.bam).index()
        bname = os.path.basename(self.bam)
        dest_bam = os.path.join(dest, bname)
        samfile = pysam.AlignmentFile(self.bam, 'rb')
        subfile = pysam.AlignmentFile(dest_bam, 'wb', template=samfile)

        # count
        n = 0
        for read in samfile.fetch():
            n += 1
            if n > size:
                break
            subfile.write(read)

        return dest_bam


    def stat(self, strandness=1):
        """
        strandness; 1 or 2; for featureCounts

        using -s 1, -s 2
        and return the mapping reads
        """
        ## -s 1
        outname = 'count.s_' + str(strandness) + '.txt'
        args = {
            'gtf': self.gtf,
            'bam_list': self.bam_sub,
            'outdir': self.outdir,
            'strandness': strandness,
            'outname': outname}

        total, assign, assign_df = FeatureCounts(**args).run()

        return assign


    def run(self, with_status=False):
        """
        Check -s 1, -s 2:
        """
        a = self.stat(strandness=1).to_list()[0] # s=1
        b = self.stat(strandness=2).to_list()[0] # s=2

        ## loginfo
        strand_status = 'Library strandness by featureCounts: -s\n'
        strand_status += 'sample: {}\n'.format(self.bam)
        strand_status += 'gtf: {}\n'.format(self.gtf)
        strand_status += 'assigned_fwd: {:>6.2f}%  -s=1  stranded\n'.format(a)
        strand_status += 'assigned_rev: {:>6.2f}%  -s=2  reversely_stranded\n'.format(b)
        log.info(strand_status)

        ## correct BED
        if a + b < 20:
            log.warning('check out gtf file: {}'.format(self.gtf))
            tag = 'undetermined'
        elif a == b:
            tag = 'read0'
        elif a > b:
            tag = 'read1'
        else:
            tag = 'read2'

        if self.cleanup is True:
            # check ?
            # os.remove(self.outdir)
            log.warning('Temporary directory in: {}'.format(self.outdir))
            if self.outdir.startswith('/tmp/'):
                # remove files within /tmp folder # linux
                shutil.rmtree(self.outdir) # dangerous !!!

        if with_status:
            return (tag, strand_status)
        else:
            return tag


class FeatureCounts(object):
    """
    Run featureCounts for GTF + BAM(s)
    mission()
    count.txt
    determin: library type: +-/-+; -+/+-;
    map reads, scale
    FPKM/RPKM;
    map pct:
    """
    def __init__(self, **kwargs):
        """
        arguments:
        gtf:
        bam_list: (index)
        outdir:
        strandness: 0=no, 1=sens, 2=anti, 3=both
        smpname: count.txt default
        threads: 4
        overwrite: False
        """
        # required args
        self.args = kwargs
        self.check_status = self.init_fc() # global config


    def init_fc(self):
        """
        Initiate the config for RNAseq analysis
        """
        args = self.args.copy()

        # required
        self.gtf = args.get('gtf', None)
        self.bam_list = args.get('bam_list', None)
        self.outdir = args.get('outdir', str(pathlib.Path.cwd()))
        assert is_path(self.outdir)

        # convert to list
        if isinstance(self.bam_list, str):
            self.bam_list = [self.bam_list]

        # optional
        self.strandness = args.get('strandness', 0)
        self.threads = args.get('threads', 4)
        self.overwrite = args.get('overwrite', False)
        self.outname = args.get('outname', 'count.txt')

        # check required args
        chk1 = isinstance(self.gtf, str) and os.path.exists(self.gtf)
        chk2 = all([os.path.exists(i) for i in self.bam_list])
        if not chk1 or not chk2:
            log.error('fc() failed, check gtf={}, bam_list={}'.format(self.gtf, self.bam_list))

        ## absolute path
        self.gtf = os.path.abspath(self.gtf)
        self.bam_list = [os.path.abspath(i) for i in self.bam_list]
        self.outdir = os.path.abspath(self.outdir)

        ## index bam files
        [Bam(i).index() for i in self.bam_list]

        # determine the output files
        self.count_txt = os.path.join(self.outdir, self.outname)
        self.summary = self.count_txt + '.summary'
        self.log_file = os.path.join(self.outdir, 'featureCounts.log')
        ## summary

        ## PE reads
        self.pe_flag = self.is_PE_bam()

        ## tool
        self.fc_exe = shutil.which('featureCounts')

        ## update args
        args = self.__dict__
        args.pop('args_input', None)
        args.pop('cmd_input', None)
        args.pop('dict', None)
        args.pop('args', None)
        ## save config
        args_txt = os.path.join(self.outdir, 'arguments.txt')
        args_json = os.path.join(self.outdir, 'arguments.json')
        args_pickle = os.path.join(self.outdir, 'arguments.pickle')
        # save to file
        Json(args).writer(args_json)
        args_logger(args, args_txt)
        chk3 = args_checker(args, args_pickle)
        chk4 = args.get('overwrite')

        # status
        return all([chk3, chk4])


    def get_cmd(self):
        """
        prepare args for featureCounts
        """
        cmd = '{} -M -O --fraction -g gene_id -t exon '.format(
            self.fc_exe)

        # for PE
        if self.pe_flag:
            cmd += ' -p -C -B '

        # args
        cmd += '-s {} -T {} -a {} -o {} {} 2> {}'.format(
            self.strandness,
            self.threads,
            self.gtf,
            self.count_txt,
            ' '.join(self.bam_list),
            self.log_file)

        return cmd


    def is_PE_bam(self):
        """
        Check whether the input bam files are Paired or Single file
        Bam().isPaired()
        """
        return all([Bam(i).isPaired() for i in self.bam_list])


    def wrap_log(self):
        """
        save output file to log,
        """
        df = pd.read_csv(self.summary, '\t', index_col=0)
        df.columns = [os.path.basename(i) for i in df.columns.tolist()]
        ## total
        total = df.sum(axis=0, skipna=True)
        ## assign
        assign = df.loc['Assigned', ] / total * 100
        ## pct
        assign_df = assign.to_frame('assigned%')
        ## minimum
        min_pct = assign.min()
        if assign.min() < 50:
            log.warning('caution: -s {}, {:.2f}% reads were assigned, {}'.format(
                self.strandness,
                assign.min(),
                self.summary))
        return (total, assign, assign_df)


    def _tmp(self):
        """
        Create a tmp filename
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=True)
        return os.path.basename(tmp.name)


    def run(self):
        """
        run featureCounts
        """
        cmd = self.get_cmd()
        if os.path.exists(self.count_txt) and self.overwrite is False:
            log.info('FeatureCounts() skipped, file exists: {}'.format(
            self.count_txt))
        else:
            run_shell_cmd(cmd)

        return self.wrap_log()

