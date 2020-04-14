#!/usr/bin/env python
# -*- coding: utf-8 -*-



"""
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
        self.args = kwargs
        self.rnaseq_type = self.mission_type()

        # for RNAseq checker
        create_dirs = kwargs.get('create_dirs', True)

        if self.rnaseq_type == 'single':
            self.init_rnaseq_single(create_dirs)
        elif self.rnaseq_type == 'multiple':
            self.init_rnaseq_multiple(create_dirs)
        elif self.rnaseq_type == 'deseq':
            self.init_rnaseq_deseq(create_dirs)
        else:
            log.error('unknown')
            pass


    def mission_type(self):
        """
        Determine the purpose the RNAseq analysis
        1. single 
        2. pair
        ...
        """
        ## check design
        design = self.args.get('design', None) # 1 single
        if not design is None:
            args_in = self.design_reader(design) # tab-sep file
            self.args.update(args_in)

        args = self.args.copy()
        fq1 = args.get('fq1', None)

        if isinstance(fq1, str):
            flag = 'single'
        elif isinstance(fq1, list):
            flag = 'multiple'
        else:
            flag = 'deseq'

        # if not design is None:
        #     args_in = self.design_reader(args['design']) # tab-sep file
        #     self.args.update(args_in) # update, global
        #     flag = 'single'
        # elif all([not i is None for i in [fq1, fq2, genome, outdir]]):
        #     flag = 'single'
        # elif isinstance(rep_list, list):
        #     flag = 'multiple'
        # elif isinstance(smp_list, list):
        #     flag = 'deseq'
        # else:
        #     raise Exception("""unknown RNAseq() arguments:
        #         single: config, design, fq1, fq2, genome, outdir; 
        #         merge: rep_list; 
        #         deseq: smp_list, the directory of merged list;
        #         """)

        return flag


    def design_reader(self, x):
        """Read the fastq and genome information
        group, name, genome, outdir, fq1, fq2
        """
        df = pd.read_csv(x, '\t', header=None, comment='#', skip_blank_lines=True,
            names=['group', 'name', 'genome', 'outdir', 'fq1', 'fq2'])
        
        # unique names
        names = df['name'].to_list()
        chk1 = len(names) == len(set(names))

        # unique genome: single, str
        genomes = df['genome'].to_list()
        chk2 = len(set(genomes)) == 1

        # outdir: single, str
        outdirs = df['outdir'].to_list()
        chk3 = len(set(outdirs)) == 1

        # fq1: str
        fq1_list = df['fq1'].to_list()
        chk4 = [os.path.exists(i) for i in fq1_list]

        # fq2: NA or str
        fq2_list = df['fq2'].to_list()

        # construct args
        args = {
            'fq1': fq1_list,
            'fq2': fq2_list,
            'genome': genomes[0],
            'outdir': outdirs[0],
            'smp_name': names,
        }

        # for check
        if not all(chk1, chk2, chk3, chk4):
            raise Exception('{} - file format not correct'.format(x))

        return args


    def init_rnaseq_single(self, create_dirs=True):
        """
        Initiate the config, directories, files for rnaseq single
        update self.
        fq1, fq2, genome, outdir: args, 
        or:
        design
        """
        args = self.args.copy() # global
        # self.rnaseq_type = 'single'

        self.fq1 = os.path.abspath(args['fq1'])
        self.fq2 = os.path.abspath(args['fq2']) if args['fq2'] else args['fq2']
        self.genome = args['genome']
        self.outdir = os.path.abspath(args['outdir'])
        self.groups = args.get('groups', 'gene')

        ## sample name
        self.fqname = args['smp_name'] if args.get('smp_name', None) else fq_name(args['fq1'])

        ## get GTF file
        self.gtf = args.get('gtf', None)
        if self.gtf is None:
            self.gtf = Genome(genome=self.genome).gene_gtf('ensembl') # ucsc version

        ## outdir
        self.configdir = os.path.join(self.outdir, self.groups, 'config')
        self.rawdir = os.path.join(self.outdir, self.groups, 'raw_data')
        self.cleandir = os.path.join(self.outdir, self.groups, 'clean_data')
        self.aligndir = os.path.join(self.outdir, self.groups, 'align')
        self.bamdir = os.path.join(self.outdir, self.groups, 'bam_files')
        self.bwdir = os.path.join(self.outdir, self.groups, 'bw_files')
        self.countdir = os.path.join(self.outdir, self.groups, 'count')
        self.reportdir = os.path.join(self.outdir, self.groups, 'report')
        self.out_prefix = os.path.join(self.outdir, self.groups, self.fqname)

        ## raw
        self.raw_fq_list = [
            os.path.join(self.rawdir, os.path.basename(self.fq1)),
            os.path.join(self.rawdir, os.path.basename(self.fq2))]

        ## clean = fq.gz
        self.clean_fq_list = [
            os.path.join(self.cleandir, file_prefix(self.fq1)[0] + '.fq.gz'),
            os.path.join(self.cleandir, file_prefix(self.fq2)[0] + '.fq.gz')]
        
        ## 
        self.trim_stat = os.path.join(self.cleandir, self.fqname + '.qc.stat')
        self.bam_raw = os.path.join(self.aligndir, self.fqname, '2.*', self.fqname + '.bam')
        self.align_stat = os.path.join(self.aligndir, self.fqname + '.align.txt')
        self.bw_fwd = os.path.join(self.bwdir, self.fqname + '.fwd.bigWig')
        self.bw_rev = os.path.join(self.bwdir, self.fqname + '.rev.bigWig')
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

        ## update all parameters, in self


    def init_rnaseq_multiple(self, create_dirs=True):
        """
        Initiate the config, for multiple RNAseq samples
        """
        args = self.args.copy() # global

        # fastq1 files
        if isinstance(args['fq1'], list):
            self.fq1 = [os.path.abspath(i) for i in args['fq1']]
        elif isinstance(args['fq1'], str):
            self.fq1 = [os.path.abspath(args['fq1'])]
        else:
            raise Exception('unknown fq1: {}'.args['fq1'])

        # fastq2 files
        if args['fq2'] is None:
            self.fq2 = [None] * len(fq1)
        elif isinstance(args['fq2'], list):
            self.fq2 = [os.path.abspath(i) for i in args['fq2']]
        elif isinstance(args['fq2'], str):
            self.fq2 = [os.path.abspath(args['fq2'])]
        else:
            raise Exception('unknown fq2: {}'.args['fq2'])

        ## required 
        self.genome = args['genome']
        self.outdir = os.path.abspath(args['outdir'])
        self.groups = args.get('groups', 'gene')

        ## sample name
        self.fqname_list = fq_name(self.fq1)
        if args.get('smp_name', None):
            if len(self.fq1) == len(args['smp_name']):
                self.fqname_list = args['smp_name']

        ## get GTF file
        self.gtf = args.get('gtf', None)
        if self.gtf is None:
            self.gtf = Genome(genome=self.genome).gene_gtf('ensembl') # ucsc version


    def init_rnaseq_deseq(self, create_dirs=True):
        pass


    def is_rnaseq_single(self, x=None):
        """
        Check if the directory is of RNAseq single sample
        outdir/group/config/arguments.pickle
        """
        args = self.args.copy()
        outdir = args.get('outdir', str(pathlib.Path.cwd()))
        groups = args.get('groups', 'gene')
        x_pickle = x if x else os.path.join(outdir, groups, 'config', 'arguments.pickle')

        tag = None
        if os.path.exits(x_pickls):
            args_config = pickle_to_dict(x_pickle)
            tag = args_config.get('rnaseq_type', None) == 'single'

        return tag


    def is_rnaseq_multi(self, x=None):
        """
        Check if the directory is of RNAseq multiple sample
        outdir/group/config/arguments.pickle
        """
        args = self.args.copy()
        outdir = args.get('outdir', str(pathlib.Path.cwd()))
        groups = args.get('groups', 'gene')
        x_pickle = x if x else outdir + '.config.pickle'

        tag = None
        if os.path.exits(x_pickls):
            args_config = pickle_to_dict(x_pickle)
            tag = args_config.get('rnaseq_type', None) == 'multiple'

        return tag


    def is_rnaseq_deseq(self, x):
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
        self.args = kwargs
        self.status = self.init_rnaseq() # update all variables: *.config, *.args
        

    def init_rnaseq(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        self.config = RNAseqConfig(**self.args) # update, global
        self.args.update(self.config.__dict__) # update, global
        assert in_attr(self.config, ['fq1', 'genome', 'outdir'])
        assert self.config.rnaseq_type == 'single'

        # save config to files
        args_pickle = os.path.join(self.config.configdir, 'arguments.pickle')
        args_json = os.path.join(self.config.configdir, 'arguments.json')
        args_txt = os.path.join(self.config.configdir, 'arguments.txt')

        # check arguments
        chk1 = args_checker(self.args, args_pickle)
        Json(self.args).writer(args_json)
        args_logger(self.args, args_txt)
        self.args['overwrite'] = self.args.get('overwrite', False)
        chk2 = self.args['overwrite'] is False

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
        raw_fq1, raw_fq2 = self.config.raw_fq_list

        # copy
        if copy is True:
            shutil.copy(self.config.fq1, raw_fq1)
            shutil.copy(self.config.fq2, raw_fq2)
        else:
            symlink(self.config.fq1, raw_fq1, absolute_path=True)
            symlink(self.config.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=False):
        """
        Trim reads:
        hiseq.qc.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        # args = self.args.copy()
        fq1, fq2 = self.config.raw_fq_list
        clean_fq1, clean_fq2 = self.config.clean_fq_list

        # update args
        args = self.config.args # all
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.cleandir
        
        if trimmed is True:
            # create symlink from rawdir
            # if raw is not gzipped, do it
            if is_gz(fq1) and is_gz(fq2):
                symlink(fq1, clean_fq1)
                symlink(fq2, clean_fq2)
            else:
                # gzip files
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
                gzip_cmd(fq2, clean_fq2, decompress=False, rm=False)
        else:
            if check_file(self.config.clean_fq_list):
                log.info('trim() skipped, file exists: {}'.format(
                    self.config.clean_fq_list))
            else:
                Trimmer(**args).run()


    def align(self):
        """
        Alignment PE reads to reference genome, using STAR
        """
        fq1, fq2 = self.config.clean_fq_list

        # update arguments
        args = self.config.args
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.aligndir
        args['aligner'] = 'STAR' # repeat,

        if check_file(self.config.bam_raw):
            log.info('align() skipped, file exists: {}'.format(
                self.config.bam_raw))
        else:
            Alignment(**args).run()


    def get_raw_bam(self):
        """
        Get the align bam file
        from bamdir/1., 2., ...
        !!! specific: 2.genome/*.bam
        """
        # bamdir = os.path.join(self.config.aligndir, '2.*')
        bamdir = self.config.align_stat.rstrip('.align.txt')
        bamlist = listfile(bamdir, '*.bam', recursive=True)
        bamlist = [b for b in bamlist if not b.endswith('.raw.bam')]
        # [spike-in]? [rRNA, genome]
        return(bamlist[-1]) # genome


    def fc_count(self):
        """
        Run FeatureCounts for the bam file
        """
        args = self.args.copy()

        # determine the strandness
        strand, strand_status = RNAseqLibrary(
            bam=self.get_raw_bam(), 
            gtf=self.config.gtf).run(with_status=True)

        # save to status
        with open(self.config.strandness_status, 'w') as w:
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
            'outdir': self.config.countdir,
            'strandness': strand_fwd,
            'outname': 'count.sens.txt'}
        _, _, assign_fwd = FeatureCounts(**args_fwd).run()

        # run anti
        args_rev = {
            'gtf': self.config.gtf, 
            'bam_list': self.get_raw_bam(),
            'outdir': self.config.countdir,
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
        args = self.config.args.copy()

        copy_raw_fq = args.get('copy_raw_fq', False)
        trimmed = args.get('trimmed', False)

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
        self.args = kwargs
        self.status = self.init_rnaseq() # update all variables: *.config, *.args


    def init_rnaseq(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        self.config = RNAseqConfig(**self.args) # update, global
        self.args.update(self.config.__dict__) # update, global
        assert in_attr(self.config, ['fq1', 'genome', 'outdir'])
        assert self.config.rnaseq_type == 'multiple'

        # save config to files
        args_pickle = os.path.join(self.config.outdir, 'arguments.pickle')
        args_json = os.path.join(self.config.outdir, 'arguments.json')
        args_txt = os.path.join(self.config.outdir, 'arguments.txt')

        # check arguments
        chk1 = args_checker(self.args, args_pickle)
        Json(self.args).writer(args_json)
        args_logger(self.args, args_txt)
        self.args['overwrite'] = self.args.get('overwrite', False)
        chk2 = self.args['overwrite'] is False

        # status
        return all([chk1, chk2])

     
    def run(self):
        """
        check each sample
        save to outdir/fq_name
        """
        args = self.args.copy()

        ## run each sample
        smp_dirs = []
        for i, fq1 in enumerate(args['fq1']):
            ## update required args
            args_i = args.copy()
            args_i['fq1'] = fq1
            args_i['fq2'] = args['fq2'][i]
            args_i['fqname'] = args['fqname_list'][i]
            args_i['smp_name'] = args_i['fqname']
            args_i['outdir'] = os.path.join(args['outdir'], args_i['fqname'])
            args_i['rnaseq_type'] = 'single'
            smp_dirs.append(args_i['outdir'])
            RNAseqSingle(**args_i).run()

        return smp_dirs


class RNAseqDeseq(object):
    """
    Run RNAseq for multiple samples
    config
    count/txt, rpkm, up, down, criteria
    de
    pdf/ma, pca, volcano
    report/html
    enrich/GO, GSEA, ...
    """
    def __init__(self, **kwargs):
        pass


class RNAseq(object):
    def __init__(self, **kwargs):
        pass


