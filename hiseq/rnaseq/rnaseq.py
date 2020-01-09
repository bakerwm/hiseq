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
from hiseq.utils.helper import *
from hiseq.qc.trimmer import Trimmer
from hiseq.align.alignment import Alignment


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
        self.rnaseq_type = self.mission()

        # for RNAseq checker
        create_dirs = kwargs.get('create_dirs', True)

        if self.rnaseq_type == 'single':
            self.init_rnaseq_single(create_dirs)
        elif self.rnaseq_type == 'merge':
            self.init_rnaseq_merge(create_dirs)
        elif self.ranseq_type = 'deseq':
            self.init_rnaseq_deseq(create_dirs)
        else:
            print('unknown')
            pass


    def mission(self):
        """
        Determine the purpose the RNAseq analysis
        1. single 
        2. pair
        ...
        """
        args = self.args.copy()

        config = args.get('config', None) # N single
        design = args.get('design', None) # 1 single
        fq1 = args.get('fq1', None)
        fq2 = args.get('fq2', None)
        genome = args.get('genome', None)
        outdir = args.get('outdir', None)
        rep_list = args.get('rep_list', None)
        input_dir = args.get('smp_list', None)

        ## check working type
        if not design is None:
            args_in = Json(args['design']).dict
            self.args.update(args_in) # update, global
            flag = 'single'
        elif all([not i is None for i in [fq1, fq2, genome, outdir]]):
            flag = 'single'
        elif isinstance(rep_list, list):
            flag = 'merge'
        elif isinstance(smp_list, list):
            flag = 'deseq'
        else:
            raise Exception("""unknown RNAseq() arguments:
                single: config, design, fq1, fq2, genome, outdir; 
                merge: rep_list; 
                deseq: smp_list, the directory of merged list;
                """)


    def init_rnaseq_single(self, create_dirs=True):
        """
        initiate the config, directories, files for rnaseq single
        update self.
        fq1, fq2, genome, outdir: args, 
        or:
        design
        """
        args = self.args.copy() # global

        self.fq1 = args['fq1']
        self.fq2 = args['fq2']
        self.genome = args['genome']
        self.outdir = args['outdir']

        ## absolute path
        self.fq1 = os.path.abspath(self.fq1)
        self.fq2 = os.path.abspath(self.fq2)
        self.outdir = os.path.abspath(self.outdir)

        ## sample name
        args['smp_name'] = args.get('smp_name', None)
        fqname = file_prefix(args['fq1'])[0]
        fqname = re.sub('[._][rR]?1$', '', fqname)
        fqname = re.sub('_1$', '', fqname)
        if args['smp_name']:
            fqname = args['smp_name']
        self.fqname = fqname

        ## outdir
        self.configdir = os.path.join(self.outdir, 'config')
        self.rawdir = os.path.join(self.outdir, 'raw_data')
        self.cleandir = os.path.join(self.outdir, 'clean_data')
        self.aligndir = os.path.join(self.outdir, 'align')
        self.bamdir = os.path.join(self.outdir, 'bam_files')
        self.bwdir = os.path.join(self.outdir, 'bw_files')
        self.peakdir = os.path.join(self.outdir, 'peak')
        self.motifdir = os.path.join(self.outdir, 'motif')
        self.qcdir = os.path.join(self.outdir, 'qc')
        self.reportdir = os.path.join(self.outdir, 'report')
        self.out_prefix = os.path.join(self.outdir, fqname)

        # self.raw_fq_list = [os.path.join(self.rawdir, i) for i in fqnames]
        # self.clean_fq_list = [os.path.join(self.cleandir, i) for i in fqnames]
        ## fastq files
        ## consider input:
        ## input: fastq, fq, fastq.gz, fq.gz
        ## raw = input
        ## clean = *.fq.gz # gzip if required.
        # fqnames = list(map(os.path.basename, [self.fq1, self.fq2]))
        self.raw_fq_list = [
            os.path.join(self.rawdir, os.path.basename(self.fq1)),
            os.path.join(self.rawdir, os.path.basename(self.fq2))]
        ## clean = fq.gz
        self.clean_fq_list = [
            os.path.join(self.cleandir, file_prefix(self.fq1)[0] + '.fq.gz'),
            os.path.join(self.cleandir, file_prefix(self.fq2)[0] + '.fq.gz')]
        ## 
        self.trim_stat = os.path.join(self.cleandir, fqname + '.qc.stat')
        self.bam_raw = os.path.join(self.aligndir, fqname, '2.*', fqname + '.bam')
        self.align_stat = os.path.join(self.aligndir, fqname + '.align.txt')
        self.bam_rmdup = os.path.join(self.bamdir, fqname + '.rmdup.bam')
        self.bam_proper_pair = os.path.join(self.bamdir, fqname + '.proper_pair.bam')
        self.peak = os.path.join(self.peakdir, fqname + '_peaks.narrowPeak')
        self.bw = os.path.join(self.bwdir, fqname + '.bigWig')
        ## qc files
        self.lendist_txt = os.path.join(self.qcdir, 'length_distribution.txt')
        self.lendist_pdf = os.path.join(self.qcdir, 'length_distribution.pdf')
        self.frip_txt = os.path.join(self.qcdir, 'FRiP.txt')

        ## update args
        args['overwrite'] = args.get('overwrite', False)
        args['threads'] = args.get('threads', 8)
        args['fqname'] = self.fqname

        ## trimming, cutadapt
        adapter3 = hiseq.utils.args.Adapter('nextera').adapters[0] #
        args['len_min']  = args.get('len_min', 20)
        args['adapter3'] = args.get('adapter3', adapter3) # nextera
        args['AD3'] = args.get('AD3', adapter3) # nextera

        ## alignment
        args['aligner'] = args.get('aligner', 'bowtie2') # bowtie alignment
        args['n_map'] = args.get('n_map', 2)
        args['align_to_chrM'] = True
        args['extra_para'] = '-X 2000'

        ## update args
        args['rawdir'] = self.rawdir
        args['cleandir'] = self.cleandir
        args['aligndir'] = self.aligndir
        args['bamdir'] = self.bamdir
        args['bwdir'] = self.bwdir
        args['peakdir'] = self.peakdir
        args['motifdir'] = self.motifdir
        args['reportdir'] = self.reportdir
        args['out_prefix'] = self.out_prefix 
        args['raw_fq_list'] = self.raw_fq_list
        args['clean_fq_list'] = self.clean_fq_list
        args['trim_stat'] = self.trim_stat
        args['bam_raw'] = self.bam_raw
        args['align_stat'] = self.align_stat
        args['bam_rmdup'] = self.bam_rmdup
        args['bam_proper_pair'] = self.bam_proper_pair
        args['peak'] = self.peak
        args['bw'] = self.bw

        # update, global
        self.args = args

        ## create directories
        if create_dirs is True:
            check_path([
                self.configdir,
                self.rawdir, 
                self.cleandir, 
                self.aligndir, 
                self.bamdir, 
                self.bwdir, 
                self.peakdir, 
                self.motifdir, 
                self.qcdir, 
                self.reportdir])





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
        self.status = self.init_atac() # update all variables: *.config, *.args
        

    def init_rnaseq(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        self.config = RNAseqConfig(**self.args) # update, global
        self.args.update(self.config.args) # update, global
        assert in_attr(self.config, ['fq1', 'genome', 'outdir'])
        assert self.config.atac_type == 'single'

        # required
        self.fq1 = os.path.abspath(self.args.get('fq1', None))
        # self.fq2 = os.path.abspath(self.args.get('fq2', None))
        self.genome = self.args.get('genome', None)
        self.outdir = os.path.abspath(self.args.get('outdir', None))
        self.args['fq1'] = self.fq1
        self.args['fq2'] = self.fq2
        self.args['outdir'] = self.outdir

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
        # args = self.args.copy()
        fq1, fq2 = self.config.raw_fq_list
        clean_fq1, clean_fq2 = self.config.clean_fq_list

        # update args
        args = self.config.args # all
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.cleandir
        
        if trimmed is True:
            print('!xxxx')
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
            print('!yyyy')
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
        args['extra_para'] = '-X 2000' # repeat, 
        args['aligner'] = 'STAR' # repeat,

        if check_file(self.config.bam_rmdup):
            log.info('align() skipped, file exists: {}'.format(
                self.config.bam_rmdup))
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
        bamlist = listfiles2('*.bam', bamdir, recursive=True)
        # [chrM, genome]
        return(bamlist[1]) # genome



