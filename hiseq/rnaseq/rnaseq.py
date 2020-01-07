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
            # self.is_atac_single(self.outdir)
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

