#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
RNAseq pipeline: level-3 (run _rn)

fastq files from: fq1/fq2

analysis-module:
"""

import os
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.rnaseq.rnaseq_r1 import RnaseqR1
from hiseq.rnaseq.rnaseq_rp import RnaseqRp
from hiseq.rnaseq.utils import rnaseq_trim, rnaseq_align_spikein, \
    rnaseq_align_rRNA, rnaseq_align_genome, rnaseq_quant, rnaseq_merge_bam, \
    rnaseq_bam2bw, qc_trim_summary, qc_align_summary, qc_bam_cor, \
    qc_genebody_enrich

from hiseq.utils.file import check_path, check_fx_paired, symlink_file, \
    file_exists, file_abspath, file_prefix, fx_name, Genome
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu, \
    read_hiseq, list_hiseq_file, run_shell_cmd
from hiseq.utils.bam import Bam
from hiseq.align.align_index import AlignIndex



class RnaseqRn(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RnaseqRnConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)


    def run_pipe_rn(self): # for rep_list > 1
        rnaseq_merge_bam(self.project_dir, 'rn')
        rnaseq_bam2bw(self.project_dir, 'rn')
        rnaseq_quant(self.project_dir, 'rn')
        qc_genebody_enrich(self.project_dir, 'rn')
        qc_bam_cor(self.project_dir, 'rn')
        

    def run_single_fx(self, i):
        # required args
        args_required = [
            'aligner', 'fq1', 'fq2', 'is_mut', 'outdir', 'gene_bed', 'gene_gtf', 
            'extra_index', 'genome', 'genome_index', 'spikein', 'spikein_index',
            'rRNA_index', 'to_rRNA', 'genome_size', 'genome_size_file',
            'threads', 'parallel_jobs', 'overwrite', 'binsize', 
            'trimmed', 'cut_to_length', 'recursive'
        ]
        args = dict((k, getattr(self, k)) for k in args_required \
            if k in self.__dict__)
        # update fq1, fq2, rep_list, ...
        fq1 = self.fq1[i]
        fq2 = None if self.fq2 is None else self.fq2[i]
        args.update({
            'fq1': fq1,
            'fq2': fq2,
        })
        RnaseqR1(**args).run()
        

    def run_multi_fx(self):
        for i in range(len(self.fq1)):
            self.run_single_fx(i)

        
    def run(self):
        # 1. save config
        Config().dump(self.__dict__, self.config_yaml)
        # 2. run CnrR1
        self.run_multi_fx()
        # 3. run CnrRn, merge
        self.run_pipe_rn()
        # 4. generate reprot
        RnaseqRp(**self.__dict__).run()

        
class RnaseqRnConfig(object):
    """
    prepare for RnaseqRn() command
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for RNAseq analysis
        """
        args_init = {
            'aligner': 'STAR',
            'is_mut': False,
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'smp_name': None,
            'norm_project': None,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': None,
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'to_rRNA': False,
            'rRNA_index': None,
            'extra_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_rn' #
        # aligner
        if self.aligner is None:
            self.aligner = 'STAR'
        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        if not isinstance(self.smp_name, list):
            self.smp_name = fx_name(self.fq1[0], fix_pe=True, fix_rep=True, fix_unmap=True)
        self.init_files()
        self.init_fq()
        self.init_index()
        

    def init_files(self):
        # path, files
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        default_dirs = {
            'align_dir': 'align',
            'bam_dir': 'bam_files',
            'bw_dir': 'bw_files',
            'count_dir': 'count',
            'qc_dir': 'qc',
            'report_dir': 'report',
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key
        # files
        align_prefix = os.path.join(self.align_dir, self.smp_name)
        count_prefix = os.path.join(self.count_dir, self.smp_name)
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
            'report_html': os.path.join(self.report_dir, 'HiSeq_report.html'),
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'bw_fwd': self.bw_dir + '/' + self.project_name + '.fwd.bigWig',
            'bw_rev': self.bw_dir + '/' + self.project_name + '.rev.bigWig',
            'align_json': align_prefix+'.align.json',
            
            # quantification files
            'strandness_json': count_prefix+'.strandness.json',
            'count_sens': count_prefix+'.sens.txt',
            'count_anti': count_prefix+'.anti.txt',     
            
            # qc
            'align_scale_json': align_prefix+'.scale.json',
            'genebody_enrich_matrix': os.path.join(self.qc_dir, '05.genebody_enrich.mat.gz'),
            'genebody_enrich_matrix_log': os.path.join(self.qc_dir, '05.genebody_enrich.log'),
            'genebody_enrich_png': os.path.join(self.qc_dir, '05.genebody_enrich.png'),
            'genebody_enrich_cmd': os.path.join(self.qc_dir, '05.genebody_enrich.cmd.sh'),
            'bam_cor_npz': os.path.join(self.qc_dir, '06.bam_cor.npz'),
            'bam_cor_counts': os.path.join(self.qc_dir, '06.bam_cor.counts.tab'),
            'bam_cor_heatmap_png': os.path.join(self.qc_dir, '06.bam_cor.cor_heatmap.png'),
            'bam_cor_pca_png': os.path.join(self.qc_dir, '06.bam_cor.cor_PCA.png')
        }
        self = update_obj(self, default_files, force=True) # key
        check_path([
            self.config_dir, self.align_dir, self.bam_dir, self.bw_dir,
            self.count_dir, self.qc_dir, self.report_dir
        ])


    def init_fq(self):
        # convert to list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        # check message
        c1 = isinstance(self.fq1, list)
        if self.fq2 is None or self.fq2 == 'None':
            self.fq2 = None # convert 'None' -> None; from yaml
            c2e = c2p = True
        else:
            c2e = isinstance(self.fq2, list)
            c2p = check_fx_paired(self.fq1, self.fq2)
        if not all([c1, c2e, c2p]):
            msg = '\n'.join([
                '='*80,
                'Input',
                '{:>14} : {}'.format('fq1', self.fq1),
                '{:>14} : {}'.format('fq2', self.fq2),
                '-'*40,
                'Status',
                '{:>14} : {}'.format('fq1 is str', c1),
                '{:>14} : {}'.format('fq1 is str', c2e),
                '{:>14} : {}'.format('fq1,fq2 paired', c2p),
                '-'*40,
                'Output: {}'.format(all([c1, c2e, c2p])),
                '='*80                
            ])
            print(msg)
            raise ValueError('fq1, fq2 not valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        # update rep_list
        self.rep_list = [
            os.path.join(self.outdir, i) for i in fx_name(self.fq1, fix_pe=True, fix_unmap=True)
        ]
       

    # update: genome_size_file
    def init_index(self):
        # get data from: genome, extra_index
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).get_fasize()
        else:
            raise ValueError('--genome or --extra-index; required')
