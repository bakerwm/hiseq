#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
RNAseq for single fx, single index

"""

import os
import sys
from hiseq.rnaseq.rnaseq_rp import RnaseqRp
from hiseq.align.align_index import AlignIndex
from utils import rnaseq_trim

# from hiseq.rnaseq.utils import rnaseq_trim, rnaseq_align_genome, \
#     rnaseq_align_spikein, rnaseq_call_peak, cal_norm_scale, rnaseq_bam_to_bw, \
#     qc_trim, qc_align, qc_lendist, qc_frip, qc_tss_enrich, qc_genebody_enrich
from hiseq.utils.file import check_path, check_fx_paired, symlink_file, \
    file_abspath, file_prefix, fx_name, Genome
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu, \
    read_hiseq


class RnaseqR1(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = RnaseqR1Config(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)
        Config().dump(self.__dict__, self.config_toml)

        
    def run_pipe(self):
        # symlink raw data
        raw_fq1, raw_fq2 = self.raw_fq_list
        symlink_file(self.fq1, raw_fq1, absolute_path=True)
        symlink_file(self.fq2, raw_fq2, absolute_path=True)
        rnaseq_trim(self.project_dir, 'r1')
#         rnaseq_align_genome(self.project_dir, '_r1')
#         rnaseq_call_peak(self.project_dir, '_r1')
#         rnaseq_bam_to_bw(self.project_dir, '_r1')
#         if isinstance(self.spikein_index, str):
#             rnaseq_align_spikein(self.project_dir, '_r1')
#         # qc
#         qc_trim(self.project_dir, '_r1')
#         qc_align(self.project_dir, '_r1')
#         qc_lendist(self.project_dir, '_r1')
#         qc_frip(self.project_dir, '_r1')
#         qc_tss_enrich(self.project_dir, '_r1')
#         qc_genebody_enrich(self.project_dir, '_r1')


    def run(self):
        # 1. run pipe: RnaseqR1
        self.run_pipe()
        # 2. generate report
        # RnaseqRp(**self.__dict__).run()



class RnaseqR1Config(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'is_mut': False,
            'aligner': 'STAR',            
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'smp_name': None,
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'rRNA_index': None,
            'align_to_rRNA': False,
            'extra_index': None,
            'threads': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'genome_size_file': None,
            'gene_bed': None,
            'gtf': None,
            'keep_tmp': False,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False,
            'extra_para': None,
            'norm_project': None,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_r1'
        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        if not isinstance(self.smp_name, str):
            self.smp_name = fx_name(self.fq1, fix_pe=True)
        self.init_fx()
        self.init_files()
        self.init_index()
        # threads
        self.threads, _ = init_cpu(self.threads, 1)


    def init_fx(self):
        # check message
        c1 = isinstance(self.fq1, str)
        if self.fq2 is None:
            c2 = c3 = True
        else:
            c2 = isinstance(self.fq2, str)
            c3 = check_fx_paired(self.fq1, self.fq2)
        if not all([c1, c2, c3]):
            msg = '\n'.join([
                '='*80,
                'Input',
                '{:>14} : {}'.format('fq1', self.fq1),
                '{:>14} : {}'.format('fq2', self.fq2),
                '-'*40,
                'Status',
                '{:>14} : {}'.format('fq1 is str', c1),
                '{:>14} : {}'.format('fq1 is str', c2),
                '{:>14} : {}'.format('fq1,fq2 paired', c3),
                '-'*40,
                'Output: {}'.format(all([c1, c2, c3])),
                '='*80                
            ])
            print(msg)
            raise ValueError('fq1, fq2 not valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)

    
    # update: genome_size_file    
    def init_index(self):
        # get data from: genome, extra_index
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).get_fasize()
        else:
            raise ValueError('--genome or --extra-index; required')

            

    def init_files(self):
        # dirs
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        default_dirs = {
            'raw_dir': 'raw_data',
            'clean_dir': 'clean_data',
            'align_dir': 'align',
            'spikein_dir': 'spikein',
            'rRNA_dir': 'rRNA',
            'bam_dir': 'bam_files',
            'bg_dir': 'bg_files',
            'bw_dir': 'bw_files',
            'count_dir': 'count',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key
        # files
        default_files = {
            'config_toml': self.config_dir + '/config.toml',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'bw_fwd': self.bw_dir + '/' + self.project_name + '.fwd.bigWig',
            'bw_rev': self.bw_dir + '/' + self.project_name + '.rev.bigWig',
            'strandness_file': self.count_dir + '/strandness.txt',
            'count_sens': self.count_dir \
                + '/' + self.project_name + '.sens.txt',
            'count_anti': self.count_dir \
                + '/' + self.project_name + '.anti.txt',
            'trim_stat_txt': self.clean_dir + '/' + self.project_name \
                + '/' + self.project_name + '.trim.stat',
            'align_scale_txt': self.align_dir + '/' + 'scale.txt',
            'align_flagstat': self.align_dir + '/' \
                + self.smp_name + '.flagstat',
            'align_stat': self.align_dir \
                + '/' + self.smp_name + '.bowtie2.stat',
            'align_toml': self.align_dir \
                + '/' + self.smp_name + '.bowtie2.toml',
            'rRNA_bam': self.rRNA_dir + '/rRNA.bam',
            'rRNA_stat': self.rRNA_dir + '/rRNA.align.stat',
            'rRNA_toml': self.rRNA_dir + '/rRNA.align.toml',
            'spikein_bam': self.spikein_dir + '/spikein.bam',
            'spikein_scale_txt': self.spikein_dir + '/scale.txt',
            'spikein_flagstat': self.spikein_dir + '/spikein.flagstat',
            'spikein_stat': self.spikein_dir + '/spikein.align.stat',
            'spikein_toml': self.spikein_dir + '/spikein.align.toml',
            'trim_summary_toml': self.qc_dir + '/00.trim_summary.toml',
            'align_summary_toml': self.qc_dir + '/01.alignment_summary.toml',
            'genebody_enrich_matrix': self.qc_dir \
                + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir \
                + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_npz': self.qc_dir + '/06.bam_cor.npz',
            'bam_cor_counts': self.qc_dir + '/06.bam_cor.counts.tab',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png'
        }
        self = update_obj(self, default_files, force=True) # key
        self.raw_fq_list = [self.raw_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            raw_fq2 = self.raw_dir + '/' + os.path.basename(self.fq2)
        else:
            raw_fq2 = None
        self.raw_fq_list.append(raw_fq2)
        self.clean_fq_list = [self.clean_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            clean_fq2 = self.clean_dir + '/' + os.path.basename(self.fq2)
        else:
            clean_fq2 = None
        self.clean_fq_list.append(clean_fq2)
        # create dirs
        dir_list = [
            self.config_dir, self.raw_dir, self.clean_dir, self.align_dir,
            self.spikein_dir, self.rRNA_dir, self.bam_dir, self.bg_dir,
            self.bw_dir, self.count_dir, self.qc_dir, self.report_dir
        ]
        check_path(dir_list)



        

        
