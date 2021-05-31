#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
RNAseq for single fx, single index

"""

import os
import sys
from hiseq.rnaseq.rnaseq_rp import RnaseqRp
from hiseq.align.align_index import AlignIndex
from hiseq.rnaseq.utils import rnaseq_trim, rnaseq_align_genome, \
    rnaseq_align_spikein, rnaseq_call_peak, cal_norm_scale, rnaseq_bam_to_bw, \
    qc_trim, qc_align, qc_lendist, qc_frip, qc_tss_enrich, qc_genebody_enrich
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
        rnaseq_align_genome(self.project_dir, '_r1')
        rnaseq_call_peak(self.project_dir, '_r1')
        rnaseq_bam_to_bw(self.project_dir, '_r1')
        if isinstance(self.spikein_index, str):
            rnaseq_align_spikein(self.project_dir, '_r1')
        # qc
        qc_trim(self.project_dir, '_r1')
        qc_align(self.project_dir, '_r1')
        qc_lendist(self.project_dir, '_r1')
        qc_frip(self.project_dir, '_r1')
        qc_tss_enrich(self.project_dir, '_r1')
        qc_genebody_enrich(self.project_dir, '_r1')


    def run(self):
        # 1. run pipe: RnaseqR1
        self.run_pipe()
        # 2. generate report
        RnaseqRp(**self.__dict__).run()



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
        self.threads, self.parallel_jobs = init_cpu(
            self.threads,
            self.parallel_jobs)


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

        # raw data
        self.raw_fq_list = [self.raw_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            raw_fq2 = self.raw_dir + '/' + os.path.basename(self.fq2)
        else:
            raw_fq2 = None
        self.raw_fq_list.append(raw_fq2)

        ## clean data
        self.clean_fq_list = [self.clean_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            clean_fq2 = self.clean_dir + '/' + os.path.basename(self.fq2)
        else:
            clean_fq2 = None
        self.clean_fq_list.append(clean_fq2)

        check_path([
            self.config_dir,
            self.raw_dir,
            self.clean_dir,
            self.align_dir,
            self.spikein_dir,
            self.rRNA_dir,
            self.bam_dir,
            self.bg_dir,
            self.bw_dir,
            self.count_dir,
            self.qc_dir,
            self.report_dir])            

    def init_files(self, create_dirs=True):
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        default_dirs = {
            'raw_dir': 'raw_data',
            'clean_dir': 'clean_data',
            'align_dir': 'align',
            'spikein_dir': 'spikein',
            'bam_dir': 'bam_files',
            'bg_dir': 'bg_files',
            'bw_dir': 'bw_files',
            'peak_dir': 'peak',
            'motif_dir': 'motif',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True)
        # files
        default_files = {
            'config_toml': self.config_dir + '/config.toml', # updated
            'report_log': self.report_dir + '/report.log',
            'report_html': self.report_dir + '/HiSeq_report.html',
            'raw_fq1': self.raw_dir + '/' + os.path.basename(self.fq1),
            'raw_fq2': self.raw_dir + '/' + os.path.basename(self.fq2),
            'clean_fq1': self.clean_dir + '/' + os.path.basename(self.fq1),
            'clean_fq2': self.clean_dir + '/' + os.path.basename(self.fq2),

            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bam_rmdup': self.bam_dir + '/' + self.smp_name + '.rmdup.bam',
            'bam_proper_pair': self.bam_dir + '/' + self.smp_name + '.proper_pair.bam',
            'bed': self.bam_dir + '/' + self.project_name + '.bed',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',

            'align_scale_json': self.bam_dir + '/scale.json',
            'trim_stat_txt': self.clean_dir + '/' + self.project_name + '/' + self.project_name + '.trim.stat',
            'trim_stat_json': self.clean_dir + '/' + self.project_name + '/' + self.project_name + '.trim.json',
            'align_flagstat': self.align_dir + '/' + self.smp_name + '.flagstat',
            'align_stat': self.align_dir + '/' + self.smp_name + '.align.stat',
            'align_json': self.align_dir + '/' + self.smp_name + '.align.json',
            'spikein_bam': self.spikein_dir + '/spikein.bam',
            'spikein_scale_json': self.spikein_dir + '/scale.json',
            'spikein_flagstat': self.spikein_dir + '/spikein.flagstat',
            'spikein_stat': self.spikein_dir + '/spikein.align.stat',
            'spikein_json': self.spikein_dir + '/spikein.align.json',

            'trim_summary_json': self.qc_dir + '/00.trim_summary.json',
            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.pdf',
            'frip_json': self.qc_dir + '/03.FRiP.json',
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
        }
        self = update_obj(self, default_files, force=True) # key
        self.raw_fq_list = [self.raw_fq1, self.raw_fq2]
        self.clean_fq_list = [self.clean_fq1, self.clean_fq2]
        # dirs
        dir_list = [
            self.config_dir, self.raw_dir, self.clean_dir, self.align_dir,
            self.spikein_dir, self.bam_dir, self.bg_dir, self.bw_dir,
            self.peak_dir, self.motif_dir, self.qc_dir, self.report_dir,
        ]
        check_path(dir_list, create_dirs=True)






        


class RNAseqR1Config(object):

    def init_files(self):
        """
        Prepare directories, files
        """
        # path, files
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
        # convert to path
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

        # raw data
        self.raw_fq_list = [self.raw_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            raw_fq2 = self.raw_dir + '/' + os.path.basename(self.fq2)
        else:
            raw_fq2 = None
        self.raw_fq_list.append(raw_fq2)

        ## clean data
        self.clean_fq_list = [self.clean_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            clean_fq2 = self.clean_dir + '/' + os.path.basename(self.fq2)
        else:
            clean_fq2 = None
        self.clean_fq_list.append(clean_fq2)

        check_path([
            self.config_dir,
            self.raw_dir,
            self.clean_dir,
            self.align_dir,
            self.spikein_dir,
            self.rRNA_dir,
            self.bam_dir,
            self.bg_dir,
            self.bw_dir,
            self.count_dir,
            self.qc_dir,
            self.report_dir])
     
        
        
        
        
        
        
        
        
        


class RNAseqR1(object):
    """
    RNAseq analysis for single replicate
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = RNAseqR1Config(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)
        Toml(self.__dict__).to_toml(self.config_toml)
        # update genome_size: 
        if self.extra_index:
            self.genome_size = AlignIndex(index=self.extra_index, aligner=self.aligner).index_size()


    # main pipeline #
    def prep_raw(self):
        """
        Cut the reads to specific length, from the 3' end
        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.raw_fq_list
        file_symlink(self.fq1, raw_fq1, absolute_path=True)
        file_symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self):
        """
        Trim reads:
        Trim(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        fq1, fq2 = self.raw_fq_list
        clean_fq1, clean_fq2 = self.clean_fq_list

        # whether to trim or not
        if self.trimmed:
            ## copy files
            file_symlink(fq1, clean_fq1)
            file_symlink(fq2, clean_fq2)
        else:
            # update args
            args_local = self.__dict__.copy()
            args_init = {
                'fq1': fq1,
                'fq2': fq2,
                'outdir': self.clean_dir,
                'library_type': None,
                'len_min': 20,
                'cut_after_trim': '7,-7',
                'cut_to_length': self.cut_to_length,
                'recursive': self.recursive,
                'parallel_jobs': 1 # do not allowed > 1
            }
            args_local.update(args_init)
            trim = Trim(**args_local)
            trim.run()

            ## copy files
            if isinstance(fq2, str):
                file_symlink(trim.clean_fq1, clean_fq1)
                file_symlink(trim.clean_fq2, clean_fq2)
            else:
                file_symlink(trim.clean_fq, clean_fq1)


    def align_rRNA(self):
        """
        Alignment reads to rRNA
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        # !!!! force: bowtie2 !!!!
        args_sp = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.rRNA_dir,
            'aligner': AlignIndex().get_aligner(index=self.rRNA_index),
            'genome_index': self.rRNA_index,
            'threads': self.threads,
            'extra_para': self.extra_para,
            'keep_tmp': self.keep_tmp,
            'overwrite': self.overwrite
        }
        # args_local.update(args_sp)
        args_local = args_sp
        align = Align(**args_local)

        # output
        tmp0 = listfile(self.rRNA_dir, '*.bam', recursive=True)
        bam = tmp[0] if len(tmp0) > 0 else None
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            align.run()

        # copy bam to bam_files
        file_symlink(bam, self.rRNA_bam)


    def align_spikein(self):
        """
        Alignment reads, spikein
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_sp = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.spikein_dir,
            'genome_index': self.spikein_index,
            'threads': self.threads,
            'extra_para': self.extra_para,
            'keep_tmp': self.keep_tmp,
            'overwrite': self.overwrite
        }
        args_local.update(args_sp)
        align = Align(**args_local)

        # output
        tmp0 = listfile(self.spikein_dir, '*.bam', recursive=True)
        bam = tmp[0] if len(tmp0) > 0 else None
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            align.run()

        # copy bam to bam_files
        file_symlink(bam, self.spikein_bam)

        self.spikein_scale = self.cal_norm_scale(self.spikein_bam)
        with open(self.spikein_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.spikein_scale))


    def align_genome(self):
        """
        Alignment PE reads to reference genome, using bowtie2
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.align_dir,
            'aligner': self.aligner,
            'genome_index': self.genome_index,
            'extra_para': self.extra_para,
            'keep_tmp': self.keep_tmp,
            'overwrite': self.overwrite
        }
        args_local.update(args_init)
        # fix, extra_index alignment
        if args_local.get('genome_index', None) == args_local.get('extra_index', None):
            args_local['extra_index'] = None
        # align = Alignment(**args_local)
        align = Align(**args_local)
        align.run()

        # output
        tmp0 = listfile(self.align_dir, '*.bam', recursive=True)
        bam = tmp0[0] if len(tmp0) > 0 else None
        file_symlink(bam, self.bam)

        # align.toml
        tmp1 = listfile(self.align_dir, '*align.toml', recursive=True)
        if len(tmp1) > 0:
            file_symlink(tmp1[0], self.align_toml)

#         self.align_scale = self.cal_norm_scale(self.bam)
#         with open(self.align_scale_txt, 'wt') as w:
#             w.write('{:.4f}\n'.format(self.align_scale))
        self.align_scale = self.get_norm_scale()


    def get_norm_scale(self):
        """
        Option-1: from mapped reads
        Option-2: from another project, RNAseqRx_dir, with the same design;
        for example:
        norm TE by genome mapped reads        
        """
        try:
            p = RNAseqReader(self.norm_project)
            if p.hiseq_type.endswith('_rx'):
                r1_list = p.get_r1_dir()
                for r1 in r1_list:
                    if os.path.basename(r1) == self.smp_name:
                        r1_dir = r1
                        break
                    else:
                        r1_dir = None
                # parse the scale txt
                if isinstance(r1_dir, str):
                    q = RNAseqReader(r1_dir)
                    q_scale = q.args.get('align_scale_txt')
                    with open(q_scale) as r:
                        out = eval(r.readline().strip())
                    # link to file
                    file_symlink(q_scale, self.align_scale_txt)
            else:
                out = None
        except:
            log.warning('could not read scale from project: {}'.format(self.norm_project))
            log.warning('force norm by current project: {}'.format(self.project_dir))
            out = None
        if not isinstance(out, float):
            try:
                ## calculate norm scale
                if file_exists(self.align_scale_txt):
                    with open(self.align_scale_txt) as r:
                        out = eval(r.readline().strip())
                else:
                    out = self.cal_norm_scale(self.bam)
            except:
                log.error('force scale=1, failed: {},'.format(self.align_scale_txt))
                out = 1
            with open(self.align_scale_txt, 'wt') as w:
                w.write('{:.4f}\n'.format(out))
        # final
        return out
            
            
    def get_bam(self, x):
        """
        Get the align bam file, from x dicectory
        align_dir/A, B, C, ...
        Add prefix/ 01, 02, ...
        """
        bamlist = listfile(x, '*.bam', recursive=True)
        return bamlist.pop() if len(bamlist) > 0 else None # last one


    def bam_to_bw(self):
        """
        Convert bam to bigWig
        """
        args = {
            'bam': self.bam,
            'outdir': self.bw_dir,
            'binsize': self.binsize,
            'strandness': 12,
            'genome': self.genome,
            'scaleFactor': self.align_scale,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size,
        }
        Bam2bw(**args).run()
    
    
    def bam_to_bg(self):
        """
        Convert bam to bedgraph

        norm
        bedtools genomecov -bg -scale $scale_factor -ibam bam > bg
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedtools')),
            'genomecov -bg -scale {}'.format(self.align_scale),
            '-ibam {}'.format(self.bam),
            '| sort -k1,1 -k2,2n > {}'.format(self.bg)
            ])

        if file_exists(self.bg) and not self.overwrite:
            log.info('bam_to_bg() skipped, file exists:{}'.format(self.bg))
        else:
            cmd_shell = self.bg_dir + '/cmd.sh'
            with open(cmd_shell, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bam_to_bg() failed')


    def bg_to_bw(self):
        """
        Create bigWig
        bedgraph to bigWig
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedGraphToBigWig')),
            '{} {} {}'.format(self.bg, self.genome_size_file, self.bw)
            ])

        if file_exists(self.bw) and not self.overwrite:
            log.info('bg_to_bw() skipped, file exists:{}'.format(self.bw))
        else:
            cmd_shell = self.bw_dir + '/cmd.sh'
            with open(cmd_shell, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bg_to_bw() failed')


    def fc_count(self):
        """
        Run FeatureCounts for the bam file
        """
        args = self.__dict__.copy()

        # check library type
        s = HiSeqLibType(bam=self.bam, gtf=self.gtf).run()

        if s == 1:
            strand_sens = 1
            strand_anti = 2
        elif s == 2:
            strand_sens = 2
            strand_anti = 1
        else:
            strand_sens = strand_anti = 0

        # args for local run
        args_local = {
            'gtf': self.gtf,
            'bam_list': self.bam,
            'outdir': self.count_dir,
            'threads': self.threads,
            'overwrite': self.overwrite
        }

        # for sense
        args_sens = {
            'prefix': os.path.basename(self.count_sens),
            'strandness': strand_sens
        }
        args_sens.update(args_local)
        fc_sens = FeatureCounts(**args_sens)
        fc_sens.run()

        # for antisense
        args_anti = {
            'prefix': os.path.basename(self.count_anti),
            'strandness': strand_anti
        }
        args_anti.update(args_local)
        fc_anti = FeatureCounts(**args_anti)
        fc_anti.run()

        # save strandness to file
        msg = '\n'.join([
            'sens: {}'.format(strand_sens),
            'anti: {}'.format(strand_anti)
        ])
        with open(self.strandness_file, 'wt') as w:
            w.write(msg + '\n')


    # quality wildtype #
    def cal_norm_scale(self, bam, norm=1000000):
        """
        scale factor

        Bam().count_reads()
        """
        bam_o = Bam(bam)
        is_pe = bam_o.isPaired()
        n_map = bam_o.getNumberOfAlignments()
        if is_pe:
            n_map = n_map/2

        if n_map > 0:
            n_scale = norm/n_map
        else:
            log.error('no mapped reads detected')
            n_scale = 1

        return n_scale


    def qc_trim_summary(self):
        """
        reads trim off

        # version-1
        #sample wildtype   output  percent
        fq_rep1      2234501 2234276 99.99%

        # version-2
        #name   total   too_short       dup     too_short2      clean   percent
        """
        if file_exists(self.trim_stat_txt):
            with open(self.trim_stat_txt, 'rt') as r:
                lines = r.readlines()
            lines = [i for i in lines if not i.startswith('#')]
            line = lines.pop() # last line
            tabs = line.strip().split('\t')
            fqname, n_wildtype = tabs[:2]
            n_output, n_pct = tabs[-2:]
            d = {
                'id': fqname,
                'input': n_wildtype,
                'output': n_output,
                'out_pct': n_pct,
                'rm_pct': round(100 - float(n_pct.strip()), 2)
            }
        else:
            d = {
                'id': self.smp_name,
                'input': 0,
                'output': 0,
                'out_pct': 100,
                'rm_pct': 0
            }
        Toml(d).to_toml(self.trim_summary_toml)


    def qc_align_summary(self):
        """
        Save summary to file
        """
        # align to genome
        align = Toml().from_toml(self.align_toml)

        # spikein
        if file_exists(self.spikein_bam):
            n_spikein = Bam(self.spikein_bam).getNumberOfAlignments()
            if is_pe:
                n_spikein = n_spikein/2
        else:
            n_spikein = 0

        # rRNA
        if file_exists(self.rRNA_bam):
            n_rRNA = Bam(self.spikein_bam).getNumberOfAlignments()
            if is_pe:
                n_rRNA = n_rRNA/2
        else:
            n_rRNA = 0

        align['nodup'] = 0
        align['spikein'] = n_spikein
        align['chrM'] = n_rRNA

        # save to file
        Toml(align).to_toml(self.align_summary_toml)


    def qc_genebody_enrich(self):
        """
        Calculate the TSS enrichment
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'scale-regions',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(self.bw),
            '-o {}'.format(self.genebody_enrich_matrix),
            '-b 2000 -a 2000 --regionBodyLength 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.genebody_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.genebody_enrich_matrix),
            '-o {}'.format(self.genebody_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.genebody_enrich_png) and not self.overwrite:
            log.info('qc_genebody_enrich() skipped, file exists: {}'.format(
                self.genebody_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.genebody_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_genebody_enrich() failed, see: ')


    def pipe_genome(self):
        """
        Run for genome
        """
        self.prep_raw()
        self.trim()

        if isinstance(self.spikein_index, str):
            self.align_spikein() # run alignment

        if isinstance(self.rRNA_index, str):
            self.align_rRNA()

        self.align_genome()
        self.bam_to_bw()
#         if self.genome:
#             self.bam_to_bw()
#         else:
#             self.bam_to_bg()
#             self.bg_to_bw()

        # count reads
        if file_exists(self.count_sens) and not self.overwrite:
            log.info('fc_count() skipped, file exists: {}'.format(
                self.count_sens))
        else:
            self.fc_count()


    def pipe_qc(self):
        """
        Quality control
        """
        self.qc_trim_summary()
        self.qc_align_summary()
        self.qc_genebody_enrich()


    def pseudo_align(self):
        """
        pseudo_align_dir: salmon.
        Quantification + DEanalysis 
        aligner: salmon, kallisto
        deseq:   DESeq2
        
        required arguments:
        aligner: salmon|kallisto
            fq1
            fq2
            outdir
            index
            'index_name',
            'smp_name',
            'threads',
            'overwrite'
        """
        pass
        
        
    def report(self):
        """
        Create report for one sample
        html
        """
        RNAseqRp(self.project_dir).report()


    def run(self):
        """
        Run all steps for RNAseq pipeline
        """
        self.pipe_genome()
        self.pipe_qc()
        self.report()
