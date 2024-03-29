#!/usr/bin/env python3

"""
piRNA analysis (small RNAseq) 


date: 2020-12-27

1. support collapsed reads, counting, RNA seqs, reads


date: 2020-12-23
Author: Ming Wang

# flowchart of piRNA analysis
1. remove structural RNAs (uniqe + multi)
2. remove miRNAs (unique + multi) 
3. remove reads not in [23-29] nt 
4. collapse reads: consider 1-23nt only, allow 1-2 at 3' differ (2019-11-26)
5. split into 4 groups: (1U+/-,10A+/-;)
6. map to TE consensus, (unique + multi), only 1U_not_10A
7. map to genome (unique + multi) 
Functions:
rename fastq reads: piR0000001-0000001
piR (piRNA), (piRNA number) {reads number}
mapping, unique + multiple/unique


version: 2020-07-28
update:
1. remove temp fastq files
2. gzip fastq files

version: 2020-07-25
update:
1. collapse reads 


date: 2020-07-23
in brief:

1. remove 3' adapters  
2. remove structural RNAs (tRNAs, rRNAs) (unique + multiple)     
3. filt by length, 23-29 nt   
4. collapse reads, only compare 1-23 nt (*)  
   collapse reads, (regular)  
   trim reads to 23nt from 3' end  

5. map to TE consensus (unique, multiple, unique + multiple)   
6. map to piRNA clusters (unique, multiple, unique + multiple)   
7. map to genome (not-TE, no-piRNAcluster)    
8. Overall, map to genome

"""

import os
import sys
import re
import gzip
import fnmatch
import binascii
import shutil
import logging
import argparse
import pybedtools
import pysam
import pandas as pd
import numpy as np
from multiprocessing import Pool
from xopen import xopen
from hiseq.trim.trimmer import Trim
from piRNA_pipe_utils import *


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


def get_args():
    """
    require arguments:
    -i fq
    -o outdir 
    --index-rRNA  
    --index-miRNA 
    --index-te  
    --index-piRNAcluster 
    --index-genome 
    --ref-te
    --ref-piRNAcluster
    ...
    """
    parser = argparse.ArgumentParser(description='piRNA pipe')
    parser.add_argument('-i', '--fq', required=True,
        help='fastq file') 
    parser.add_argument('-o', '--outdir', required=True,
        help='directory to save results')
    parser.add_argument('-w', '--workflow', type=int, default=1,
        help='workflow, 1:4, 1:TE->piRC->Genome; 2:TE->Genome; \
        3:piRC->TE->Genome; 4:piRC->Genome;  default: [1]')
    parser.add_argument('-c', '--collapse', action='store_true',
        help='collapse fastq reads')
    parser.add_argument('-t', '--trimmed', action='store_true',
        help='Input file was clean data, no need to trim adapters')
    parser.add_argument('--subject', default=None,
        help='The target file for overlap')
    parser.add_argument('-p', '--threads', type=int, default=4,
        help='Number of threads to run')
    parser.add_argument('-j', '--parallel-jobs', type=int, default=1,
        dest='parallel_jobs',
        help='number of jobs to run in parallel')
    parser.add_argument('-g', '--genome', default='dm6',
        help='reference genome, default: [dm6]')
    parser.add_argument('-ir', '--index-rRNA', default='smRNA',
        help='bowtie_index for rRNA, default: ')
    parser.add_argument('-im', '--index-miRNA', default='hairpin',
        help='bowtie_index for miRNA (hairpin), default:')
    parser.add_argument('-ite', '--index-te', default='transposon',
        help='bowtie_index for te (TE consensus), default:') 
    parser.add_argument('-ip', '--index-piRNAcluster', default='piRNAcluster',
        help='bowtie_index for piRNA cluster, default:') 
    parser.add_argument('-ig', '--index-genome', default='genome',
        help='bowtie_index for genome, default:') 
    parser.add_argument('-te', '--te-fa', default='te.fa',
        help='TE sequence in fasta format') 
    parser.add_argument('-pi', '--piRNAcluster-fa', default='piRNAcluster.fa',
        help='piRNA cluster sequence in fasta format')

    return parser.parse_args()


class pipeConfig(object):
    """Prepare data for piRNA analysis

    1. index: structureRNA, miRNA, TE, piRNA_cluster, genome, ...
    2. outputdir
    3. pipeline version
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = {
            'outdir': None,
            'smp_name': None,
            'genome': 'dm6',
            'threads': 4,
            'parallel_jobs': 1,
            'collapse': True,
            'workflow': 1,
            'trimmed': True,
            'force_overlap': False
        }
        self = update_obj(self, args_local, force=False)
        if self.outdir == None:
            self.outdir = str(pathlib.Path.cwd())
        check_path(self.outdir)

        if self.smp_name == None:
            self.smp_name = fq_name(self.fq, pe_fix=True)
        self.fq_name = self.smp_name + '.fq.gz' # output file

        self.outdir = os.path.join(self.outdir, self.smp_name)
        if self.force_overlap and file_exists(self.subject):
            self.outdir = os.path.join(self.outdir, 'overlap')
        else:
            self.subject = None
        self.outdir = file_abspath(self.outdir)

        # update files
        self.init_dirs()
        self.init_files()


    def init_files(self):
        index_dir = ''.join([
            '/home/wangming/data/genome/',
            '{}'.format(self.genome),
            '/bowtie_index'
            ])
        args_f = {
            'index_dir': index_dir,
            'smRNA_index': index_dir + '/smRNA',
            'miRNA_index': index_dir + '/hairpin',
            'te_index': index_dir + '/te',
            'piRC_index': index_dir + '/piRNA_cluster',
            'genome_index': index_dir + '/' + self.genome,
            'config_toml': self.config_dir + '/config.toml'
        }
        self = update_obj(self, args_f, force=True)


    def init_dirs(self):
        """The directory structure
        00.total
        01.collapse
        02.smRNA
        03.miRNA
        04.size_select
        05.TE (could be: piRNA_cluster, ...)
        06.genome
        07.unmap
        08.stat
        """
        arg_dirs = {
            'config_dir': self.outdir + 'config',
            'raw_dir': self.outdir + '/00.raw_data', 
            'overlap_dir': self.outdir + '/01.overlap',
            'clean_dir': self.outdir + '/02.clean_data',
            'collapse_dir': self.outdir + '/03.collapse',
            'smRNA_dir': self.outdir + '/04.smRNA',
            'miRNA_dir': self.outdir + '/05.miRNA',
            'size_dir': self.outdir + '/06.size_select',
            'size_ex_dir': self.outdir + '/06.size_exclude',
            'te_dir': self.outdir + '/07.te', # TE/piR_C/genome
            'piRC_dir': self.outdir + '/08.piRNA_cluster', # piR_C (optional)
            'genome_dir': self.outdir + '/09.genome', # genome (optional)
            'unmap_dir': self.outdir + '/10.unmap',
            'stat_dir': self.outdir + '/11.stat',
            'report_dir': self.outdir + '/12.report'
        }
        self = update_obj(self, arg_dirs, force=True)

        check_path([
            self.raw_dir,
            self.overlap_dir,
            self.clean_dir,
            self.collapse_dir,
            self.smRNA_dir,
            self.miRNA_dir,
            self.size_dir,
            self.size_ex_dir,
            self.te_dir,
            self.piRC_dir,
            self.genome_dir,
            self.unmap_dir,
            self.stat_dir,
            self.report_dir
            ])


class pipe(object):
    """Run piRNA analysis pipeline

    pipeline for piRNA analysis

    length distribution
    alignment stat
    1U, 10A stat
    TE mapping
    piRNA cluster mapping

    overlap TE/piRNA_cluster
    
    Overview:
    1. remove structural RNAs (tRNAs, rRNAs) (unique + multiple)     
    2. filt by length, 23-29 nt   
    3. collapse reads, only compare 1-23 nt (*)  
       collapse reads, (regular)  
       trim reads to 23nt from 3' end
    4. map to TE consensus (unique, multiple, unique + multiple)   
    5. map to piRNA clusters (unique, multiple, unique + multiple)   
    6. map to genome (not-TE, no-piRNAcluster)    
    7. overall, map to genome
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args_local = pipeConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)


    def prep_raw(self, fq):
        """Copy/symlink raw fastq files
        split files into 1U,10A, add stat
        """
        log.info('00.Copy raw data')
        fq_raw = os.path.join(self.raw_dir, self.fq_name)
        file_symlink(fq, fq_raw)
        self.fx_not_empty(fq_raw, 'prep_raw', exit=True)
        return fq_raw


    def run_overlap(self, fq):
        """Check single fastq, overlap with subject
        Input fastq
        """
        log.info('01.Overlap with other RNAs')
        fq_overlap = os.path.join(self.overlap_dir, self.fq_name)
        if self.subject:
            overlap_fq(fq, self.subject, self.overlap_dir)
        else:
            file_symlink(fq, fq_overlap)
        self.fx_not_empty(fq_overlap, 'run_overlap', exit=True)
        return fq_overlap


    def run_trim(self, fq, trimmed=False):
        """Cut the adapters from the 3' end
        Trim reads:
        Trim(fq1, outdir, fq2, cut_after_trim='9,-6').run()
        if not:
            copy/links
        """
        log.info('02.Trim adapters')
        fq_clean = os.path.join(self.clean_dir, self.fq_name)
        if trimmed:
            file_symlink(fq, fq_clean)
        else:
            args_local = {
                'fq1': fq,
                'outdir': self.clean_dir,
                'library_type': None,
                'len_min': 15,
                'cut_to_length': 0,
                'recursive': False,
                'parallel_jobs': 1
            }
            trim = Trim(**args_local)
            trim.run()
            file_symlink(trim.clean_fq, fq_clean)
        self.fx_not_empty(fq_clean, 'trim', exit=True)
        return fq_clean


    def run_collapse(self, fq):
        """Collapse fastq files
        save read count in id line
        """
        log.info('03.Collapse fastq')
        fq_collapse = os.path.join(self.collapse_dir, self.fq_name)
        if self.collapse:
            out_fq = collapse_fx(fq, fq_collapse)
            file_symlink(out_fq, fq_collapse)
        else:
            file_symlink(fq, fq_collapse)
        self.fx_not_empty(fq_collapse, 'run_collapse', exit=True)
        return fq_collapse


    def run_smRNA(self, fq):
        """Map reads to structural RNAs, remove map reads
        unique, multiple
        """
        log.info('04.Map to small RNA')
        # alignment, unique + multiple, k=1
        fq_smRNA, fq_not_smRNA = pipe_align(fq, self.smRNA_index, self.smRNA_dir, 
            self.threads, self.genome, remove_1u10a=False)
        self.fx_not_empty(fq_not_smRNA, 'run_smRNA', exit=True)
        return fq_not_smRNA


    def run_miRNA(self, fq):
        """Map reads to miRNA, remove
        unique, multiple
        """
        log.info('05.Map to miRNA')
        fq_miRNA, fq_not_miRNA = pipe_align(fq, self.miRNA_index, self.miRNA_dir, 
            self.threads, self.genome, remove_1u10a=False)
        self.fx_not_empty(fq_not_miRNA, 'run_miRNA', exit=True)
        return fq_not_miRNA


    def run_size_select(self, fq):
        """Split fastq files by size
        23-29: piRNA candidates
        """
        log.info('06.Size select')        
        fq_size_in = os.path.join(self.size_dir, self.fq_name)
        fq_size_ex = os.path.join(self.size_ex_dir, self.fq_name)
        fq_in, fq_ex = splilt_fq_size(fq, outdir=self.size_dir,
            min=23, max=29, gzipped=True)
        file_symlink(fq_in, fq_size_in)
        file_symlink(fq_ex, fq_size_ex)
        self.fx_not_empty(fq_size_in, 'run_size_select', exit=True)
        return fq_size_in


    def run_TE(self, fq):
        """Map reads to TE consensus
        """
        log.info('07.map to TE')
        # fq_te = os.path.join(self.te_dir, self.fq_name)
        fq_te, fq_not_te = pipe_align(fq, self.te_index, self.te_dir, 
            self.threads, self.genome, unique_multi='all')
        self.fx_not_empty(fq_not_te, 'run_TE', exit=True)
        return fq_not_te


    def run_piRC(self, fq):
        """Map reads to piRNA cluster
        mapping
        """
        log.info('08.Map to piRNA cluster')
        # fq_piRC = os.path.join(self.piRC_dir, self.fq_name)
        fq_piRC, fq_not_piRC = pipe_align(fq, self.piRC_index, self.piRC_dir, 
            self.threads, self.genome, unique_multi='all')
        self.fx_not_empty(fq_not_piRC, 'run_PiRC', exit=False)
        return fq_not_piRC


    def run_genome(self, fq):
        """map reads to genome
        non-TE reads to genome
        """
        log.info('09.Map to genome')
        # fq_genome = os.path.join(self.genome_dir, self.fq_name)
        fq_unmap = os.path.join(self.unmap_dir, self.fq_name)
        fq_genome, fq_not_genome = pipe_align(fq, self.genome_index, self.genome_dir, 
            self.threads, self.genome, unique_multi='all')
        file_symlink(fq_not_genome, fq_unmap)
        self.fx_not_empty(fq_not_genome, 'run_genome', exit=False)
        return fq_not_genome


    ################################
    # if reads empty:
    # Break pipe
    # prep_raw, trim, run_collapse, run_smRNA, run_miRNA, run_size_select
    #
    # Skip pipe
    # run_TE, run_piRC, run_genome
    ################################
    def fx_not_empty(self, fx, name=None, exit=False):
        """Check the fx empty or not;
        count reads
        """
        chk = check_file(fx, emptycheck=True)
        if not chk:
            log.error('Stop at {}(), file not exists: {}'.format(name, fx))
            if exit:
                raise Exception('exit pipeline at {}().'.format(name))


    ################################
    # return files, folders
    # group=overlap,collapse,...
    ################################
    def get_sub_dir(self, u1a10=False, group=None):
        """Return the working directories
        group: overlap
        """
        # level-1 00.raw_dat to 09.unmap
        sub_dir = listfile(self.outdir, include_dir=True)
        # add group
        if isinstance(group, str):
            sub_dir = [os.path.join(i, group) for i in sub_dir]
        # u1a10
        if u1a10:
            sub_dir = [os.path.join(i, '1U_10A') for i in sub_dir]
        # existence
        sub_dir = [i for i in sub_dir if os.path.exists(i)]
        return sub_dir


    def get_sub_fx(self, u1a10=False, group=None):
        """Return the 00.raw_data to 09.unmap fastq files
        group: overlap, collapse
        """
        sub_fx = []
        sub_dirs = self.get_sub_dir(u1a10, group)
        for d in sub_dirs:
            fx = listfile(d, '*.gz', recursive=False)
            if len(fx) > 0:
                sub_fx.extend(fx)
        return sub_fx


    def get_sub_stat(self, u1a10=False, group=None):
        """Return the 00.raw_data to 09.unmap stat file: fq_stat.toml
        group: overlap, collapse
        """
        stat_list = []
        sub_dirs = self.get_sub_dir(u1a10, group)
        for d in sub_dirs:
            s = listfile(d, 'fq_stat.toml', recursive=False)
            if len(s) > 0:
                stat_list.extend(s)
        return stat_list


    ################################
    # Check overlap between subject file
    # group=overlap,collapse,...
    ################################
    # # worker: for single fx
    # def run_overlap_single(self, i):
    #     """Check single fastq, overlap with subject
    #     Input fastq
    #     """
    #     fx = self.overlap_fq_list[i]
    #     outdir = os.path.join(os.path.dirname(fx), 'overlap')
    #     overlap_fq(fx, self.subject, outdir)


    # def run_overlap(self):
    #     """Check fastq files overlap with subject
    #     Create overlap dir
    #     output overlap files
    #     """
    #     # all fastq files, to compare
    #     self.overlap_fq_list = self.get_sub_fx(u1a10=False) # level-1: 00.raw_data 
    #     self.overlap_fq_list = [i for i in self.overlap_fq_list \
    #         if check_file(i, emptycheck=True)]
    #     # multi threads
    #     if len(self.overlap_fq_list) > 1 and self.parallel_jobs > 1:
    #         with Pool(processes=self.parallel_jobs) as pool:
    #             pool.map(self.run_overlap_single, range(len(self.overlap_fq_list)))
    #     else:
    #         for i in range(len(self.overlap_fq_list)):
    #             self.run_overlap_single(i)


    ################################
    # stat: 1U_10A
    ################################
    # worker: for single fx
    def run_1u10a_single(self, i):
        """Split fx, by 1U, 10A
        input: index
        """
        fx = self.u1a10_fq_list[i]
        split_fq_1u10a(fx, remove=False)


    def run_1u10a(self):
        """Check the 1U10A content for fastq files
        level-1: 00.raw_data
        level-2: 00.raw_data/overlap/
        """
        fq_a = self.get_sub_fx(u1a10=False, group=None) # level-1
        fq_b = self.get_sub_fx(u1a10=False, group='overlap')
        self.u1a10_fq_list = fq_a + fq_b
        # remove empty files
        self.u1a10_fq_list = [i for i in self.u1a10_fq_list \
            if check_file(i, emptycheck=True)]
        # multi threads
        if len(self.u1a10_fq_list) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_1u10a_single, range(len(self.u1a10_fq_list)))
        else:
            for i in range(len(self.u1a10_fq_list)):
                self.run_1u10a_single(i)


    ################################
    # count reads: stat_fq
    ################################
    # worker: for single dir
    def run_count_fx_single(self, i):
        """Stat fq for single file/dir
        wrap fq_stat.toml for dir
        """
        i_dir = self.count_dir_list[i]
        count_fx_dir(i_dir, collapsed=True) # save to fname.fq_stat.toml


    def run_count_fx(self):
        """Count reads, RNA seqs in each categories
        level-1: 00.raw_data
        level-1: 00.raw_data/1U_10A
        level-2: 00.raw_data/overlap
        level-2: 00.raw_data/overlap/1U_10A
        """
        dir_a = self.get_sub_dir(u1a10=False, group=None) # level-1
        dir_b = self.get_sub_dir(u1a10=True, group=None) # level-1,u1a10
        dir_c = self.get_sub_dir(u1a10=False, group='overlap') # level-2
        dir_d = self.get_sub_dir(u1a10=True, group='overlap') # level-2,u1a10
        self.count_dir_list = dir_a + dir_b + dir_c + dir_d
        # multi threads
        if len(self.count_dir_list ) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_count_fx_single, 
                    range(len(self.count_dir_list)))
        else:
            for i in range(len(self.count_dir_list)):
                self.run_count_fx_single(i)


    def read_fx_toml(self, x):
        """Convert toml to pd.DataFrame
        dict -> pd
        num_reads
        num_seqs
        """
        df = None
        if isinstance(x, str):
            if x.endswith(".toml"):
                d = Toml().from_toml(x)
                if len(d) > 0:
                    df = pd.DataFrame.from_dict(d, 'index')
            else:
                log.error('x is not .toml file')
        else:
            log.error('x is not file')
        return df


    def run_qc_f1(self):
        """Summary the reads in level-1
        f1: level-1: 00.raw_dir
        f2: level-1: 00.raw_data/1U_10A
        v1: level-2: 00.raw_dir/overlap
        v2: level-2: 00.raw_dir/overlap/1U_10A
        """
        dir_list = self.get_sub_dir(u1a10=False, group=None)
        s_list = [i + '/fx_stat.toml' for i in dir_list]
        s_list = [i for i in s_list if os.path.exists(i)]
        s_dir = [os.path.dirname(i) for i in s_list]
        s_frames = [self.read_fx_toml(i) for i in s_list]
        df = pd.concat(s_frames, axis=0)
        df.columns = ['sample', 'num_reads', 'num_seqs']
        df['group'] = [pathlib.Path(i).parts[-1] for i in s_dir]
        df['u1a10'] = 'all'
        df['overlap'] = 'all'
        return df


    def run_qc_f2(self):
        """Summary the reads in level-1
        f1: level-1: 00.raw_dir
        f2: level-1: 00.raw_data/1U_10A
        v1: level-2: 00.raw_dir/overlap
        v2: level-2: 00.raw_dir/overlap/1U_10A
        """
        dir_list = self.get_sub_dir(u1a10=True, group=None)
        s_list = [i + '/fx_stat.toml' for i in dir_list]
        s_list = [i for i in s_list if os.path.exists(i)]
        s_dir = [os.path.dirname(i) for i in s_list]
        s_frames = [self.read_fx_toml(i) for i in s_list]
        # for each file
        s_frames = []
        for i in s_list:
            dfx = self.read_fx_toml(i)
            dfx.columns = ['sample', 'num_reads', 'num_seqs']
            dfx['group'] = pathlib.Path(i).parts[-3]
            s_frames.append(dfx)
        df = pd.concat(s_frames, axis=0)
        df[['sample', 'u1a10']] = df['sample'].str.split('.', 1, expand=True)
        df['overlap'] = 'all'
        return df


    def run_qc_fn(self):
        """Summary the reads in level-1
        f1: level-1: 00.raw_dir
        f2: level-1: 00.raw_data/1U_10A
        v1: level-2: 00.raw_dir/overlap
        v2: level-2: 00.raw_dir/overlap/1U_10A
        """
        df1 = self.run_qc_f1()
        df2 = self.run_qc_f2()
        df3 = pd.concat([df1, df2], axis=0)
        df3 = df3.astype({'num_reads': 'int32', 'num_seqs': 'int32'})
        # reads
        stat_a = os.path.join(self.stat_dir, 'fx_stat.reads.txt')
        df3a = df3.drop(['num_seqs'], axis=1)
        df3a = df3a.pivot_table(columns='u1a10', values='num_reads', 
            index=['sample', 'group', 'overlap'])
        df3a.reset_index(inplace=True)
        df3a = df3a.fillna(0)
        df3a.to_csv(stat_a, index=False)
        # RNAs
        stat_b = os.path.join(self.stat_dir, 'fx_stat.seqs.txt')
        df3b = df3.drop(['num_reads'], axis=1)
        df3b = df3b.pivot_table(columns='u1a10', values='num_seqs', 
            index=['sample', 'group', 'overlap'])
        df3b.reset_index(inplace=True)
        df3b = df3b.fillna(0)
        df3b.to_csv(stat_b, index=False)


    def run_qc_v1(self):
        """Summary the reads in level-1
        f1: level-1: 00.raw_dir
        f2: level-1: 00.raw_data/1U_10A
        v1: level-2: 00.raw_dir/overlap
        v2: level-2: 00.raw_dir/overlap/1U_10A
        """
        dir_list = self.get_sub_dir(u1a10=False, group='overlap')
        s_list = [i + '/fx_stat.toml' for i in dir_list]
        s_list = [i for i in s_list if os.path.exists(i)]
        s_dir = [os.path.dirname(i) for i in s_list]
        s_frames = [self.read_fx_toml(i) for i in s_list]
        df = pd.concat(s_frames, axis=0)
        df.columns = ['sample', 'num_reads', 'num_seqs']
        df['group'] = [pathlib.Path(i).parts[-2] for i in dir_list]
        df['u1a10'] = 'all'
        df['overlap'] = 'overlap'
        df = df.astype({'num_reads':'int32', 'num_seqs': 'int32'})
        return df


    def run_qc_v2(self):
        """Summary the reads in level-1
        f1: level-1: 00.raw_dir
        f2: level-1: 00.raw_data/1U_10A
        v1: level-2: 00.raw_dir/overlap
        v2: level-2: 00.raw_dir/overlap/1U_10A
        """
        dir_list = self.get_sub_dir(u1a10=True, group='overlap')
        s_list = [i + '/fx_stat.toml' for i in dir_list]
        s_list = [i for i in s_list if os.path.exists(i)]
        s_dir = [os.path.dirname(i) for i in s_list]
        # s_frames = [self.read_fx_toml(i) for i in s_list]
        s_frames = []
        for i in s_list:
            dfx = self.read_fx_toml(i)
            dfx.columns = ['sample', 'num_reads', 'num_seqs']
            dfx['group'] = pathlib.Path(i).parts[-4]
            s_frames.append(dfx)
        df = pd.concat(s_frames, axis=0)
        df[['sample', 'u1a10']] = df['sample'].str.split('.', 1, expand=True)
        df['overlap'] = 'overlap'
        df = df.astype({'num_reads':'int32', 'num_seqs': 'int32'})
        return df


    def run_qc_vn(self):
        """Summary the reads in level-1
        f1: level-1: 00.raw_dir
        f2: level-1: 00.raw_data/1U_10A
        v1: level-2: 00.raw_dir/overlap
        v2: level-2: 00.raw_dir/overlap/1U_10A
        """
        df1 = self.run_qc_v1()
        df2 = self.run_qc_v2()
        df3 = pd.concat([df1, df2], axis=0)
        df3 = df3.astype({'num_reads': 'int32', 'num_seqs': 'int32'})
        # reads
        stat_a = os.path.join(self.stat_dir, 'fx_stat.overlap.reads.txt')
        df3a = df3.drop(['num_seqs'], axis=1)
        df3a = df3a.pivot_table(columns='u1a10', values='num_reads', 
            index=['sample', 'group', 'overlap'])
        df3a.reset_index(inplace=True)
        df3a = df3a.fillna(0)
        df3a.to_csv(stat_a, index=False)
        # RNAs
        stat_b = os.path.join(self.stat_dir, 'fx_stat.overlap.seqs.txt')
        df3b = df3.drop(['num_reads'], axis=1)
        df3b = df3b.pivot_table(columns='u1a10', values='num_seqs', 
            index=['sample', 'group', 'overlap'])
        df3b.reset_index(inplace=True)
        df3b = df3b.fillna(0)
        print(df3b)
        df3b.to_csv(stat_b, index=False)


    ################################
    # Post analysis
    ################################
    def run_post_analysis(self):
        """For all post analysis
        overlap
        1U_10A
        stat_fq
        ...
        """
        log.info('10.Post analysis, 1U_10A')
        ##### self.run_overlap() # skipped, time
        self.run_1u10a()
        self.run_count_fx()
        self.run_qc_fn()
        ##### self.run_qc_vn() # skipped, time


    ################################
    # Main port:
    ################################
    def run(self):
        # fq_raw = self.prep_raw(self.fq)
        # fq_overlap = self.run_overlap(fq_raw)
        # fq_clean = self.run_trim(fq_overlap, trimmed=True)
        # fq_collapse = self.run_collapse(fq_clean)
        # fq_not_smRNA = self.run_smRNA(fq_collapse)
        # fq_not_miRNA = self.run_miRNA(fq_not_smRNA)
        # fq_size = self.run_size_select(fq_not_miRNA)
        ## run 00.raw_data -> 05.miRNA
        fq_size = self.run_size_select(
            self.run_miRNA(
                self.run_smRNA(
                    self.run_collapse(
                        self.run_trim(
                            self.run_overlap(
                                self.prep_raw(self.fq)
                            ), 
                            trimmed=self.trimmed
                        )
                    )
                )
            )
        )

        if self.workflow == 1:
            # workflow-1: smRNA->miRNA->size->TE->piRC->genome
            fq_not_genome = self.run_genome(
                self.run_piRC(
                    self.run_TE(fq_size)))
        elif self.workflow == 2:
            # workflow-2: smRNA->miRNA->size->TE->genome
            fq_not_genome = self.run_genome(
                self.run_TE(fq_size))
        elif self.workflow == 3:
            # workflow-1: smRNA->miRNA->size->piRC->TE->genome
            fq_not_genome = self.run_genome(
                self.run_TE(
                    self.run_piRC(fq_size)))
        elif self.workflow == 4:
            # workflow-1: smRNA->miRNA->size->piRC->TE->genome
            fq_not_genome = self.run_genome(
                    self.run_piRC(fq_size))
        else:
            raise Exception('workflow expect [1:4], default: 1')

        ## for statistics
        self.run_post_analysis()


def overlap_fq(query, subject, outdir=None, threads=4):
    """Check overlap between fastq files, by sequence
    using command tool: 
    1. seqkit common -s small.fq.gz big.fq.gz
    2. seqkit grep -s -f <(seqkit seq -s small.fq.gz) big.fq.gz # by seq
    """
    if not isinstance(outdir, str):
        outdir = os.path.dirname(query)
    if not check_file([query, subject], emptycheck=True):
        log.warning('file not exists, or empty: {}, {}'.format(query, subject))
    else:
        check_path(outdir)
        out_fq = os.path.join(outdir, os.path.basename(query))
        out_log = os.path.join(outdir, 'cmd.log')
        out_cmd = os.path.join(outdir, 'cmd.sh')
        ## slower than seqkit common
        # cmd = ' '.join([
        #     'seqkit grep -s -i -j {}'.format(threads),
        #     '-f <(seqkit seq -s {})'.format(query),
        #     '{}'.format(subject),
        #     '2>{}'.format(out_log),
        #     '| pigz -p {} >{}'.format(threads, out_fq),
        # ])
        cmd = ' '.join([
            'seqkit common -s -j {}'.format(threads),
            '{} {}'.format(query, subject),
            '2> {}'.format(out_log),
            '| pigz -p {} > {}'.format(threads, out_fq)
        ])
        with open(out_cmd, 'wt') as w:
            w.write(cmd + '\n')

        if file_exists(out_fq):
            log.info('overlap_fq() skipped, file exists: {}'.format(out_fq))
        else:
            # os.system(cmd)
            run_shell_cmd(cmd)

        return out_fq


def main():
    args = vars(get_args())
    if args['subject']:
        args_local = args.copy()
        args['force_overlap'] = True
        pipe(**args_local).run()
    # force outdir, subject=None
    args_local2 = args.copy()
    args['force_overlap'] = False
    pipe(**args_local2).run()


if __name__ == '__main__':
    main()


## EOF