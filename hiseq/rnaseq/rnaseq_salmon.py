#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
# # building index
# salmon index -t data/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz -i data/dm6_index -k 31
# salmon index -t data/dm6_cdna.fa.gz -i data/dm6_cdna -p 8 -k 31

# quantify
#index="data/dm6_cdna"
# index="data/Dm.BDGP6.28.100.te_piRC.salmon-0.14.1/dm6_100"
# index="data/Dm.BDGP6.28.100.te_piRC.salmon-0.14.1/Dm.BDGP6.28.100.te_piRC"
index="$HOME/data/genome/dm6/salmon_index/genome"
# index="$HOME/data/genome/dm6/salmon_index/te_piRC"
outdir="results/dm6_txome"
for r1 in data/raw_data/*1.fq.gz
do
    r2=${r1/_1.fq/_2.fq}
    prefix=$(basename ${r1/_1.fq.gz})
    echo $prefix
    salmon quant -i $index -l A -p 8 --gcBias -o ${outdir}/$prefix -1 $r1 -2 $r2
done

# DE analysis (DESeq2)
#
# tximport
# tximeta
"""


import os
import sys
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.align.align import Align
from hiseq.rnaseq.rnaseq_rd import RnaseqRd
from hiseq.rnaseq.rnaseq_rp import RnaseqRp
from hiseq.rnaseq.utils import deseq_salmon
from hiseq.utils.file import (
    check_fx_args, check_path, check_fx_paired, symlink_file, file_exists,
    file_abspath, file_prefix, fx_name, Genome
)
from hiseq.utils.utils import (
    log, update_obj, Config, get_date, init_cpu, read_hiseq, list_hiseq_file,
    run_shell_cmd
)
from hiseq.utils.bam import Bam
from hiseq.align.align_index import AlignIndex, fetch_index


class RnaseqSalmonPipe(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RnaseqSalmonPipeConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)


    def run_single_rx(self, i):
        i = {
            k:v for k,v in i.items() if k in \
            ['mut_fq1', 'mut_fq2', 'wt_fq1', 'wt_fq2']
        }
        # update outdir
        wt_name = fx_name(i['wt_fq1'][0], fix_pe=True, fix_rep=True)
        mut_name = fx_name(i['mut_fq1'][0], fix_pe=True, fix_rep=True)
        smp_name = '{}.vs.{}'.format(wt_name, mut_name)
        i.update({
            'build_design': False,
            'design': None,
            'outdir': os.path.join(self.outdir, smp_name)
        })
        if len(self.fq_groups) > 1:
            i['parallel_jobs'] = 1 # force
        args_local = self.__dict__.copy()
        args_local.update(i)
        RnaseqSalmon(**args_local).run()


    def run_multiple_rx(self):
        # load fq groups
        self.fq_groups = Config().load(self.design)
        if len(self.fq_groups) == 0:
            raise ValueError('no data in design: {}'.format(self.design))
        # run multiple in parallel
        if self.parallel_jobs > 1 and len(self.fq_groups) > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_rx, self.fq_groups.values())
        else:
            for i in list(self.fq_groups.values()):
                self.run_single_rx(i)


    def run(self):
        if self.build_design:
            RnaseqRd(**self.__dict__).run()
        else:
            if AlignIndex(self.salmon_index).is_valid():
                self.run_multiple_rx()
            else:
                log.error('salmon_index not valid: {}'.format(self.salmon_index))


class RnaseqSalmonPipeConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'build_design': False,
            'design': None, # required
            'fq_dir': None,
            'mut': None,
            'wt': None,
            'mut_fq1': None,
            'mut_fq2': None,
            'wt_fq1': None,
            'wt_fq2': None,
            'salmon_index': None,
            'genome': None,
            'threads': 4,
            'parallel_jobs': 1
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_ra'
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)
        if not isinstance(self.design, str):
            raise ValueError('--design, expect str, got {}'.format(
                type(self.design).__name__))
        # update salmon_index
        if self.salmon_index is None and isinstance(self.genome, str):
            si = fetch_index(self.genome, aligner='salmon')
            if AlignIndex(si, 'salmon').is_valid():
                self.salmon_index = si


class RnaseqSalmon(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RnaseqSalmonConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)


    def run_salmon_rn(self):
        args = {
            'aligner': 'salmon',
            'fq1': self.mut_fq1 + self.wt_fq1,
            'fq2': self.mut_fq2 + self.wt_fq2 if self.is_paired else None,
            'outdir': self.quant_dir,
            'extra_index': self.salmon_index,
            'threads': self.threads,
        }
        Align(**args).run()


    def run(self):
        # 1. save config
        Config().dump(self.__dict__, self.config_yaml)
        # quant
        self.run_salmon_rn() # multiple
        #deseq
        deseq_salmon(self.project_dir, 'rx')
        #report
        RnaseqRp(self.project_dir).run()


class RnaseqSalmonConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'salmon',
            'mut_fq1': None,
            'mut_fq2': None,
            'wt_fq1': None,
            'wt_fq2': None,
            'salmon_index': None,
            'genome': None, # for deseq analysis
            'smp_name': None, # mut.vs.wt
            'mut_name': None,
            'wt_name': None,
            'outdir': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'trimmed': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_rx' # salmon
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)
        self.mut_fq1, self.mut_fq2, self.mut_is_paired = self.init_fq(
            self.mut_fq1, self.mut_fq2)
        self.wt_fq1, self.wt_fq2, self.wt_is_paired = self.init_fq(
            self.wt_fq1, self.wt_fq2)
        # paired-end
        self.is_paired = self.mut_is_paired and self.wt_is_paired
        self.init_index()
        self.init_name()
        self.init_files()


    def init_fq(self, fq1, fq2):
        # convert to list
        if isinstance(fq1, str):
            fq1 = [fq1]
        if isinstance(fq2, str):
            fq2 = [fq2]
        if not isinstance(fq1, list):
            raise ValueError('fq1 require list, got {}'.format(
                type(fq1).__name__))
        if not check_fx_args(fq1, fq2):
            raise ValueError('fq1, fq2 not valid')
        return (file_abspath(fq1), file_abspath(fq2), check_fx_paired(fq1, fq2))


    def init_index(self):
        # get data from: genome, extra_index
        if not AlignIndex(self.salmon_index, aligner='salmon').is_valid():
            raise ValueError('--salmon-index, not valid: {}'.format(
                self.salmon_index))
        # index_name #updated
        self.index_name = AlignIndex(index=self.salmon_index).index_name()


    def init_name(self):
        """
        mut_name, wt_name, smp_name
        mut_dir, wt_dir
        """
        if not isinstance(self.mut_name, str):
            # fix mut_name
            self.mut_name = fx_name(self.mut_fq1[0], fix_pe=True, fix_rep=True)
        if not isinstance(self.wt_name, str):
            # fix wt_name
            self.wt_name = fx_name(self.wt_fq1[0], fix_pe=True, fix_rep=True)
        if not isinstance(self.smp_name, str):
            # fix smp_name
            self.smp_name = '{}.vs.{}'.format(self.mut_name, self.wt_name)


    def init_files(self):
        # dirs
        self.project_name = self.smp_name
        self.project_dir = self.outdir
        # self.project_dir = os.path.join(self.outdir, self.project_name)
        default_dirs = {
            'config_dir': 'config',
            'quant_dir': 'quant',
            'deseq_dir': 'deseq',
            'enrich_dir': 'enrich',
            'report_dir': 'report',
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key
        # files
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
            'report_html': os.path.join(self.report_dir, 'HiSeq_report.html'),
            'quant_json': os.path.join(self.quant_dir, 'quant.json'),
            'deseq_fix_xls': os.path.join(self.deseq_dir, 'transcripts_deseq2.fix.xls')
        }
        self = update_obj(self, default_files, force=True) # key
        # dirs
        dir_list = [
            self.config_dir, self.quant_dir, self.deseq_dir, self.enrich_dir,
            self.report_dir
        ]
        # update rep dirs
        self.mut_dirs = [os.path.join(self.quant_dir, i) for i in fx_name(
            self.mut_fq1, fix_pe=self.is_paired)]
        self.wt_dirs = [os.path.join(self.quant_dir, i) for i in fx_name(
            self.wt_fq1, fix_pe=self.is_paired)]
        # quant.sf files
        self.init_quant_files()
        check_path(dir_list, create_dirs=True)


    def init_quant_files(self):
        """
        default: outdir/smp_name/index_name/quant.sf
        """
        # mut
        self.mut_quant = []
        for fq1 in self.mut_fq1:
            smp_name = fx_name(fq1, fix_pe=self.is_paired, fix_rep=False)
            sf = os.path.join(self.quant_dir, smp_name, self.index_name, 'quant.sf')
            self.mut_quant.append(sf)
        # wt
        self.wt_quant = []
        for fq1 in self.wt_fq1:
            smp_name = fx_name(fq1, fix_pe=self.is_paired, fix_rep=False)
            sf = os.path.join(self.quant_dir, smp_name, self.index_name, 'quant.sf')
            self.wt_quant.append(sf)
        # save to json
        d ={
            'mut_quant': self.mut_quant,
            'wt_quant': self.wt_quant,
            'genome': self.genome,
        }
        Config().dump(d, self.quant_json)


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. support fastq input',
        '$ hiseq rnaseq_salmon -b -d design.json -r data --wt control --mut treatment',
        '$ hiseq rnaseq_salmon -d design.json -g dm6 -o results',
        '2. for specific index',
        '$ hiseq rnaseq_salmon -d design.json -o results -x path-to-salmon-index',
    ])
    parser = argparse.ArgumentParser(
        prog='rnaseq_salmon',
        description='run rnaseq_salmon',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-b', '--build-design', dest='build_design',
        action='store_true',
        help='generate design.yaml, with --fq-dir, or fq1/fq2')
    parser.add_argument('-d', '--design', required=True,
        help='The file saving fastq files config; generated by rnaseq_rd.py')
    # for design
    parser.add_argument('-r', '--fq-dir', dest='fq_dir', default=None,
        help='directory of fastq files, for --build-design')
    parser.add_argument('--mut', nargs='+', default=None,
        help='keyword of mut fastq file, auto-find read1/2')
    parser.add_argument('--wt', nargs='+', dest='wt', default=None,
        help='keyword of wt fastq file, auto-find read1/2')
    parser.add_argument('--se', dest='as_se', action='store_true',
        help='choose only fastq1 for PE reads')
     # details
    parser.add_argument('--mut-fq1', nargs='+', dest='mut_fq1', default=None,
        help='filepath or keyword of mut fastq file, read1 of PE')
    parser.add_argument('--mut-fq2', nargs='+', dest='mut_fq2', default=None,
        help='filepath or keyword of mut fastq file, read2 of PE')
    parser.add_argument('--wt-fq1', nargs='+', dest='wt_fq1',
        default=None,
        help='filepath or keyword of wt fastq file, read1 of PE')
    parser.add_argument('--wt-fq2', nargs='+', dest='wt_fq2',
        default=None,
        help='filepath or keyword of wt fastq file, read2 of PE')

    # required arguments
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-x', dest='salmon_index', required=False,
        help='The path to salmon index')
    parser.add_argument('-g', '--genome', default=None,
        # choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg38, mm10, default: hg38')

    # optional arguments - 0
    parser.add_argument('-mn', '--mut-name', dest='mut_name', default=None,
        help='Name of mutant samples')
    parser.add_argument('-wn', '--wt-name', dest='wt_name', default=None,
        help='Name of wildtype/control samples')
    parser.add_argument('-n', '--smp-name', dest='smp_name', default=None,
        help='Name of the samples, default: from fq1 prefix')

    # optional arguments - 1
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('-p', '--threads', default=4, type=int,
        help='Number of threads for each job, default [4]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    return parser


def main():
    args = vars(get_args().parse_args())
    RnaseqSalmonPipe(**args).run()


if __name__ == '__main__':
    main()

#
