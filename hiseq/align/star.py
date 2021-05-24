#!/usr/bin/env python3

"""
Align reads to index, using STAR

Basic usage: 

STAR 
"""

import os
import sys
import re
import shutil
import argparse
from Levenshtein import distance
from hiseq.utils.seq import Fastx
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd
from hiseq.utils.file import check_fx, check_file, check_path, file_abspath, file_exists, fx_name, symlink_file, remove_file
from hiseq.align.align_index import AlignIndex
from hiseq.align.utils import check_fx_args, AlignReader


def parse_star(x):
    """Wrapper STAR log
    for PE: only proper paired reads

    *final.Log.out, (changed to *.log, in this script)

                                 Started job on |       Sep 12 11:08:57
                             Started mapping on |       Sep 12 11:11:27
                                    Finished on |       Sep 12 11:11:29
       Mapping speed, Million of reads per hour |       18.00

                          Number of input reads |       10000
                      Average input read length |       73
                                    UNIQUE READS:
                   Uniquely mapped reads number |       47
                        Uniquely mapped reads % |       0.47%
                          Average mapped length |       51.66
                       Number of splices: Total |       5
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       3
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       2
                      Mismatch rate per base, % |       2.14%
                         Deletion rate per base |       0.04%
                        Deletion average length |       1.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       83
             % of reads mapped to multiple loci |       0.83%
        Number of reads mapped to too many loci |       19
             % of reads mapped to too many loci |       0.19%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.02%
                 % of reads unmapped: too short |       98.31%
                     % of reads unmapped: other |       0.18%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%

    total unique multiple map unmap
    """
    total = 0
    mapped = 0
    unmapped = 0
    multi = 0
    unique = 0
    out = 0
    is_paired = False
    warn_chunkmbs = False
    if check_file(x, check_empty=True):
        with open(x) as r:
            for line in r:
                p = line.strip().split()
                if len(p) < 2:
                    continue
                n = p[-1]
                if not re.match('^[0-9]+$', n):
                    continue
                n = eval(n)
                if 'Number of input reads' in line:
                    total = n
                elif 'Uniquely mapped reads number' in line:
                    unique = n
                elif 'Number of reads mapped to multiple loci' in line:
                    multi = n
                elif 'Number of reads mapped to too many loci' in line:
                    multi2 = n
        # if unique_only
        # multi = 0, multi2 >0
        if multi == 0 and multi2 > 0:
            multi = multi2
    else:
        log.error('file not exists, {}'.format(x))
    # msg
    if warn_chunkmbs:
        log.warning('{}\nset --chunkmbs 128, to fix the errors'.format(
            warn_chunkmbs))
    return {
        'total': total,
        'map': unique + multi,
        'unique': unique,
        'multi': multi,
        'unmap': total - unique - multi,
        }


class StarConfig(object):
    """Check args, prepare files for STAR
    arguments
    output files
    parser ?!
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'STAR',
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'index': None,            
            'index_name': None,
            'smp_name': None,
            'extra_para': None,
            'smp_name': None,
            'threads': 1,
            'overwrite': False,
            'n_map': 20,
            'unique_only': False,
            'keep_tmp': False,
            'keep_unmap': True,
            'large_insert': False,
            'genomeLoad': 'LoadAndRemove',
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'STAR_r1'
        self.init_fx()
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        # index name
        if not AlignIndex(self.index, self.aligner).is_valid():
            raise ValueError('index not valid, {}'.format(self.index))
        if self.index_name is None:
            self.index_name = AlignIndex(self.index).index_name()
        # samll genome
        self.small_genome = AlignIndex(
            self.index, self.aligner).index_size() < 10000000
        self.seed_max = 5 if self.small_genome else 50
        self.n_map = 1 if self.unique_only else self.n_map
        if self.n_map < 1:
            self.n_map = 20 # default, see unique_only
        if self.smp_name is None:
            self.smp_name = fx_name(self.fq1, fix_pe=True)
        # update files
        self.init_files()


    def init_fx(self):
        """Make sure, fx
        1. str
        2. fq1 exists
        3. fq2 None or exists
        """
        if not check_fx(self.fq1):
            raise ValueError('--fq1, not exists, or empty')
        if self.fq2 is not None and not check_fx(self.fq2):
            raise ValueError('--fq2, not exists, or empty')
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('--fq1, --fq2 faild, not properly paired')
        # format
        self.fx_format = Fastx(self.fq1).format # fasta/q
        self.is_paired = file_exists(self.fq2) # fq2


    def init_files(self):
        self.project_dir = os.path.join(self.outdir, self.smp_name, 
            self.index_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        # output files
        prefix = os.path.join(self.project_dir, self.smp_name)
        default_files = {
#             'project_dir': self.project_dir,
            'config_toml': os.path.join(self.config_dir, 'config.toml'),
            'cmd_shell': os.path.join(self.project_dir, 'cmd.sh'),
            'bam': prefix + '.bam',
            'sam': prefix + '.sam',
            'unmap': prefix + '.unmap.' + self.fx_format,
            'unmap1': prefix + '.unmap.1.' + self.fx_format, # 
            'unmap2': prefix + '.unmap.2.' + self.fx_format, #
            'align_log': prefix + '.align.log',
            'align_stat': prefix + '.align.stat',
            'align_json': prefix + '.align.json',
            'align_flagstat': prefix + '.flagstat',
            'bam_raw': prefix + 'Aligned.sortedByCoord.out.bam',
            'log_raw': prefix + 'Log.final.out',
            'unmap_raw': prefix + 'Unmapped.out',
            'unmap1_raw': prefix + 'Unmapped.out.mate1',
            'unmap2_raw': prefix + 'Unmapped.out.mate2',
        }
        self = update_obj(self, default_files, force=True)
        self.align_prefix = prefix
        check_path([self.project_dir, self.config_dir], create_dirs=True)


class Star(object):
    """Alignment, using STAR
    Single index, SE/PE
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = StarConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'STAR' # force changed
        self.get_cmd()        
        Config().dump(self.__dict__, self.config_toml)


    def get_cmd(self):
        """The command line
        ########################
        # STAR unique mapping  #
        ########################
        --outFilterMultimapNmax 1
        
        --seedPerWindowNmax 5, for small genome

        Parameters
        ----------
        For sharing memory in STAR
        by Devon Ryan: https://www.biostars.org/p/260069/#260077 
        by Dobin: https://github.com/alexdobin/STAR/pull/26
        
        --genomeLoad

        NoSharedMemory: each job use its own copy
        LoadAndExit:   load genome to memory, do not run alignment
        Remove:        remove genome from memory, do not run alignment
        LoadAndRemove: load genome to memory, remove genome after alignment
        LoadAndKeep:   load genome to memory, keep it after run

        pratice: for general usage
        1. NoSharedMemory: (?) what if other jobs using it?

        pratice: for general usage
        1. LoadAndRemove: (?) 

        pratice: for multiple alignment together.
        1. LoadAndExit
        2. (loop over samples): LoadAndKeep
        3. Remove
        """
        args_extra = self.extra_para if self.extra_para else ''
        args_reader = 'zcat' if self.fq1.endswith('.gz') else '-'
        args_fq2 = self.fq2 if self.is_paired else ''
        # command
        cmd_main = ' '.join([
            '{}'.format(shutil.which('STAR')),
            '--genomeLoad {}'.format(self.genomeLoad), #
            '--runMode alignReads',
            '--genomeDir {}'.format(self.index),
            '--readFilesCommand {}'.format(args_reader),            
            '--readFilesIn {} {}'.format(self.fq1, args_fq2),
            '--outFileNamePrefix {}'.format(self.align_prefix),
            '--runThreadN {}'.format(self.threads),
            '--limitBAMsortRAM 10000000000',
            '--outSAMtype BAM SortedByCoordinate',
            '--outFilterMismatchNoverLmax 0.07',
            '--seedSearchStartLmax 20',
            '--outReadsUnmapped Fastx', # self.unmap1,
            '--outFilterMultimapNmax {}'.format(self.n_map),
            '--seedPerWindowNmax {}'.format(self.seed_max),
            args_extra,
            ])
        # fix file names
        self.cmd = cmd_main


    def fix_out_files(self):
        """Update the filenames of STAR output
        bam: *Aligned.sortedByCoord.out.bam -> *.bam # mv
        log: *Log.final.out -> *.log # copy
        log: *Log.out -> *.out
        unmap: *Unmapped.out.mate1 -> *.unmap.1.fastq
               *Unmapped.out.mate1 -> *.unmap.1.fastq
        """
        ## default output of STAR
        ## bam, log, unmap files
        # new files
        symlink_file(self.bam_raw, self.bam)
        symlink_file(self.log_raw, self.align_log)
        if self.is_paired:
            symlink_file(self.unmap1_raw, self.unmap1)
            symlink_file(self.unmap2_raw, self.unmap2)
        else:
            symlink_file(self.unmap1_raw, self.unmap1)
            self.unmap2 = None
        # remove old files
        del_list = []
        if not self.keep_unmap:
            del_list.extend([
                self.unmap1, self.unmap2, self.unmap,
                self.unmap1_raw, self.unmap2_raw])
        remove_file(del_list, ask=False)


    def run(self):
        if file_exists(self.bam) and not self.overwrite:
            log.info('Star() skipped, file exists: {}'.format(self.bam))
        else:
            # save cmd
            with open(self.cmd_shell, 'wt') as w:
                w.write(self.cmd + '\n')
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('Star() failed, check {}'.format(self.align_log))
        # for output files
        self.fix_out_files()
        if not check_file(self.bam, check_empty=True):
            log.error('Star() failed, not bam: {}'.format(self.bam))
        # log
        df = parse_star(self.align_log)
        df.update({
            'name': self.smp_name,
            'index': self.index_name,
            'unique_only': self.unique_only,
            })
        Config().dump(df, self.align_json)
        return (self.bam, self.unmap1, self.unmap2)


def get_args():
    """Parsing arguments for STAR
    """
    example = '\n'.join([
        'Examples:',
        '$ python star.py -1 f1.fq -x genome -o output',
        '# add extra para',
        '$ python star.py -1 f1.fq -2 f2.fq -x genome -o output',
        '# unique reads, update index_name',
        '$ python star.py -1 f1.fq -x genome -o output -u -in 01.genome',
    ])    
    parser = argparse.ArgumentParser(
        prog='run_star',
        description='run STAR program',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
                        help='Fasta/q file, read1 of PE, or SE read')
    parser.add_argument('-2', '--fq2', required=False, default=None,
                        help='Fasta/q file, read2 of PE, or SE read, optional')
    parser.add_argument('-x', '--index', required=True,
                        help='The alignment index for STAR')
    parser.add_argument('-o', '--outdir', default=None,
                        help='Directory saving results, default: [cwd]')
    parser.add_argument('-in', '--index-name', default=None, dest='index_name',
                        help='The name of the index')
    parser.add_argument('-n', '--smp-name', default=None, dest='smp_name',
                        help='The name of the sample')
    parser.add_argument('-p', '--threads', default=1, 
                        help='Number of threads, default: [1]')
    parser.add_argument('-w', '--overwrite', action='store_true',
                        help='Overwrite the exist files')
    parser.add_argument('-u', '--unique-only', action='store_true',
                        dest='unique_only', 
                        help='Report unique mapped reads only')
    parser.add_argument('-l', '--largs-insert', action='store_true',
                        dest='large_insert',
                        help='For large insert, use: -X 1000 --chunkmbs 128')
    parser.add_argument('--clean', dest='keep_tmp', action='store_false',
                        help='Clean temp files')
    parser.add_argument('-X', '--extra-para', dest='extra_para', default=None,
                        help='Add extra parameters, eg: "-X 2000"')
    return parser


def main():
    args = vars(get_args().parse_args())
    # update: keep_tmp, keep_unmap
    args['keep_unmap'] = args['keep_tmp']
    Star(**args).run()


if __name__ == '__main__':
    main()
    
    