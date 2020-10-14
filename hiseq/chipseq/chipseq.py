"""
ChIP-seq pipeline

Working mode:

ChIPseqSingle(): 
1. IP / Input
  - raw_data
  - clean_data
  - align
  - bam_files
  - bw_files
  - peak
  - qc (multiqc, fastqc)

2. bam_files
3. bw_files (IP - Input)
4. peak_files (IP - Input)
5. report

# for replicates
# merge

Output:

- raw_data 
- clean_data
- align  
- bam_files 
- bw_files
- peak 
- motif 
- report
  - reads num
  - peaks num 
  - peaks annotation 
  - qc 

"""


import os
import re
from multiprocessing import Pool
from Levenshtein import distance
from hiseq.utils.helper import *
from hiseq.trim.trimmer import Trimmer
from hiseq.align.alignment import Alignment
from hiseq.peak.call_peak import Macs2
from hiseq.atac.atac_utils import *
from hiseq.fragsize.fragsize import BamPEFragSize
import hiseq
import collections


def print_dict(d):
    d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


def check_fq(fq):
    """
    Make sure
    fq: str or list, or None
    """
    if fq is None:
        # raise ValueError('fq1 required, got None')
        pass
    elif isinstance(fq, list):
        fq = file_abspath(fq)
    elif isinstance(fq, str):
        fq = [file_abspath(fq)]
    else:
        log.error('fq failed, Nont, str, list expected, got {}'.format(type(fq).__name__))

    return fq


def fq_paired(fq1, fq2):
    """
    Make sure fq1 and fq2, proper paired
    """
    fq1 = check_fq(fq1)
    fq2 = check_fq(fq2)

    if isinstance(fq1, str) and isinstance(fq2, str):
        return distance(fq1, fq2) == 1
    elif isinstance(fq1, list) and isinstance(fq2, list):
        return [distance(i, j) == 1 for i, j in zip(fq1, fq2)]
    else:
        log.warning('fq not paired: {}, {}'.format(fq1, fq2))
        return False


def init_cpu(threads=1, parallel_jobs=1):
    """
    threads, CPUs
    """
    n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()

    max_jobs = int(n_cpu / 4.0)
    ## check parallel_jobs (max: 1/4 of n_cpus)
    if parallel_jobs > max_jobs: 
        log.warning('Too large, change parallel_jobs from {} to {}'.format(
            parallel_jobs, max_jobs))
        parallel_jobs = max_jobs

    ## check threads
    max_threads = int(0.8 * n_cpu / parallel_jobs)
    if threads * parallel_jobs > 0.8 * n_cpu:
        log.warning('Too large, change threads from {} to {}'.format(
            threads, max_threads))
        threads = max_threads

    return (threads, parallel_jobs)


def update_obj(obj, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in obj.__dict__:
                delattr(obj, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(obj, k) or force:
                    setattr(obj, k, v)

        return obj


class ChIPseqR1Config(object):
    """
    Prepare directories for R1 single replicates

    require: fq1/fq2, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()


    def init_args(self):
        """
        required arguments for ChIPseq analysis
        """
        args_init = {
            'is_ip': True,
            'is_trimmed': False,
            'smp_name': None,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': str(pathlib.Path.cwd()),
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0
        }
        self = update_obj(self, args_init, force=False)
        self.chipseq_type = 'chipseq_r1' # 

        # check fastq files
        if self.fq1 is None:
            raise ValueError('ChIPseqR1Config failed, fq1 should not be NoneType')
        elif isinstance(self.fq1, str):
            # self.fq1, self.fq2 = self.check_fq_paired(self.fq1, self.fq2)
            if not self.fq2 is None:
                if not fq_paired(self.fq1, self.fq2):
                    raise ValueError('fq1, fq2 not paired properly: {}, {}'.format(self.fq1, self.fq2))
        else:
            raise ValueError('ChIPseqR1Config failed, fq1 expect str, got {}'.format(type(self.fq1).__name__))

        # smp_name
        self.smp_name = getattr(self, 'smp_name', None)
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, self.parallel_jobs)

        self.outdir = file_abspath(self.outdir)


    def init_dirs(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')

        auto_files = {
            ## dirs
            'raw_dir': os.path.join(self.project_dir, 'raw_data'),
            'clean_dir': os.path.join(self.project_dir, 'clean_data'),
            'align_dir': os.path.join(self.project_dir, 'align'),
            'bam_dir': os.path.join(self.project_dir, 'bam_files'),
            'bw_dir': os.path.join(self.project_dir, 'bw_files'),
            'peak_dir': os.path.join(self.project_dir, 'peak'),
            'motif_dir': os.path.join(self.project_dir, 'motif'),
            'qc_dir': os.path.join(self.project_dir, 'qc'),
            'qc_lendist_dir': os.path.join(self.project_dir, 'qc', 'lendist'),
            'qc_frip_dir': os.path.join(self.project_dir, 'qc', 'FRiP'),
            'report_dir': os.path.join(self.project_dir, 'report'),
            ## files
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'align_stat': os.path.join(self.project_dir, 'align', self.smp_name + '.align.txt'),
            # 'bam_raw': from align_dir, to-do
            'bam_rmdup': os.path.join(self.project_dir, 'bam_files', self.smp_name + '.rmdup.bam'),
            'bam_proper_pair': os.path.join(self.project_dir, 'bam_files', self.smp_name + '.proper_pair.bam'),
            'bam': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bam'),
            'bed': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bed'),
            'peak': os.path.join(self.project_dir, 'peak', self.project_name + '_peaks.narrowPeak'),
            'bw': os.path.join(self.project_dir, 'bw_files', self.project_name + '.bigWig'),
            # qc
            'lendist_txt': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.txt'),
            'lendist_pdf': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.pdf'),
            'frip_txt': os.path.join(self.project_dir, 'qc', 'FRiP', 'FRiP.txt'),
            }
        self = update_obj(self, auto_files, force=True) # key

        ## raw data
        self.raw_fq_list = [os.path.join(self.raw_dir, os.path.basename(self.fq1))]
        fq2_raw = None if self.fq2 is None else os.path.join(self.raw_dir, os.path.basename(self.fq2))
        self.raw_fq_list.append(fq2_raw)

        ## clean data
        self.clean_fq_list = [os.path.join(self.clean_dir, fq_name(self.fq1, pe_fix=False) + '.fq.gz')]
        fq2_clean = None if self.fq2 is None else os.path.join(self.clean_dir, fq_name(self.fq2, pe_fix=False) + '.fq.gz')
        self.clean_fq_list.append(fq2_clean)

        if create_dirs:
            check_path([
                self.project_dir, 
                self.config_dir, 
                self.raw_dir,
                self.clean_dir,
                self.align_dir,
                self.bam_dir, 
                self.bw_dir, 
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir,
                os.path.join(self.qc_dir, 'lendist'),
                os.path.join(self.qc_dir, 'FRiP') ])


class ChIPseqRnConfig(object):
    """
    Prepare directories/args for multiple replicates

    require: rep_list/fq1-list, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()


    def init_args(self):
        """
        required arguments for ChIPseq analysis
        """
        args_init = {
            'is_ip': True,
            'smp_name': None,
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': str(pathlib.Path.cwd()),
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0
        }
        self = update_obj(self, args_init, force=False)
        self.chipseq_type = 'chipseq_rn'

        # check fastq files/build rep_list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]

        if isinstance(self.fq1, list):
            self.rep_list = [os.path.join(self.outdir, fq_name(i, pe_fix=True)) for i in self.fq1]

        if isinstance(self.rep_list, list):
            self.smp_name = fq_name_rmrep(self.rep_list)
            self.smp_name = self.smp_name.pop() # to str !!! to-do, check unique
        else:
            raise ValueError('ChIPseqRnConfig failed, fq1/rep_list expect list, got {}'.format(type(self.rep_list).__name__))

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, self.parallel_jobs)

        self.outdir = file_abspath(self.outdir)


    def init_dirs(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        self.rep_list = file_abspath(self.rep_list)

        # path, files
        self.project_name = self.smp_name #
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        auto_files = {
            ## dirs
            'align_dir': os.path.join(self.project_dir, 'align'),
            'bam_dir': os.path.join(self.project_dir, 'bam_files'),
            'bw_dir': os.path.join(self.project_dir, 'bw_files'),
            'peak_dir': os.path.join(self.project_dir, 'peak'),
            'motif_dir': os.path.join(self.project_dir, 'motif'),
            'qc_dir': os.path.join(self.project_dir, 'qc'),
            'qc_lendist_dir': os.path.join(self.project_dir, 'qc', 'lendist'),
            'qc_frip_dir': os.path.join(self.project_dir, 'qc', 'FRiP'),
            'qc_idr_dir': os.path.join(self.project_dir, 'qc', 'IDR'),
            'qc_cor_dir': os.path.join(self.project_dir, 'qc', 'cor'),
            'qc_overlap_dir': os.path.join(self.project_dir, 'qc', 'overlap'),
            'report_dir': os.path.join(self.project_dir, 'report'),
            ## files
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'bam': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bam'),
            'peak': os.path.join(self.project_dir, 'peak', self.project_name + '_peaks.narrowPeak'),
            'bw': os.path.join(self.project_dir, 'bw_files', self.project_name + '.bigWig'),
            ## qc files
            'lendist_txt': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.txt'),
            'lendist_pdf': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.pdf'),
            'frip_txt': os.path.join(self.project_dir, 'qc', 'FRiP', 'FRiP.txt'),
            'cor_npz': os.path.join(self.project_dir, 'qc', 'cor', 'cor.bam.npz'),
            'cor_counts': os.path.join(self.project_dir, 'qc', 'cor', 'cor.bam.counts.tab'),
            'idr_txt': os.path.join(self.project_dir, 'qc', 'IDR', 'dir.txt'),
            'idr_log': os.path.join(self.project_dir, 'qc', 'IDR', 'dir.log'),
            'peak_overlap_pdf': os.path.join(self.project_dir, 'qc', 'overlap', 'peak_overlap_pdf')        
            }
        self = update_obj(self, auto_files, force=True) # key

        if create_dirs:
            check_path([
                self.project_dir, 
                self.config_dir, 
                self.align_dir,
                self.bam_dir, 
                self.bw_dir, 
                self.peak_dir, 
                self.qc_dir, 
                self.report_dir,
                os.path.join(self.qc_dir, 'lendist'),
                os.path.join(self.qc_dir, 'FRiP'),
                os.path.join(self.qc_dir, 'cor'),
                os.path.join(self.qc_dir, 'IDR'),
                os.path.join(self.qc_dir, 'overlap')])


class ChIPseqRxConfig(object):
    """
    Prepare directories/args for IP/Input

    require: ip_rep_list, input_rep_list, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()


    def init_args(self):
        """
        required arguments for ChIPseq analysis
        """
        args_init = {
            'chipseq_type': 'chipseq_rx',
            'ip_dir': None,
            'input_dir': None,
            'genome': None,
            'outdir': str(pathlib.Path.cwd()),
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0
        }
        self = update_obj(self, args_init, force=False)
        self.chipseq_type = 'chipseq_rx'

        # check ip/input
        if self.ip_dir is None or self.input_dir is None:
            raise ValueError('ip_dir, input_dir required, ip:{}, \
                input:{}'.format(ip_dir, input_dir))

        # check smp_name
        if self.smp_name is None:
            self.ip_name = ChIPseqReader(self.ip_dir).args.get('smp_name', None)
            self.input_name = ChIPseqReader(self.input_dir).args.get('smp_name', None)
            self.smp_name = '{}.vs.{}'.format(self.ip_name, self.input_name)

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, 
            self.parallel_jobs)

        self.outdir = file_abspath(self.outdir)


    def init_dirs(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        # path, files
        self.project_name = self.smp_name #
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.bam_dir = os.path.join(self.project_dir, 'bam_files')
        self.bw_dir = os.path.join(self.project_dir, 'bw_files')
        self.peak_dir = os.path.join(self.project_dir, 'peak')
        self.motif_dir = os.path.join(self.project_dir, 'motif')
        self.qc_dir = os.path.join(self.project_dir, 'qc')
        self.report_dir = os.path.join(self.project_dir, 'report')
        auto_files = {
            ## config files
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            ## output files
            'ip_bam': os.path.join(self.bam_dir, self.ip_name + '.bam'),
            'ip_bw':  os.path.join(self.bw_dir, self.ip_name + '.bigWig'),
            'input_bam': os.path.join(self.bam_dir, self.input_name + '.bam'),
            'input_bw': os.path.join(self.bw_dir, self.input_name + '.bigWig'),
            'peak': os.path.join(self.peak_dir, self.ip_name + '_peaks.narrowPeak'),
            'ip_over_input_bw': os.path.join(self.bw_dir, self.ip_name + '.ip_over_input.bigWig'),  
            'bw': os.path.join(self.project_dir, 'bw_files', self.ip_name + '.bigWig')}
        self = update_obj(self, auto_files, force=True) # key

        if create_dirs:
            check_path([
                self.project_dir, 
                self.config_dir, 
                self.bam_dir, 
                self.bw_dir, 
                self.peak_dir, 
                self.qc_dir, 
                self.report_dir])


class ChIPseqConfig(object):
    """
    Global config for chipseq analysis
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()
        self.init_mission()


    def init_args(self):
        """
        required arguments for ChIPseq analysis
        """
        args_init = {
            'build_design': False,
            'design': None,
            'is_ip': True,
            'ip': None,
            'input': None,
            'fq1': None,
            'fq2': None,
            'rep_list': None,
            'ip_dir': None,
            'input_dir': None,
            'smp_name': None,
            'genome': None,
            'outdir': str(pathlib.Path.cwd()),
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0
        }
        self = update_obj(self, args_init, force=False)

        # # update fq from design
        # if file_exists(self.design) and not self.build_design:
        #     self.design = file_abspath(self.design)
        #     args_design = Json(self.design).reader()
        #     self = update_obj(self, args_design, force=True)

        # print('!AAAA-1')
        # print_dict(self.__dict__)
        # sys.exit()
        if self.outdir is None:
          self.outdir = str(pathlib.Path.cwd())
        
        # 1st level
        if not self.fq1 is None:
            self.fq1, self.fq2 = self.check_fq_paired(self.fq1, self.fq2)

        # 2nd level
        if isinstance(self.rep_list, list):
            self.rep_list = file_abspath(self.rep_list)
        elif isinstance(self.rep_list, str):
            self.rep_list = [file_abspath(self.rep_list)]
        else:
            pass
            # raise ValueError('rep_list failed, None, str, list expected, got {}'.format(type(self.rep_list).__name__))


    def check_fq(self, fq):
        """
        Make sure
        fq: str or list, or None
        """
        if fq is None:
            # raise ValueError('fq1 required, got None')
            pass
        elif isinstance(fq, list):
            fq = file_abspath(fq)
        elif isinstance(fq, str):
            fq = [file_abspath(fq)]
        else:
            log.error('fq failed, Nont, str, list expected, got {}'.format(type(fq).__name__))

        return fq


    def check_fq_paired(self, fq1, fq2=None):
        """
        Make sure fq1 and fq2, proper paired
        """
        fq1 = self.check_fq(fq1)
        fq2 = self.check_fq(fq2)

        # fq1 = fq2
        if isinstance(fq2, list):
            if not len(fq1) == len(fq2):
                raise ValueError('fq1, fq2, files not paired correctly.')

            # check
            q_tag = 0
            for q1, q2 in zip(fq1, fq2):
                if not distance(fq1, fq2) == 1:
                    log.error('PE fq files failed, {}, {}'.format(q1, q2))
                    q_tag += 1

            if q_tag:
                raise ValueError('fastq files not paired correctly.')

        return (fq1, fq2)


    def init_dirs(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        self.project_dir = os.path.join(self.outdir)
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.config_txt = os.path.join(self.config_dir, 'arguments.txt')
        self.config_pickle = os.path.join(self.config_dir, 'arguments.pickle')
        self.config_json = os.path.join(self.config_dir, 'arguments.json')

        if create_dirs:
            check_path(self.config_dir)


    def init_mission(self):
        """
        Determine the type of ChIPseq analysis
        1. build_design
        2. Single: IP/Input pair
        3. Multiple: IP/Input pair
        """
        create_dirs = getattr(self, 'create_dirs', True)

        self.chipseq_type = None

        # 1st level
        if self.build_design:
            self.chipseq_type = 'build_design'
            # self.init_build_design()
        elif file_exists(self.design):
            self.chipseq_type = 'chipseq_rx_from_design'
        elif isinstance(self.ip_dir, str) and isinstance(self.input_dir, str):
            self.chipseq_tyep = 'chipseq_rx'
        # elif isinstance(self.ip, list) and isinstance(self.input, list):
        elif all([not i is None for i in [self.ip, self.input]]):
            self.chipseq_type = 'chipseq_rx'
        elif isinstance(self.rep_list, list) or isinstance(self.fq1, list):
            self.chipseq_type = 'chipseq_rn'
        elif isinstance(self.fq1, str):
            self.chipseq_type = 'chipseq_r1'
        else:
            raise ValueError('unknown chipseq type')

        # check
        if self.chipseq_type is None:
            raise ValueError('unknown chipseq_type')


class ChIPseqR1(object):
    """
    Run ChIPseq basic for single fastq file (ip/input)

    align - bam - rmdup - proper_paired - peak - bw - qc - report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        obj_local = ChIPseqR1Config(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


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
        hiseq.trim.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        args_trim = self.__dict__.copy()
        # args = self.args.copy()
        fq1, fq2 = self.raw_fq_list
        clean_fq1, clean_fq2 = self.clean_fq_list

        # update args
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.clean_dir,
            'library_type': 'Nextera'
        }
        args_trim.update(args_init)
        
        if trimmed is True:
            ## fq1
            if is_gz(fq1):
                symlink(fq1, clean_fq1, absolute_path=True)
            else:
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
            ## fq2
            if not fq2 is None:
                if is_gz(fq2):
                    symlink(fq2, clean_fq2, absolute_path=True)
                else:
                    gzip_cmd(fq2, clean_fq2, decompress=False, rm=False)
        else:
            if check_file(self.clean_fq_list[0]):
                log.info('trim() skipped, file exists: {}'.format(
                    self.clean_fq_list))
            else:
                trimmer = Trimmer(**args_trim)
                trimmer.run()
                fq1, fq2 = trimmer.out_files
                symlink(fq1, clean_fq1)
                symlink(fq2, clean_fq2)


    def align(self):
        """
        Alignment PE reads to reference genome, using bowtie2

        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'fq': fq1,
            'outdir': self.align_dir,
            'extra_para': '--very-sensitive-local', # --no-mixed -X 800
            'aligner': 'bowtie2'
        }
        args_local.update(args_init)

        # output
        if check_file(args_local.get('bam_rmdup', None)):
            log.info('align() skipped, file exists: {}'.format(
                args_local.get('bam_rmdup', None)))
        else:
            Alignment(**args_local).run()


    def get_raw_bam(self):
        """
        Get the align bam file
        from bam_dir/1., 2., ...
        !!! specific: 2.genome/*.bam
        """
        bamlist = listfile(self.align_dir, '*.bam', recursive=True)
        bamlist = sorted(bamlist)
        bamlist = [i for i in bamlist if not i.endswith('.raw.bam')] # remove raw.bam
        # [chrM, genome]
        # 
        return(bamlist[-1]) # the last one


    def get_bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        bam_raw = self.get_raw_bam()

        # rmdup
        if rmdup is True:
            Bam(bam_raw).rmdup(self.bam_rmdup)
        else:
            shutil.copy(bam_raw, self.bam_rmdup)

        # bam_proper_pair = bam
        if check_file(self.bam):
            log.info('bam_rmdup() skipped, file exists: {}'.format(
                self.bam))
        else:
            # save proper pair reads
            # Bam(self.bam_rmdup).proper_pair(self.bam)
            Bam(self.bam_rmdup).sort(self.bam)
            # index
            Bam(self.bam).index()


    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 0,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }

        Bam2bw(**args_local).run()


    def call_peak(self):
        """
        Call peaks using MACS2
        """
        args_peak = self.__dict__.copy()
        args_peak['genome_size'] = getattr(self, 'genome_size', 0)
        bam = self.bam
        bed = os.path.splitext(bam)[0] + '.bed'
        Bam(bam).to_bed(bed)
        args_peak.pop('genome', None)
        args_peak.pop('outdir', None)
        if check_file(self.peak):
            log.info('call_peak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, self.genome, self.peak_dir, self.project_name, atac=False).callpeak()


    def qc_lendist(self):
        """
        Create length distribution, txt
        """
        if os.path.exists(self.lendist_txt) and self.overwrite is False:
            log.info('lendist() skipped: file exists: {}'.format(self.lendist_txt))
        else:
            # _ = frag_length(self.bam, self.lendist_txt)
            BamPEFragSize(self.bam).saveas(self.lendist_txt)


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            frip, n, total = peak_FRiP(self.peak, self.bam, threads=1) # threads > 1, not allowed,

            hd = ['FRiP', "peak_reads", "total_reads", "id"]
            n = list(map(str, [frip, n, total]))
            # n.append('self.config.fqname')
            n.append(self.project_name)
            with open(self.frip_txt, 'wt') as w:
                w.write('\t'.join(hd) + '\n')
                w.write('\t'.join(n) + '\n')


    def qc_mito(self):
        """
        Mito reads, percentage
        """
        pass


    def qc_tss(self):
        """
        TSS enrichment
        require:
        BAM, peak, TSS, ...
        """
        pass


    def qc(self):
        """
        QC for ChIPseq sample
        FRiP
        Fragment size (length distribution)
        TSS
        """
        self.qc_frip()
        self.qc_lendist()
        self.qc_mito()
        self.qc_tss()


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'chipseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'ChIPseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        run_shell_cmd(cmd)
        # try:
        #     run_shell_cmd(cmd)
        # except:
        #     log.warning('report() failed.')


    def run(self):
        """
        Run all steps for ChIPseq pipeline
        """
        # init dir
        args = self.__dict__.copy()

        trimmed = args.get('trimmed', False)
        copy_raw_fq = args.get('copy_raw_fq', False)
        rmdup = args.get('rmdup', True) # remove PCR dup

        # 1. copy raw data
        self.prep_raw(copy_raw_fq)

        # 2. trim
        self.trim(trimmed)

        # 3. align
        self.align()

        # 4. rmdup
        self.get_bam_rmdup(rmdup)

        # 5. bw
        self.bam_to_bw()

        # 5. peak 
        self.call_peak()
        
        # 6. motif
        # self.motif()

        # 7. qc
        self.qc()

        # 8.report
        self.report()


class ChIPseqRn(object):
    """
    Run ChIPseq for multiple replicates
    input: rep_list/fq1
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = ChIPseqRnConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)
        # self.chipseq_type = 'chipseq_rn'

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def get_bam_list(self):
        """
        get the proper_mapped.bam files
        """
        return [ChIPseqReader(i).args.get('bam', None) for i in self.rep_list]


    def get_peak_list(self):
        """
        get the .narrowPeak files
        """
        return [ChIPseqReader(i).args.get('peak', None) for i in self.rep_list]


    def merge_bam(self):
        """
        Merge replicates, BAM
        """
        self.bam_list = self.get_bam_list()

        cmd = ' '.join([
            'samtools merge {}'.format(self.bam + '.tmp'),
            ' '.join(self.bam_list),
            '&& samtools sort -o {} {}'.format(self.bam, self.bam + '.tmp'),
            '&& samtools index {}'.format(self.bam)])

        if os.path.exists(self.bam):
            log.info('merge_bam() skipped, file exists: {}'.format(
                self.bam))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('merge_bam() failed.')


    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 0,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }

        Bam2bw(**args_local).run()


    def call_peak(self):
        """
        Call peaks using MACS2
        ...
        """
        args_peak = self.__dict__.copy()
        args_peak['genome_size'] = getattr(self, 'genome_size', 0)

        
        bam = self.bam
        bed = os.path.splitext(bam)[0] + '.bed'
        Bam(bam).to_bed(bed)
        genome = args_peak.pop('genome', None)
        output = args_peak.pop('peak_dir', None)
        prefix = args_peak.pop('smp_name', None)

        if check_file(self.peak):
            log.info('call_peak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, genome, output, prefix, atac=False).callpeak()


    def get_bam_cor(self, window=500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp
        
        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
            --outRawCounts *counts.tab -b bam
        """
        args_local = {
            'bam': self.get_bam_list(),
            'outdir': self.qc_cor_dir,
            'threads': self.threads,
            'overwrite': self.overwrite,
            'binsize': self.binsize
        }
        Bam2cor(**args_local).run()


    def get_peak_overlap(self):
        """
        Compute the overlaps between overlaps
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_overlap_dir,
            'overwrite': self.overwrite
        }   

        BedOverlap(**args_local).run()


    def get_peak_idr(self):
        """
        Calculate IDR for replicates
        1 vs 1
        peak files
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_idr_dir,
            'overwrite': self.overwrite
        }

        PeakIDR(**args_local).run()


    def get_peak_frip(self):
        """
        Save all FRiP.txt file to one
        """
        # get list
        self.frip_list = []
        for fq in self.fq1:
            fq_dir = os.path.join(self.outdir, fq_name(fq, pe_fix=True))
            self.frip_list.append(
                ChIPseqReader(fq_dir).args.get('frip_txt', None))

        if not os.path.exists(self.frip_txt):
            with open(self.frip_txt, 'wt') as w:
                for f in self.frip_list:
                    with open(f) as r:
                        for line in r:
                            w.write(line)


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            frip, n, total = peak_FRiP(self.peak, self.bam, threads=1) # threads > 1, not allowed,

            hd = ['FRiP', "peak_reads", "total_reads", "id"]
            n = list(map(str, [frip, n, total]))
            # n.append('self.config.fqname')
            n.append(self.project_name)
            with open(self.frip_txt, 'wt') as w:
                w.write('\t'.join(hd) + '\n')
                w.write('\t'.join(n) + '\n')


    def get_align_txt(self):
        """
        Copy align.txt files to align/
        """
        pass


    def get_tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        pass


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'chipseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'chipseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        print('!AAAA-4', 'report')
        run_shell_cmd(cmd) 


    def pick_fq_samples(self, i):
        """
        Pick samples for each group
        """
        if isinstance(i, int):
            # fq1
            if isinstance(self.fq1, list):
                fq1 = self.fq1[i]
            else:
                fq1 = None

            # fq2
            if isinstance(self.fq2, list):
                fq2 = self.fq2[i]
            else:
                fq2 = None

            # rep_list
            rep_list = None

            args_fq = {
                'fq1': fq1,
                'fq2': fq2,
                'rep_list': rep_list}

            return(args_fq)
        else:
            raise ValueError('int expected, {} found'.format(type(g).__name__))


    def run_fq_single(self, i):
        """
        Run for single group
        """
        args_tmp = self.__dict__.copy()
        # required args
        args_required = ['align_to_chrM', 'aligner', 'fq1', 'fq2', 'genome', 
            'genome_size', 'is_ip', 'is_trimmed', 'outdir', 'overwrite', 
            'parallel_jobs', 'threads']
        args_local = dict((k, args_tmp[k]) for k in args_required 
            if k in args_tmp)

        args_input = self.pick_fq_samples(i)
        args_init = {
            'build_design': None,
            'design': None
        }
        args_local.update(args_input) # update fq1/rep_list/group
        args_local.update(args_init) #

        obj_local = ChIPseqR1Config(**args_local)
        ChIPseqR1(**obj_local.__dict__).run()


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # run each fq
        if isinstance(self.fq1, list):
            i_list = list(range(len(self.fq1)))
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_fq_single, i_list)

            # # # alternative
            # for i in i_list:
            #     self.run_fq_single(i)

        if len(self.rep_list) > 1:
            print('!AAA', self.rep_list)
            # run
            self.merge_bam()
            self.bam_to_bw()
            self.call_peak()
            # qc
            self.get_bam_cor()
            self.get_peak_overlap()
            self.get_peak_idr()
            self.get_peak_frip()
            self.qc_frip()
            self.get_align_txt()
            self.get_tss_enrich()
            self.report()
        elif len(self.rep_list) == 1:
            log.warning('merge() skipped, Only 1 replicate detected')
            # copy files: bam, bw, peak
            rep_dict = ChIPseqReader(self.rep_list[0]).args
            symlink(rep_dict.get('bam', None), self.bam)
            symlink(rep_dict.get('bw', None), self.bw)
            symlink(rep_dict.get('peak', None), self.peak)
        else:
            log.error('merge() failed, no rep detected')
            raise ValueError('merge() failed, no rep detected')


class ChIPseqRx(object):
    """
    Run ChIPseq for ip and input
    input: ip_dir, input_dir
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = ChIPseqRxConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)

        ## check args
        self.ip_args = ChIPseqReader(self.ip_dir).args
        self.input_args = ChIPseqReader(self.input_dir).args


    def get_bam_files(self):
        """
        Copy bam files
        """
        symlink(self.ip_args.get('bam', None), self.ip_bam)
        symlink(self.input_args.get('bam', None), self.input_bam)


    def get_bw_files(self):
        """
        Copy bw files
        """
        symlink(self.ip_args.get('bw', None), self.ip_bw)
        symlink(self.input_args.get('bw', None), self.input_bw)

        # ip over input, log2
        bwCompare(self.ip_bw, self.input_bw, self.ip_over_input_bw, 'log2',
            threads=self.threads, binsize=10)


    def call_peak(self):
        """
        Call peaks using MACS2
        ip, input
        ...
        """
        args_global = self.__dict__.copy()
        args_required = ['']
        args_local = dict((k, args_global[k]) for k in args_required if
            k in args_global)

        peak = Macs2(
            ip=self.ip_args.get('bam', None),
            control=self.input_args.get('bam', None),
            genome=self.genome,
            output=getattr(self, 'peak_dir', None),
            prefix=self.ip_args.get('smp_name', None),            
            genome_size=self.genome_size,
            **args_local)

        ## call peaks
        peak.callpeak()
        peak.bdgcmp(opt='ppois')
        peak.bdgcmp(opt='FE')
        peak.bdgcmp(opt='logLR')

        # ## annotation
        # peak.broadpeak_annotation()


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'chipseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'atac_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        run_shell_cmd(cmd) 


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # run
        self.get_bam_files()
        self.get_bw_files()
        self.call_peak()
        # qc
        self.report()


class ChIPseq(object):
    """
    Main port for ChIPseq analysis
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = ChIPseqConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def run_chipseq_design(self):
        """
        Create design
        """
        ChIPseqDesign(**self.__dict__).save()


    def run_chipseq_r1(self):
        """
        Run single
        """
        ChIPseqR1(**self.__dict__).run()


    def run_chipseq_rn(self):
        """
        Run multiple samples, rep_list
        """
        ChIPseqRn(**self.__dict__).run()


    def run_chipseq_rx_from_design(self):
        """
        Run ChIPseq all, require design

        multiple projects
        """
        design = Json(self.design).reader()

        for k, v in design.items():
            project_name = k
            args_local = v
            self = update_obj(self, args_local, force=True)
            self.run_chipseq_rx()


    def run_chipseq_rx(self):
        """
        main:
        Run whole pipeline

        required args:
        input_fq, ip_fq
        """
        # for ip fastq files
        ip_args = self.__dict__
        ip_local = {
            'is_ip': True,
            'fq1': self.ip,
            'fq2': self.ip_fq2,
            'build_design': None,
            'rep_list': None}
        ip_args.update(ip_local)
        ip = ChIPseqRn(**ip_args)

        # for input fastq
        input_args = self.__dict__
        input_local = {
            'is_ip': True,
            'fq1': self.input,
            'fq2': self.input_fq2,
            'build_design': None,
            'rep_list': None}
        input_args.update(input_local)
        input = ChIPseqRn(**input_args)

        # run
        ip.run()
        input.run()

        # for pipeline
        rx_args = self.__dict__
        rx_local = {
            'ip_dir': ip.project_dir,
            'input_dir': input.project_dir,
            'fq1': None,
            'fq2': None,
            'ip': None,
            'ip_fq2': None,
            'input': None,
            'input_fq2': None
        }
        rx_args.update(rx_local)
        rx = ChIPseqRx(**rx_args)
        rx.run()


    def run(self):
        """
        Run all
        """
        print('!AAAA-1', self.chipseq_type)
        if self.chipseq_type == 'build_design':
            self.run_chipseq_design()
        elif self.chipseq_type == 'chipseq_rx_from_design':
            self.run_chipseq_rx_from_design()
        elif self.chipseq_type == 'chipseq_r1':
            ChIPseqR1(self.__dict__).run()
        elif self.chipseq_type == 'chipseq_rn':
            ChIPseqRn(self.__dict__).run()
        elif self.chipseq_type == 'chipseq_rx':
            self.run_chipseq_rx()
        else:
            raise ValueError('unknown chipseq_type: {}'.format(self.chipseq_type))


class ChIPseqReader(object):
    """
    Read config.pickle from the local directory
    x
    |--config/
    |    |--arguments.pickle
    """
    def __init__(self, x):
        self.x = x
        self.get_config() # update
        self.chipseq_type = self.args.get('chipseq_type', None)

        # auto
        self.is_chipseq_r1 = self.chipseq_type == 'chipseq_r1'
        self.is_chipseq_rn = self.chipseq_type == 'chipseq_rn'
        self.is_chipseq_rx = self.chipseq_type == 'chipseq_rx'


    def get_config(self):
        """
        Read config 
        """
        config_pickle = os.path.join(self.x, 'config', 'arguments.pickle')
        if file_exists(config_pickle):
            with open(config_pickle, 'rb') as fh:
                self.args = pickle.load(fh)
        else:
            self.args = None


class ChIPseqDesign(object):
    """
    Prepare ip/input samples for ChIPseq analysis
    append/update

    format:
    ip: fq1, fq2
    input: fq1, fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'ip': None,
            'ip_fq2': None,
            'input': None,
            'input_fq2': None,
            'design': None,
            'append': True
        }
        self = update_obj(self, args_init, force=False)

        ## check
        self.design = file_abspath(self.design)
        if not isinstance(self.design, str):
            raise ValueError('--design, required')

        ## check fastq files
        log.info('Check fastq file existence:')
        self.ip = self.fq_input(self.ip)
        self.input = self.fq_input(self.input)
        
        if self.ip_fq2:
            self.ip_fq2 = self.fq_input(self.ip_fq2)
        if self.input_fq2:
            self.input_fq2 = self.fq_input(self.input_fq2)

        ## pairing
        self.fq_paired(self.ip, self.ip_fq2)
        self.fq_paired(self.input, self.input_fq2)


    def fq_input(self, fq):
        """
        Check input file, str or list
        
        convert to absolute path
        """
        if isinstance(fq, str):
            fq = [fq] # convert to list

        if isinstance(fq, list):            
            chk1 = file_exists(fq)
            for v, f in zip(chk1, fq):
                log.info('{} : {}'.format(v, f))

            if all(chk1):
                return file_abspath(fq)
            else:
                log.error('fastq files not exists')
        else:
            raise ValueError('unknown fastq, list expected, got {}'.format(type(fq).__name__))


    def fq_paired(self, fq1, fq2):
        """
        Check fq1, fq2, (list)
        """
        if fq2 is None:
            chk = True
        elif isinstance(fq1, str) and isinstance(fq2, str):
            chk = distance(fq1, fq2) == 1 and all(file_exists([fq1, fq2]))
        elif isinstance(fq1, list) and isinstance(fq2, list):
            chk = all([fq_paired(i, j) for i, j in zip(fq1, fq2)]) and all(file_exists(fq1 + fq2))
        else:
            chk = False

        if not chk:
            log.warning('fastq pairing failed')

        return chk


    def save(self, design=None):
        """
        Save args in Json format
        """
        # design: init
        if file_exists(self.design):
            design_dict = Json(self.design).reader()
        else:
            design_dict = {}

        # local
        args_local = {
            'design': self.design,
            'ip': self.ip,
            'ip_fq2': self.ip_fq2,
            'input': self.input,
            'input_fq2': self.input_fq2}

        # update design
        if self.append:
            key = 'chipseq_{:03d}'.format(len(design_dict) + 1)

            if args_local in list(design_dict.values()):
                log.warning('design exists, skipped ...')
            else:
                design_dict[key] = args_local
        else:
            key = 'chipseq_001'
            design_dict = dict((key, args_local))

        # save to file
        Json(design_dict).writer(self.design)


