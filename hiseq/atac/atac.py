"""
ATAC-seq pipeline


Working mode:

1. ATACseqSingle(): for single sample/fastq  
2. ATACseqMultiple(): for multiple sample/fastq files
3. ATACseqSample(): for single sample (could have multiple replicates)
4. ATACseq():

bam -> rmdup -> proper_paired -> peak/bw/...

##
AtacSingle()
AtacMultiple()
AtacReplicate()
Atac()

Input config file:

genome    outdir    name    fq1    fq2

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
 
####

compare replicates:

- bam_files
- bw_files
- peak
- motif
- report
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
import hiseq


def len2(x):
    n = 0
    if x is None:
        pass
    elif isinstance(x, str):
        n += 1
    elif isinstance(x, list):
        n += len(x)
    elif isinstance(x, dict):
        n += len(x)
    else:
        print('unknown type')
    return n


def print_df(x):
    # for k, v in x.items():
    #     print('{:30s} : {}'.format(k, v))
    for k in sorted(x.keys()):
        print('{:30s} : {}'.format(k, x[k]))


class Atac(object):
    """
    port for ATACseq
    """
    def __init__(self, **kwargs):
        """
        options:
        1. design: design.json,
        2. single/multiple: --fq1, --fq2, --genome, --outdir
        ...
        """
        self.update(kwargs, force=True) # fresh new
        local_config = AtacConfig(**self.__dict__)
        self.update(local_config.__dict__, force=True)
        self.pickle = None # terminate `pickle` option


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)    


    def run(self):
        """
        Run all RNAseq analysis
        """
        args = self.__dict__.copy()
        if self.atacseq_type == 'build_design':
            # RNAseqBuildDesign(**args).run()
            ATACseqFqDesign(**args)
        elif self.atacseq_type == 'atacseq_from_design':
            AtacFromDesign(**args).run()
        elif self.atacseq_type == 'atacseq_snrn':
            AtacSnRn(**args).run()
        elif self.atacseq_type == 'atacseq_s1rn':
            AtacS1Rn(**args).run()
        elif self.atacseq_type == 'atacseq_s1r1':
            AtacS1R1(**args).run()
            pass


class AtacConfig(object):
    """
    init args
    prepare files, folders
    determine mission
    ...
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True) # update from options
        self.init_args() # init args
        self.init_mission()


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)


    def init_args(self):
        """
        Check argument, defaults, conflicts
        Expect PE reads for ATACseq
        Also support SE reads
        """
        args_init = {
            'build_design': False,
            'design': None,
            'fq1': None,
            'fq2': None,
            'rep_list': None,
            'group': None,
            'genome': None,
            'outdir': str(pathlib.Path.cwd()),
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'unique_only': True,
            'genome_size': 0,
            'binsize': 50
        }
        self.update(args_init, force=False)

        # 1st level: create design.json
        if self.build_design:
            if self.fq1 is None:
                log.warning('-1, -2 required for build-design')

        # outdir
        self.outdir = file_abspath(self.outdir) # absolute path

        # RAM/CPU
        self.init_cpu()


    def init_cpu(self):
        """
        threads, CPUs
        """
        ## check number of threads, parallel_jobs
        ## parallel jobs * threads
        n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()

        max_jobs = int(n_cpu / 4.0)
        ## check parallel_jobs (max: 1/4 of n_cpus)
        if self.parallel_jobs > max_jobs: 
            log.warning('Too large, change parallel_jobs from {} to {}'.format(
                self.parallel_jobs, max_jobs))
            self.parallel_jobs = max_jobs

        ## check threads
        max_threads = int(0.8 * n_cpu / self.parallel_jobs)
        if self.threads * self.parallel_jobs > 0.8 * n_cpu:
            log.warning('Too large, change threads from {} to {}'.format(
                self.threads, max_threads))
            self.threads = max_threads 


    def init_fq(self):
        """
        Make sure
        fq1: None or list
        fq2: None or list
        exists
        """
        ####################################################
        ## 1. for fastq1 files                             #
        if self.fq1 is None:
            pass
        elif isinstance(self.fq1, list):
            self.fq1 = file_abspath(self.fq1)
        elif isinstance(self.fq1, str):
            self.fq1 = [file_abspath(self.fq1)]
        else:
            log.error('failed fq1, str, list expected, {} found'.format(type(self.fq1).__name__))

        ####################################################
        ## 2. for fastq2 files                             #
        if self.fq2 is None:
            pass
        elif isinstance(self.fq2, list):
            self.fq2 = file_abspath(self.fq2)
        elif isinstance(self.fq2, str):
            self.fq2 = [file_abspath(self.fq2)]
        else:
            log.error('failed fq2, None, list expected, {} found'.format(type(self.fq2).__name__))

        ####################################################
        ## 3. for rep_list                                 #
        if self.rep_list is None:
            pass
        elif isinstance(self.rep_list, list):
            self.rep_list = file_abspath(self.rep_list)
        elif isinstance(self.rep_list, str):
            self.fq1 = [file_abspath(self.rep_list)]
        else:
            log.error('failed, rep_list, None, list expected, {} found'.format(type(self.rep_list).__name__))

        ## required:
        if self.fq1 is None and self.rep_list is None:
            raise ValueError('fq1 or rep_list, required')
        n_smp = len(self.fq1) if isinstance(self.fq1, list) else len(self.rep_list)

        ####################################################
        ## 4. for group                                    #
        if self.group is None:
            if isinstance(self.fq1, list):
                self.group = fq_name_rmrep(self.fq1, pe_fix=True)
            elif isinstance(self.rep_list, list):
                self.group = [os.path.basename(i) for i in self.rep_list]
            else:
                raise ValueError('fq1 or rep_list required, when group=None')
                pass
        elif isinstance(self.group, list):
            if len(self.group) == 1:
                log.warning('assign all fq1/rep_list into one group')
                self.group = self.group * n_smp
            elif len(self.group) == n_smp:
                pass
            else:
                log.warning(self.group)
                log.warning(self.fq1)
                log.warning(self.rep_list)
                raise ValueError('group not match fq/rep_list')
        elif isinstance(self.group, str):
            log.warning('assign all fq1/rep_list into one group')
            self.group = [self.group] * n_smp
        else:
            raise ValueError('group: list expected, {} found'.format(type(self.group).__name__))

  
    def init_mission(self):
        """
        Determine the type of the ATACseq analysis:
        1. build_design
        2. atacSmpNRepN: n sample, n replicates
            - group > 1, 
        3. atacSmp1RepN: 1 sampele, n replicates
            - group = 1, fq1 > 1
        4. atacSmp1Rep1: 1 sample, 1 replicates
            - group = 1, fq1 = 1
        """
        create_dirs = getattr(self, 'create_dirs', True) # 

        atacseq_type = None
        # 1st level:
        if self.build_design:
            atacseq_type = 'build_design'
            self.init_build_design() # fresh start
        elif file_exists(self.design):
            atacseq_type = 'atacseq_from_design'
            self.init_atac_design(create_dirs)
        else:
            self.init_fq() # update, fq, rep_list, group
            if len2(set(self.group)) > 1:
                atacseq_type = 'atacseq_snrn'
                self.init_atac_snrn(create_dirs)
            elif len(set(self.group)) == 1:
                if len2(self.fq1) > 1 or len2(self.rep_list) > 1:
                    atacseq_type = 'atacseq_s1rn'
                    self.init_atac_s1rn(create_dirs)
                elif len2(self.fq1) == 1 or len2(self.rep_list) == 1:
                    atacseq_type = 'atacseq_s1r1'
                    self.init_atac_s1r1(create_dirs)
                else:
                    log.warning('unknown atacseq_type, group, fq1 required')
            else:
                log.warning('group required')

        # check
        if atacseq_type is None:
            raise ValueError('check, group, fq1, exit...')
        else:
            self.atacseq_type = atacseq_type


    def init_build_design(self, create_dirs=True):
        """
        Create design for ATACseq
        gather replicates into groups

        args:
        --build-design
        --design
        --fq1
        --fq2
        --group (optional)
        """
        print('!BBBB-1, init_build_design')
        if self.design is None:
            raise ValueError('--design required')
        if self.fq1 is None and self.rep_list is None:
            raise ValueError('fq1 or rep_list, required')
        pass


    def init_atac_design(self, create_dirs=True):
        """
        Read fq1/fq2 from design.json file
        
        assign fq1/fq2, groups
        """
        print('!BBBB-2, init_atac_design')
        if self.design is None:
            raise ValueError('--design is required')

        # path, files
        self.project_dir = self.outdir
        self.config_dir = os.path.join(self.project_dir, 'config')
        auto_files = {
            'report_dir': os.path.join(self.project_dir, 'report'),
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json')
        }
        self.update(auto_files)

        if create_dirs is True:
            check_path([
                self.project_dir,
                self.config_dir,
                self.report_dir])


    def init_atac_snrn(self, create_dirs=True):
        """
        input: group>1, fq1/rep_list

        for n sample, n replicates

        group > 1, rep_list/fq1 >= 1

        create dirs for groups/
        """
        print('!BBBB-3, init_atac_snrn')
        # group name
        group = getattr(self, 'group', None)
        if group is None:
            raise ValueError('group required')
        elif isinstance(group, list):
            if not len(group) > 1:
                raise ValueError('require multiple unique values in group')
        else:
            raise ValueError('group, expected list')

        # rep_list, >=1
        rep_list = getattr(self, 'rep_list', None)
        fq1 = getattr(self, 'fq1', None)
        if rep_list is None and fq1 is None:
            raise ValueError('--fq1, --rep-list, required')
        elif not rep_list is None:
            if not len(group) == len(rep_list):
                raise ValueError('Number of group and rep_list not identical')
        elif not fq1 is None:
            rep_list = [os.path.join(self.outdir, i) for i in fq_name(fq1, pe_fix=True)]
            if not len(group) == len(fq1):
                raise ValueError('Number of group and fq1 not identical')
        else:
            pass

        # path, files
        self.project_dir = self.outdir
        self.config_dir = os.path.join(self.project_dir, 'config')
        auto_files = {
            'report_dir': os.path.join(self.project_dir, 'report'),
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'rep_list': rep_list,
        }
        self.update(auto_files)

        if create_dirs is True:
            check_path([
                self.project_dir,
                self.config_dir,
                self.report_dir])


    def init_atac_s1rn(self, create_dirs=True):
        """
        input: rep_list, group>1

        for 1 sample, n replicates

        group = 1, rep_list
        
        ## pre-version: atac_merge
        """
        print('!BBBB-4, init_atac_s1rn')
        # group name
        group = getattr(self, 'group', None)
        if group is None:
            raise ValueError('group required')
        elif isinstance(group, list):
            if not len(set(group)) == 1:
                raise ValueError('require one unique value in group')
        else:
            raise ValueError('group, expected list')

        # rep_list, for replicates, same sample
        # rep_list: dirs of replicates
        self.rep_list = file_abspath(self.rep_list)
        # rep_list = getattr(self, 'rep_list', None)
        if self.rep_list is None and self.fq1 is None:
            raise ValueError('rep_list, required')

        # path, files
        self.project_name = group[0] # first one
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
            'report_dir': os.path.join(self.project_dir, 'report'),
            ## files
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'bam': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bam'),
            'peak': os.path.join(self.project_dir, 'peak', self.project_name + '_peaks.narrowPeak'),
            'bw': os.path.join(self.project_dir, 'bw_files', self.project_name + '.bigWig'),
            'lendist_txt': os.path.join(self.project_dir, 'qc', 'length_distribution.txt'),
            'lendist_pdf': os.path.join(self.project_dir, 'qc', 'length_distribution.pdf'),
            'frip_txt': os.path.join(self.project_dir, 'qc', 'FRiP.txt'),
            'cor_npz': os.path.join(self.project_dir, 'qc', 'cor.bam.npz'),
            'cor_counts': os.path.join(self.project_dir, 'qc', 'cor.bam.counts.tab'),
            'idr_txt': os.path.join(self.project_dir, 'qc', 'dir.txt'),
            'idr_log': os.path.join(self.project_dir, 'qc', 'dir.log'),
            'peak_overlap_pdf': os.path.join(self.project_dir, 'qc', 'peak_overlap_pdf')        
            }
        self.update(auto_files, force=True) # key

        if create_dirs:
            check_path([
                self.project_dir, 
                self.config_dir, 
                self.align_dir,
                self.bam_dir, 
                self.bw_dir, 
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir])


    def init_atac_s1r1(self, create_dirs=True):
        """
        input: fq1, group=1
        for 1 sample, 1 replicates

        group = 1, fq1 = 1
        
        ## pre-version: atac_single
        """
        print('!BBBB-5, init_atac_s1r1')
        # group name
        group = getattr(self, 'group', None)
        if group is None:
            raise ValueError('group required')
        elif isinstance(group, list):
            if not len(group) == 1:
                raise ValueError('require one unique value in group')
        else:
            raise ValueError('group, expected list')

        # fq1 = 1
        if isinstance(self.fq1, list) and len2(self.fq1) == 1:
            self.fq1 = self.fq1.pop()
        else:
            raise ValueError('fq1, required')

        if isinstance(self.fq2, list):
            self.fq2 = self.fq2.pop()

        # smp_name / error: init_atac_snrn
        smp_name = getattr(self, 'smp_name', None)
        if smp_name is None:
            smp_name = fq_name(self.fq1, pe_fix=True) # the first one
        self.smp_name = smp_name

        # path, files
        self.project_name = smp_name
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
            'report_dir': os.path.join(self.project_dir, 'report'),
            ## files
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'align_stat': os.path.join(self.project_dir, 'align', smp_name + '.align.txt'),
            # 'bam_raw': from align_dir, to-do
            'bam_rmdup': os.path.join(self.project_dir, 'bam_files', smp_name + '.rmdup.bam'),
            'bam_proper_pair': os.path.join(self.project_dir, 'bam_files', smp_name + '.proper_pair.bam'),
            'bam': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bam'),
            'bed': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bed'),
            'peak': os.path.join(self.project_dir, 'peak', self.project_name + '_peaks.narrowPeak'),
            'bw': os.path.join(self.project_dir, 'bw_files', self.project_name + '.bigWig'),
            'lendist_txt': os.path.join(self.project_dir, 'qc', 'length_distribution.txt'),
            'lendist_pdf': os.path.join(self.project_dir, 'qc', 'length_distribution.pdf'),
            'frip_txt': os.path.join(self.project_dir, 'qc', 'FRiP.txt'),
            }
        self.update(auto_files, force=True) # key

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
                self.report_dir])


class AtacFromDesign(object):
    """
    Run ATACseq, parse data from design
    group=n, fq1/rep_list=n

    ## ATACseqSnRn
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True)
        self.init_args()


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)

    
    def init_args(self):
        self_local = AtacConfig(**self.__dict__) # update
        self.update(self_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def run_group_single(self, g):
        """
        Run for single group
        """
        # input
        args_input = self.design_args.get(g, {})

        # for local
        args_local = self.__dict__.copy()
        args_init = {
            'design_args': None,
            'group': [g],
            'design': None,
            'build_design': None
        }
        args_local.update(args_input) # update fq1/rep_list/group
        args_local.update(args_init) # remove 

        config_local = AtacConfig(**args_local)

        AtacS1Rn(**config_local.__dict__).run()


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # n samples
        self.design_args = Json(self.design).reader()
        group_list = list(self.design_args.keys())

        # with Pool(processes=self.parallel_jobs) as pool:
        #     pool.map(self.run_group_single, group_list)

        # alternative
        for g in group_list:
            self.run_group_single(g)


class AtacSnRn(object):
    """
    Run ATACseq, parse data: group=n, fq1/rep_list=n
    ## ATACseqSnRn
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True)
        self.init_args()


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)

    
    def init_args(self):
        self_local = AtacConfig(**self.__dict__) # update
        self.update(self_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def pick_group_samples(self, g):
        """
        Pick samples for each group
        """
        if isinstance(g, str):
            g_index = [a for a, b in enumerate(self.group) if b == g]
            # fq1
            if isinstance(self.fq1, list):
                fq1 = [self.fq1[i] for i in g_index]
            else:
                fq1 = None

            # fq2
            if isinstance(self.fq2, list):
                fq2 = [self.fq2[i] for i in g_index]
            else:
                fq2 = None

            # rep_list
            if isinstance(self.rep_list, list):
                rep_list = [self.rep_list[i] for i in g_index]
            else:
                rep_list = None

            return {'fq1': fq1, 'fq2': fq2, 'rep_list': rep_list, 'group': g}
        elif isinstance(g, list):
            [self.pick_group_samples(i) for i in g]
        else:
            raise ValueError('group: str or list expected, {} found'.format(type(g).__name__))


    def run_group_single(self, g):
        """
        Run for single group
        """
        args_input = self.pick_group_samples(g)
        # for local
        args_local = self.__dict__.copy()
        args_init = {
            'build_design': None,
            'design': None,
            'group': [g]            
        }
        args_local.update(args_input) # update fq1/rep_list/group
        args_local.update(args_init) # remove 

        config_local = AtacConfig(**args_local)
        AtacS1Rn(**config_local.__dict__).run()


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # n samples
        group_list = self.group
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.run_group_single, group_list)

        # # alternative
        # for g in group_list:
        #     self.run_group_single(g)


class AtacS1Rn(object):
    """
    Run ATACseq, parse data group = 1, fq1 = n

    ## ATACseqS1Rn
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True)
        self.init_args()


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)

    
    def init_args(self):
        self_local = AtacConfig(**self.__dict__) # update
        self.update(self_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def get_bam_list(self):
        """
        get the proper_mapped.bam files
        """
        a = []
        for fq in self.fq1:
            fq_dir = os.path.join(self.outdir, fq_name(fq, pe_fix=True))
            a.append(AtacReader(fq_dir).args.get('bam'))
            # smp_name = fq_name(fq, pe_fix=True)
            # bam = os.path.join(self.outdir, smp_name, 'bam_files', smp_name + '.proper_pair.bam')
            # a.append(bam)

        return a


    def get_peak_list(self):
        """
        get the .narrowPeak files
        """
        p = []
        for fq in self.fq1:
            fq_dir = os.path.join(self.outdir, fq_name(fq, pe_fix=True))
            p.append(AtacReader(fq_dir).args.get('peak', None))
        
        return p


    def pick_fq_samples(self, i):
        """
        Pick samples for each group
        """
        if isinstance(i, int):
            # fq1
            if isinstance(self.fq1, list):
                fq1 = [self.fq1[i]]
            else:
                fq1 = None

            # fq2
            if isinstance(self.fq2, list):
                fq2 = [self.fq2[i]]
            else:
                fq2 = None

            # rep_list
            rep_list = None

            # group
            if isinstance(self.group, list):
                group = [self.group[i]]

            args_fq = {
                'fq1': fq1, 
                'fq2': fq2, 
                'rep_list': rep_list,
                'group': group}

            return(args_fq)
        else:
            raise ValueError('int expected, {} found'.format(type(g).__name__))


    def run_fq_single(self, i):
        """
        Run for single group
        """
        args_input = self.pick_fq_samples(i)
        # for local
        args_local = self.__dict__.copy()
        args_init = {
            'build_design': None,
            'design': None#,
            #'outdir': os.path.join(self.outdir, fq_name(args_input['fq1']).pop())
        }
        args_local.update(args_input) # update fq1/rep_list/group
        args_local.update(args_init) # remove 

        config_local = AtacConfig(**args_local)
        AtacS1R1(**config_local.__dict__).run()


    ##############################################
    ## merge replicates
    def mergebam(self):
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
            log.info('mergebam() skipped, file exists: {}'.format(
                self.bam))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('mergebam() failed.')


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


    def callpeak(self):
        """
        Call peaks using MACS2
        """
        args_peak = self.__dict__.copy()
        args_peak['genome_size'] = getattr(self, 'genome_size', 0)
        
        bam = self.bam
        bed = os.path.splitext(bam)[0] + '.bed'
        Bam(bam).to_bed(bed)
        genome = args_peak.pop('genome', None)
        output = args_peak.pop('peak_dir', None)
        prefix = args_peak.pop('group', None).pop()

        if check_file(self.peak):
            log.info('callpeak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, genome, output, prefix, atac=True, **args_peak).callpeak()


    def get_bam_cor(self, window=500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp
        
        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
            --outRawCounts *counts.tab -b bam
        """
        multiBamSummary = shutil.which('multiBamSummary')

        # run
        cmd = ' '.join([
            '{} bins --binSize {}'.format(multiBamSummary, window),
            '--smartLabels -o {}'.format(self.cor_npz),
            '--outRawCounts {}'.format(self.cor_counts)])

        if os.path.exists(self.cor_counts):
            log.info('bam_cor() skipped, file.exsits: {}'.format(
                self.cor_counts))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('bam_cor() failed.')


    def get_peak_overlap(self):
        """
        Compute the overlaps between overlaps
        """
        self.peak_list = self.get_peak_list()

        pkg_dir = os.path.dirname(hiseq.__file__)
        peak_overlapR = os.path.join(pkg_dir, 'bin', 'atac_peak_overlap.R')

        # run
        cmd = ' '.join([
            'Rscript',
            peak_overlapR,
            self.qc_dir] + self.peak_list)

        if os.path.exists(self.peak_overlap_pdf):
            log.info('get_peak_overlap() skipped, file exists: {}'.format(
                self.peak_overlap_pdf))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('get_peak_overlap() failed.')


    def get_rep_idr(self):
        """
        Calculate IDR for replicates
        1 vs 1
        peak files
        """
        idr = shutil.which('idr') # command

        # run
        cmd = ' '.join([
            'sort -k8,8nr -o {} {}'.format(self.peak_list[0], self.peak_list[0]),
            '&& sort -k8,8nr -o {} {}'.format(self.peak_list[1], self.peak_list[1]),
            '&& {} --input-file-type narrowPeak --rank p.value --plot'.format(idr),
            '--output-file {}'.format(self.idr_txt),
            '--log-output-file {}'.format(self.idr_log),
            '--samples',
            ' '.join(self.peak_list)])

        if os.path.exists(self.idr_txt):
            logging.info('rep_idr() skipped, file exists: {}'.format(
                self.idr_txt))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('rep_idr() failed.')


    def get_frip(self):
        """
        Save all FRiP.txt file to one
        """
        # get list
        self.frip_list = []
        for fq in self.fq1:
            fq_dir = os.path.join(self.outdir, fq_name(fq, pe_fix=True))
            self.frip_list.append(
                AtacReader(fq_dir).args.get('frip_txt', None))

        if not os.path.exists(self.frip_txt):
            with open(self.frip_txt, 'wt') as w:
                for f in self.frip_list:
                    with open(f) as r:
                        for line in r:
                            w.write(line)


    def get_align_txt(self):
        """
        Copy align.txt files to align/
        """
        pass


    def tss_enrich(self):
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
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atac_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'atac_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.outdir,
            self.report_dir)
        
        if check_file(atac_report_html):
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('report() failed.')


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # n samples
        i_list = list(range(len(self.fq1)))
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.run_fq_single, i_list)

        # # alternative
        # for i in i_list:
        #     self.run_fq_single(i)

        # run
        self.mergebam()
        self.bam_to_bw()
        self.callpeak()
        # qc
        self.get_bam_cor()
        self.get_peak_overlap()
        self.get_rep_idr()
        self.get_frip()
        # self.get_align_txt()
        # self.tss_enrich()
        self.report()


class AtacS1R1(object):
    """
    Run ATACseq, parse data group = 1, fq1 = 1

    bam -> rmdup -> proper_paired -> peak/bw/...

    ## ATACseqS1Rn
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True)
        self.init_args()


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)

    
    def init_args(self):
        self_local = AtacConfig(**self.__dict__) # update
        self.update(self_local.__dict__, force=True)

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
            if check_file(self.clean_fq_list):
                log.info('trim() skipped, file exists: {}'.format(
                    self.clean_fq_list))
            else:
                Trimmer(**args_trim).run()


    def align(self):
        """
        Alignment PE reads to reference genome, using bowtie2

        extra: -X 2000, ...
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'fq': fq1,
            'outdir': self.align_dir,
            'extra_para': '-X 2000',
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
            Bam(self.bam_rmdup).proper_pair(self.bam)
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


    def callpeak(self):
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
            log.info('callpeak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, self.genome, self.peak_dir, self.project_name, atac=True, **args_peak).callpeak()


    def qc_lendist(self):
        """
        Create length distribution plot
        """
        # create plot
        pkg_dir = os.path.dirname(hiseq.__file__)
        lendistR = os.path.join(pkg_dir, 'bin', 'atac_qc_lendist.R')
        cmd = 'Rscript {} {} {}'.format(
            lendistR,
            self.lendist_txt,
            self.lendist_pdf)

        if check_file(self.lendist_txt):
            log.info('qc_lendist() skipped, file exists: {}'.format(
                self.lendist_txt))
        else:
            # compute length distribution
            frag_length(self.bam, self.lendist_txt)
            # create plot
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('qc_lendist() failed.')


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            print("!XXXX " + self.bam)
            frip, n, total = cal_FRiP(self.peak, 
                self.bam)

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
        QC for ATACseq sample
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
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atac_report_single.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'atac_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.outdir,
            self.report_dir)
        
        if check_file(atac_report_html):
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('report() failed.')


    def run(self):
        """
        Run all steps for ATACseq pipeline
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
        self.callpeak()
        
        # 6. motif
        # self.motif()

        # 7. qc
        self.qc()

        # 8.report
        self.report()


class AtacReader(object):
    """
    Read config.txt/pickle
    """
    def __init__(self, x):
        self.x = x
        self.get_config() # update config_args (dict)
        self.atacseq_type = self.args.get('atacseq_type', None)

        # auto
        self.is_atac_s1r1 = self.atacseq_type == 'atacseq_s1r1'
        self.is_atac_s1rn = self.atacseq_type == 'atacseq_s1rn'
        self.is_atac_snrn = self.atacseq_type == 'atacseq_snrn'


    def get_config(self):
        """
        locate config.pickle file
        """
        config_pickle = os.path.join(self.x, 'config', 'arguments.pickle')
        if file_exists(config_pickle):
            with open(config_pickle, 'rb') as fh:
                self.args = pickle.load(fh)
        else:
            self.args = None


class ATACseqFqDesign(object):
    """
    Prepare fq samples for RNAseq analysis
    append/update

    format:
    - fq1, fq2, group
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.init_args() # init fq
        self.save() # to file


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)


    def init_args(self):
        args_default = {
            'fq1': None,
            'fq2': None,
            'rep_list': None,
            'group': None,
            'design': None,
            'append': True
        }
        self.update(args_default, force=False)

        # 1. check fq
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.rep_list = file_abspath(self.rep_list)
        
        # 2. check fq/rep_list
        if not self.fq1 is None:
            if not self.rep_list is None:
                log.info('ignore: rep_list')
            self.rep_list = None
            # for paired reads
            if not self.fq2 is None:
                tags = [self.fq_pair(a, b) for a, b in zip(self.fq1, self.fq2)]
                if not all(tags):
                    raise ValueError('fq1, fq2 not matched')

            # for group
            if self.group is None:
                n_group = fq_name_rmrep(self.fq1)
                self.group = n_group.pop()

        elif not self.rep_list is None:
            if not self.fq1 is None:
                log.info('ignore: fq1, fq2')
            self.fq1 = self.fq2 = None
            tags = [AtacReader(i).is_atac_single() for i in self.rep_list]
            if not all(tags):
                raise ValueError('rep_list, should be ATAC_single dir')

            # for group
            if self.group is None:
                n_group = fq_name_rmrep(self.rep_list)
                self.group = n_group.pop()
        else:
            raise ValueError('fq1 or rep_list required')

        # 3. design
        if self.design is None:
            self.design = self._tmp()
            log.warning('design, create a temp file: ' + self.design)
        # convert to absolute path
        self.design = file_abspath(self.design)

        # 4. group
        if isinstance(self.group, list):
            log.warning('choose the 1st one as group: {}'.format(self.group.pop()))
            self.group = self.group.pop()
        elif isinstance(self.group, str):
            pass
        else:
            raise ValueError('group, str expected, {} found'.format(type(self.group).__name__))

        if self.group is None:
            if not self.fq1 is None:
                self.group = fq_name_rmrep(self.fq1)
            else:
                self.group = fq_name_rmrep(self.rep_list)


    def fq_pair(self, fq1, fq2, n_diff=1):
        """
        read1 and read2 
        should be in the same name, _1, _2
        """
        if isinstance(fq1, str) and isinstance(fq2, str):
            if distance(fq1, fq2) <= n_diff:
                return True


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json',
            delete=False)
        return tmp.name


    def save(self):
        """
        save config
        """
        # json file
        if not file_exists(self.design):
            args_pre = {}
        else:
            args_pre = Json(self.design).reader()

        # save only specific args
        dd = {
            'fq1': self.fq1,
            'fq2': self.fq2,
            'rep_list': self.rep_list,
            'design': self.design,
            'group': self.group
        }

        # new data
        # append or update
        if self.append:
            # whether exists
            if self.group in list(args_pre.keys()):
                log.warning('design exists, skipped ... {}'.format(self.group))
            else:
                args_pre[self.group] = dd # append
        else:
            args_pre = {} # empty
            args_pre[self.group] = dd #

        # save
        Json(args_pre).writer(self.design)

