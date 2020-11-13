"""
CUT&RUN pipeline

workflow

1. trimming
double-trimming
a. full length adapter 
b. <6bp adapters

2. Alignment
dovetail alignment using bowtie2

--dovetail 

--local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700

## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt

optional: 
spike-in
chrM

3. remove duplicates
remove or not

picardCMD MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam 
O=$projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam 
REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.rmDup.txt



4. peak calling


5. cut metrix


Quality control

1. fragment size (expected < 120bp)
2. adapter content pct (10-15%)
3. read duplication (10-15%)
4. alignment pct (>90%)
5. peaks ()
6. enrichment of motif 
7. 


Directory structure:

raw_data
clean_data
align
spike-in
  - align
  - scale
with_dup
  - bam_files
  - peak
  - motif
rm_dup
  - bam_files
  - peak
  - motif
qc
  - fastqc
  - library size
  - adapter pct
  - align pct
  - frag size 
  - dup pct
  - replicate cor
  - peaks
  - peaks overlap 
  - motif
  - overall: library size, duplicate rate%, unique mapping%

scale_factor
  - bam to bw

report


2020-11-12

+ Add arg: keep_tmp (default: False)

"""


import os
import re
import glob
from multiprocessing import Pool
from Levenshtein import distance
from hiseq.utils.helper import *
from hiseq.trim.trimmer import Trimmer
from hiseq.align.alignment import Alignment, AlignIndex
from hiseq.peak.call_peak import Macs2
from hiseq.atac.atac_utils import *
from hiseq.fragsize.fragsize import BamPEFragSize
import hiseq
import collections


def print_dict(d):
    d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


def file_symlink(src, dest, absolute_path=False):
    """
    Create symlink for dest files
    
    file, directories
    """
    if isinstance(src, str) and isinstance(dest, str):
        src = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        dest = os.path.abspath(os.path.expanduser(os.path.expandvars(dest)))

        # target exists
        if os.path.exists(dest):
            log.info('symlink() skipped, dest exists: {}'.format(dest))
        elif os.path.isdir(src):
            if absolute_path:
                os.symlink(src, dest)
            else:
                rel_src = os.path.relpath(src, dest)
                os.symlink(rel_src, dest)

        elif os.path.isfile(src):  
            src_filename = os.path.basename(src)
            src_dirname = os.path.dirname(src)
            dest_dirname = os.path.dirname(dest)

            if absolute_path:
                os.symlink(src, dest)
            else:
                rel_src_dirname = os.path.relpath(src_dirname, dest_dirname)
                rel_src = os.path.join(rel_src_dirname, src_filename)
                os.symlink(rel_src, dest)

        else:
            log.error('src path should be file or directory')

    else:
        log.error('src path should be string, got: {}'.format(
            type(src).__name__))


def path_remove(x, ask=True):
    """
    Remove directories
    """
    del_list = []
    undel_list = []
    if isinstance(x, str):
        if os.path.isdir(x) and file_exists(x):
            del_list.append(x)
        else:
            undel_list.append(x)
    elif isinstance(x, list):
        for f in x:
            if os.path.isdir(f) and file_exists(x):
                del_list.append(f)
            else:
                undel_list.append(f)
    elif isinstance(x, dict):
        for f in list(x.values()):
            if os.path.isdir(f) and file_exists(x):
                del_list.append(f)
            else:
                undel_list.append(f)
    else:
        log.info('Nothing removed, str, list, dict expected, got {}'.format(
            type(x).__name__))

    # remove files
    if len(del_list) > 0:
        del_msg = ['{:>6s}: {}'.format('remove', i) for i in del_list]
        undel_msg = ['{:>6s}: {}'.format('skip', i) for i in undel_list]
        msg = '\n'.join(del_msg + undel_msg)
        log.info('Removing files: \n' + msg)
        
        if ask:
            ask_msg = input('Removing the files? [Y|n]： ')
        else:
            ask_msg = 'Y'

        # remove
        if ask_msg.lower() in ['y', 'yes']:
            for f in del_list:
                # os.remove(f)
                shutil.rmtree(f)
            log.info('{} files removed'.format(len(del_list)))
        else:
            log.info('Nothing removed, skipped')


def file_remove(x, ask=True):
    """
    Remove files, directories
    """
    del_list = []
    undel_list = []
    if isinstance(x, str):
        if os.path.isfile(x) and file_exists(x):
            del_list.append(x)
        else:
            undel_list.append(x)
    elif isinstance(x, list):
        for f in x:
            if os.path.isfile(f) and file_exists(x):
                del_list.append(f)
            else:
                undel_list.append(f)
    elif isinstance(x, dict):
        for f in list(x.values()):
            if os.path.isfile(f) and file_exists(x):
                del_list.append(f)
            else:
                undel_list.append(f)
    else:
        log.info('Nothing removed, str, list, dict expected, got {}'.format(
            type(x).__name__))

    # remove files
    if len(del_list) > 0:
        del_msg = ['{:>6s}: {}'.format('remove', i) for i in del_list]
        undel_msg = ['{:>6s}: {}'.format('skip', i) for i in undel_list]
        msg = '\n'.join(del_msg + undel_msg)
        log.info('Removing files: \n' + msg)
        
        if ask:
            ask_msg = input('Removing the files? [Y|n]： ')
        else:
            ask_msg = 'Y'

        # remove
        if ask_msg.lower() in ['y', 'yes']:
            for f in del_list:
                os.remove(f)
            log.info('{} files removed'.format(len(del_list)))
        else:
            log.info('Nothing removed, skipped')


def file_copy(src, dest, force=False):
    """
    copy files to dest
    """
    if src is None:
        log.error('file_copy() skipped, src is NoneType')
    elif os.path.isfile(src):
        if os.path.isdir(dest):
            dest_file = os.path.join(dest, os.path.basename(src))
            shutil.copy(src, dest_file)
        elif os.path.isfile(dest):
            if force:
                shutil.copy(src, dest)
            else:
                log.error('file_copy() skipped, dest exists: {}'.format(dest))
        elif os.path.exists(os.path.dirname(dest)):
            shutil.copy(src, dest)
        else:
            log.error('file_copy() skipped, dest is not file or dir')


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


def dict2pickle(d, p, update=False):
    """Save dict as pickle
    d is dict
    p is pickle file
    """
    if isinstance(d, dict) and isinstance(p, str):
        if file_exists(p):
            with open(p, 'rb') as r:
                d_old = pickle.load(r)
        else:
            d_old = {}

        # update
        d_old.update(d)
        d_new = d_old if update else d

        # save to file
        with open(p, 'wb') as w:
            pickle.dump(d_new, w, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        log.error('dict2pickle() failed, dict, str expect')


def pickle2dict(p):
    """
    Read pickle, save as dict
    """
    d = {}
    if file_exists(p):
        with open(p, 'rb') as r:
            d = pickle.load(r)
    else:
        log.error('pickle file illegal')

    return d


def check_fq(fq):
    """
    Make sure
    fq: str or list, or None
    """
    if fq is None:
        pass
    elif isinstance(fq, list):
        fq = file_abspath(fq)
    elif isinstance(fq, str):
        fq = [file_abspath(fq)]
    else:
        log.error('fq failed, Nont, str, list expected, got {}'.format(
            type(fq).__name__))

    return fq


def fq_paired(fq1, fq2):
    """
    Make sure fq1 and fq2, proper paired
    """
    # convert fq to list
    fq1 = check_fq(fq1)
    fq2 = check_fq(fq2)

    flag = False
    if isinstance(fq1, list) and isinstance(fq2, list):
        if len(fq1) == len(fq2):
            flag = all([distance(i, j) == 1 for i, j in zip(fq1, fq2)])

    return flag


class Align(object):
    """
    Alignment reads for CUT&RUN, CUT&TAG
    suppose paired end reads
    --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700
    
    Input: fq, genome_index, outputdir, stat/log
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_index()
        self.init_files()        
        self.prep_cmd()


    def init_args(self):
        """
        ATACseq, max_fragment = 2000
        CUT&RUN, max_fragment = 700
        CUT&TAG, max_fragment = 700
        """
        default_args = {
            'fq1': None,
            'fq2': None,
            'index': None,
            'genome': None,
            'extra_index': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'threads': 1,
            'overwrite': False,
            'max_fragment': 700,
            'remove_unmap': False,
            'hiseq_type': 'hiseq_alignment',
            'keep_tmp': False
            }
        self = update_obj(self, default_args, force=False)

        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
            log.warning('set pwd as outdir: {}'.format(self.outdir))
        self.outdir = file_abspath(self.outdir)
        # check_path(self.outdir)

        # fq1
        if isinstance(self.fq1, str):
            if check_file(self.fq1):
                self.fq1 = file_abspath(self.fq1)
                self.fq_check = 0
            else:
                log.error('--fq1, file not exists, {}'.format(self.fq1))
                self.fq_check = 1
        else:
            log.error('--fq1, str expected, got {}'.format(
                type(self.fq1).__name__))
            self.fq_check = 1

        # fq2
        if isinstance(self.fq2, str):
            if check_file(self.fq2):
                self.fq2 = file_abspath(self.fq2)
                self.fq_check = 0
            else:
                log.error('--fq2, file not exists, {}'.format(self.fq1))
                self.fq_check = 1
        else:
            log.error('--fq2, str expected, got {}'.format(
                type(self.fq2).__name__))
            self.fq_check = 1

        # prefix
        self.smp_name = fq_name(self.fq1, pe_fix=True)


    def init_index(self):
        """
        Determine the index
        """
        if isinstance(self.extra_index, str):
            self.index = self.extra_index
            self.index_check = 0 if AlignIndex(
                index=self.extra_index, aligner=self.aligner).is_index() else 1
        elif isinstance(self.index, str):
            self.index_check = 0 if AlignIndex(
                index=self.index, aligner=self.aligner).is_index() else 1
        elif isinstance(self.genome, str):
            self.index = AlignIndex(aligner=self.aligner).search(
                genome=self.genome, group = 'genome')
            self.index_check = 1 if self.index is None else 0
        else:
            log.error('-g, --index, failed')
            self.index_check = 1

        # index name
        if isinstance(self.index, str):
            self.index_name = AlignIndex(index=self.index).index_name()
        else:
            self.index_name = None


    def init_files(self):
        """
        default files for output
        bam, sam, log, unmap, ...
        """
        # subdir
        subdir = os.path.join(self.outdir, self.smp_name, self.index_name)
        check_path(subdir)

        # output files
        prefix = os.path.join(subdir, self.smp_name)
        default_files = {
            'subdir': subdir,
            'config_pickle': os.path.join(subdir, 'config.pickle'),
            'cmd_shell': os.path.join(subdir, 'cmd.txt'),
            'bam': prefix + '.bam',
            'sam': prefix + '.sam',
            'unmap': prefix + '.unmap.fastq',
            'unmap1': prefix + '.unmap.1.fastq',
            'unmap2': prefix + '.unmap.2.fastq',
            'align_log': prefix + '.bowtie2.log',
            'align_stat': prefix + '.bowtie2.stat',
            'align_json': prefix + '.bowtie2.json',
            'align_flagstat': prefix + '.flagstat',
            'rmdup_bam': prefix + '.rmdup.bam',
            'rmdup_matrix': prefix + '.rmdup.matrix.txt'
        }
        self = update_obj(self, default_files, force=True)


    def prep_cmd(self):
        """
        Alignment for CUT&RUN, CUT&TAG, ATACseq, ...
        """
        self.cmd = ' '.join([
            '{}'.format(shutil.which('bowtie2')),
            '-p {}'.format(self.threads),
            '--local --very-sensitive --no-mixed --no-discordant -I 10',
            '-X {}'.format(self.max_fragment),
            '-x {}'.format(self.index),
            '-1 {} -2 {}'.format(self.fq1, self.fq2),
            '--un-conc {}'.format(self.unmap),
            '1> {} 2> {}'.format(self.sam, self.align_log),
            '&& samtools view -@ {} -bS'.format(self.threads),
            '<(samtools view -H {} ;'.format(self.sam),
            "samtools view -F 2048 {} | grep 'YT:Z:CP')".format(self.sam),
            '| samtools sort -@ {} -o {} -'.format(self.threads, self.bam),
            '&& samtools index {}'.format(self.bam),
            '&& samtools flagstat {} > {}'.format(self.bam, self.align_flagstat),
            '&& [[ -f {} ]] && rm {}'.format(self.sam, self.sam)
            ])
        
        # save cmd txt
        cmd_txt = self.subdir + '/cmd.txt'
        with open(cmd_txt, 'wt') as w:
            w.write(self.cmd + '\n')


    def parse_align(self, to_json=None):
        """
        Parsing the log of bowtie2
        """
        """
        Wrapper bowtie2 log
        Bowtie2:

        SE:
        10000 reads; of these:
          10000 (100.00%) were unpaired; of these:
            166 (1.66%) aligned 0 times
            2815 (28.15%) aligned exactly 1 time
            7019 (70.19%) aligned >1 times
        98.34% overall alignment rate

        PE:
        100000 reads; of these:
          100000 (100.00%) were paired; of these:
            92926 (92.93%) aligned concordantly 0 times
            5893 (5.89%) aligned concordantly exactly 1 time
            1181 (1.18%) aligned concordantly >1 times
            ----
            92926 pairs aligned concordantly 0 times; of these:
              1087 (1.17%) aligned discordantly 1 time
            ----
            91839 pairs aligned 0 times concordantly or discordantly; of these:
              183678 mates make up the pairs; of these:
                183215 (99.75%) aligned 0 times
                101 (0.05%) aligned exactly 1 time
                362 (0.20%) aligned >1 times
        8.39% overall alignment rate

        unique, multiple, unmap, map, total
        """
        self.unique_only = False

        dd = {}
        se_tag = 1 #
        with open(self.align_log, 'rt') as r:
            for line in r:
                value = line.strip().split(' ')[0]
                if '%' in value:
                    continue
                if line.strip().startswith('----'):
                    continue
                value = int(value)

                ## paired tag
                if 'were paired; of these' in line:
                    dd['total'] = value
                    se_tag = 0
                elif 'aligned concordantly 0 times' in line:
                    dd['unmap'] = value
                    se_tag = 0
                elif 'aligned concordantly exactly 1 time' in line:
                    dd['unique'] = value
                    se_tag = 0
                elif 'aligned concordantly >1 times' in line:
                    dd['multiple'] = value
                    se_tag = 0
                elif 'reads; of these' in line and se_tag:
                    dd['total'] = value
                elif 'aligned 0 times' in line and se_tag:
                    dd['unmap'] = value
                elif 'aligned exactly 1 time' in line and se_tag:
                    dd['unique'] = value
                elif 'aligned >1 times' in line and se_tag:
                    dd['multiple'] = value
                else:
                    pass

        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']

        # save fqname, indexname,
        dd['fqname'] = self.smp_name
        dd['index_name'] = self.index_name
        self.log_dict = dd

        # save dict to plaintext file
        with open(self.align_stat, 'wt') as w:
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 'fqname', 
                'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')

        ## save to json
        if to_json:
            Json(dd).writer(self.align_json)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        if self.fq_check > 0:
            raise ValueError('--fq1, --fq2, illegal')
        if self.index_check > 0:
            raise ValueError('--index, --genome, illegal')

        # prepare dirs
        check_path(self.subdir)

        # save config
        dict2pickle(self.__dict__, self.config_pickle)

        # run cmd
        if file_exists(self.bam) and not self.overwrite:
            log.info('align() skipped, file exists:{}'.format(self.bam))
        else:
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('align() failed, check {}'.format(self.align_log))

        # read log
        if file_exists(self.align_log):
            self.parse_align(to_json=True)

        # temp files
        del_list = [self.sam, self.unmap1, self.unmap2]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)


class CallPeak(object):
    """
    Call peaks using MACS2, or SEACR

    MACS2: 
    macs2 callpeak -t ip.bam -c input.bam \
      -g hs -f BAMPE -n macs2_peak_q0.1 \
      --outdir outdir -q 0.1 
      --keep-dup all 
      2>log.txt
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        default_args = {
            'method': 'macs2',
            'ip': None,
            'input': None,
            'genome': None,
            'genome_size': None,
            'prefix': None,
            'outdir': None,
            'overwrite': False,
            'hiseq_type': 'macs2_call_peak',
            'keep_tmp': False
            }
        self = update_obj(self, default_args, force=False)

        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
            log.warning('set pwd as outdir: {}'.format(self.outdir))
        self.outdir = file_abspath(self.outdir)
        check_path(self.outdir)

        # bam files
        if not file_exists(self.ip):
            raise ValueError('ip, required')

        self.ip = file_abspath(self.ip)
        self.input = file_abspath(self.input)

        # genome
        g = {
          'dm6': 'dm',
          'dm3': 'dm',
          'mm9': 'mm',
          'mm10': 'mm',
          'hg19': 'hs',
          'hg38': 'hs',
          'GRCh38': 'hs'
        }

        if isinstance(self.genome, str):
            self.genome_size = g.get(self.genome, None)
            if self.genome_size is None:
                log.error('genome, unknown, {}'.format(self.genome))

        elif isinstance(self.genome_size, int):
            if self.genome_size < 1:
                log.error('genome_size, failed, {}'.format(self.genome_size))

        else:
            pass

        self.ip_name = os.path.splitext(os.path.basename(self.ip))[0]
        if self.input is None:
            self.input_name = ''
        else:
            self.input_name = os.path.splitext(os.path.basename(self.input))[0]
        
        # prefix
        if isinstance(self.prefix, str):
            pass
        else:
            self.prefix = self.ip_name

        self.prefix_top001 = self.prefix + '.top0.01'


    def init_files(self):
        # files
        default_files = {
            'config_pickle': self.outdir + '/config.pickle',
            'macs2_peak': self.outdir + '/' + self.prefix + '_peaks.narrowPeak',
            'macs2_log': self.outdir + '/' + self.prefix + '.callpeak.log',
            'macs2_peak_xls': self.outdir + '/' + self.prefix + '_peaks.xls',
            'ip_bed': self.outdir + '/' + self.ip_name + '.frag.bed',
            'ip_bg': self.outdir + '/' + self.ip_name + '.frag.bg',
            'input_bed': self.outdir + '/' + self.input_name + '.frag.bed',
            'input_bg': self.outdir + '/' + self.input_name + '.frag.bg',
            'seacr_peak': self.outdir + '/' + self.prefix + '.stringent.bed',
            'seacr_peak_top': self.outdir + '/' + self.prefix_top001 + '.stringent.bed'
        }
        self = update_obj(self, default_files, force=True) # key

        ## gsize
        self.genome_size_file = Genome(genome=self.genome).get_fasize()
        

    def run_macs2(self):
        input_arg = '' if self.input is None else '-c {}'.format(self.input)

        cmd = ' '.join([
            '{} callpeak'.format(shutil.which('macs2')),
            '-t {} {}'.format(self.ip, input_arg),
            '-g {} -f BAMPE'.format(self.genome_size),
            '-n {}'.format(self.outdir + '/' + self.prefix),
            '-q 0.1 --keep-dup all',
            '--nomodel --extsize 150',
            '2> {}'.format(self.macs2_log)
            ])

        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.macs2.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        # run
        if file_exists(self.macs2_peak) and not self.overwrite:
            log.info('run_macs2() skipped, file exists:{}'.format(
                self.macs2_peak))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('run_macs2() failed, check {}'.format(
                    self.macs2_log))


    def bampe_to_bg(self, bam, bg):
        """
        Convert bam to bed
        bedpe, sort -n (by name)
        """
        bam_sorted = os.path.splitext(bg)[0] + '.sorted_by_name.bam'
        bed = os.path.splitext(bg)[0] + '.bed'

        cmd = ' '.join([
            'samtools sort -@ 4 -n -o {} {}'.format(bam_sorted, bam),
            '&& bedtools bamtobed -bedpe -i {}'.format(bam_sorted),
            "| awk '$1==$4 && $6-$2 < 1000 {print $0}'",
            '| cut -f 1,2,6',
            '| sort -k1,1 -k2,2n > {}'.format(bed),
            '&& bedtools genomecov -bg -i {}'.format(bed),
            '-g {}'.format(self.genome_size_file),
            '> {}'.format(bg)
        ])

        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.bam2bg.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        # run
        if file_exists(bg) and not self.overwrite:
            log.info('bampe_to_bg() skipped, file exists:{}'.format(bg))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('bampe_to_bg() failed, check {}'.format(bg))

        # temp files
        del_list = [bam_sorted, bed]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)


    def run_seacr(self):
        """
        bash $seacr ip.bedgraph input.bedgraph \
        non stringent output
    
        see: https://github.com/FredHutch/SEACR
        1. bam-to-bed
        bedtools bamtobed -bedpe -i in.bam > out.bed
        awk '$1==$4 && $6-$2 < 1000 {print $0}' out.bed > out.clean.bed
        cut -f1,2,6 out.clean.bed | sort -k1,1 -k2,2n -k3,3n > out.fragments.bed
        bedtools genomecov -bg -i out.fragments.bed -g genome > out.fragments.bg
        
        2. call peak
        bash seacr out.fragments.bg IgG.bg non stringent output
        bash seacr out.fragments.bg 0.01 non 0.01 non stringent top0.01.peaks
        """
        self.bampe_to_bg(self.ip, self.ip_bg)
        cmd = ' '.join([
            'bash {}'.format(shutil.which('SEACR_1.3.sh')),
            '{} 0.01'.format(self.ip_bg),
            'non stringent {}'.format(self.outdir + '/' + self.prefix_top001)
        ])

        if not self.input is None:
            self.bampe_to_bg(self.input, self.input_bg)

            cmd = ' '.join([
                '&& {}'.format(cmd),
                'bash {}'.format(shutil.which('SEACR_1.3.sh')),
                '{} {}'.format(self.ip_bg, self.input_bg),
                'non stringent {}'.format(self.outdir + '/' + self.prefix)            
            ])

        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.SEACR.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        # run
        if file_exists(self.seacr_peak_top) and not self.overwrite:
            log.info('run_seacr() skipped, file exists:{}'.format(
                self.seacr_peak))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('run_seacr() failed, check {}'.format(
                    self.seacr_peak))


    def run(self):
        # save config
        dict2pickle(self.__dict__, self.config_pickle)

        # generate cmd
        if self.method.lower() == 'macs2':
            self.run_macs2()
        elif self.method.lower() == 'seacr':
            self.run_seacr()
        else:
            raise ValueError('unknown method: {macs2|seacr}, got {}'.format(
                self.method))

        # remove files
        del_list = []


class CnRConfig(object):
    """
    Global config for CnR analysis
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()
        self.init_fq()
        self.init_mission()


    def init_args(self):
        """
        required arguments for ChIPseq analysis
        """
        args_init = {
            'build_design': False,
            'design': None,
            'rep_list': None,
            'ip_dir': None,
            'input_dir': None,
            'ip': None,
            'input': None,
            'ip_fq2': None,
            'input_fq2': None,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'spikein': None,
            'spikein_index': None,
            'extra_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'trimmed': True,
            'keep_tmp': False
        }
        self = update_obj(self, args_init, force=False)

        if self.outdir is None:
          self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # aligner
        if self.aligner is None:
            self.aligner = 'bowtie2'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, 
            self.parallel_jobs)


    def init_files(self, create_dirs=True):
        """
        Prepare directories, files
        """
        # path, files
        self.project_dir = self.outdir
        self.config_dir = self.project_dir + '/config'

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json'
        }
        self = update_obj(self, default_files, force=True) # key

        if create_dirs:
            check_path(self.config_dir)


    def init_fq(self):
        """
        Check fastq files
        IP, Input
        """
        chk1 = fq_paired(self.ip, self.ip_fq2)
        chk2 = fq_paired(self.input, self.input_fq2)


    def init_mission(self):
        """
        Determine the type of CnR analysis
        1. build_design
        2. Single: IP/Input pair
        3. Multiple: IP/Input pair
        4. x, ip vs input
        """
        self.hiseq_type = None

        # 1st level
        if self.build_design:
            self.hiseq_type = 'build_design'
        elif file_exists(self.design):
            self.hiseq_type = 'hiseq_from_design'
        elif isinstance(self.ip_dir, str) and isinstance(self.input_dir, str):
            self.hiseq_type = 'hiseq_rx'
        elif all([not i is None for i in [self.ip, self.input]]):
            self.hiseq_type = 'hiseq_rx'
        elif isinstance(self.rep_list, list) or isinstance(self.fq1, list):
            self.hiseq_type = 'hiseq_rn'
        elif isinstance(self.fq1, str):
            self.hiseq_type = 'hiseq_r1'
        else:
            raise ValueError('unknown CnR type')

        # check
        if self.hiseq_type is None:
            raise ValueError('unknown CnR')


class CnRxConfig(object):
    """
    Prepare directories for Rx: IP vs Input

    require: ip_dir, input_dir
    or:
    ip,ip-fq2, input,input-fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_fq()
        # self.init_index()
        self.init_files()


    def init_args(self):
        """
        required arguments for CnR analysis
        """
        args_init = {
            'outdir': None,
            'ip_dir': None,
            'input_dir': None,
            'ip': None,
            'ip_fq2': None,
            'input': None,
            'input_fq2': None,
            'smp_name': None,
            'genome': None,
            'spikein': None,
            'spikein_index': None,
            'aligner': 'bowtie2',
            'extra_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'gene_bed': None,
            'keep_tmp': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'hiseq_rx' # 

        # aligner
        if self.aligner is None:
            self.aligner = 'bowtie2'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, 
            self.parallel_jobs)

        self.outdir = file_abspath(self.outdir)


    def init_fq(self):
        """
        Support fastq files
        ip, input
        """
        if isinstance(self.ip, list) and isinstance(self.input, list):
            self.ip = file_abspath(self.ip)
            self.ip_fq2 = file_abspath(self.ip_fq2)
            self.input = file_abspath(self.input)
            self.input_fq2 = file_abspath(self.input_fq2)

            # file exists
            if not all(file_exists(self.ip)):
                raise ValueError('--ip, file not exists: {}'.format(
                    self.ip))

            if not all(file_exists(self.input)):
                raise ValueError('--input, file not exists: {}'.format(
                    self.ip))

            # check, ip,input name
            ip_names = set(map(fq_name_rmrep, self.ip))
            input_names = set(map(fq_name_rmrep, self.input))
            if len(ip_names) > 1:
                raise ValueError('--ip failed, filename differ: {}'.format(
                    ip_names))

            if len(input_names) > 1:
                raise ValueError('--input failed, filename differ: {}'.format(
                    input_names))

            # check paired
            ip_pe = fq_paired(self.ip, self.ip_fq2)
            input_pe = fq_paired(self.input, self.input_fq2)
            if not ip_pe:
                raise ValueError('--ip, --ip-fq2, not paired: {}, {}'.format(
                    self.ip, self.ip_fq2))

            if not input_pe:
                raise ValueError('--input, --input-fq2, not paired: \
                    {}, {}'.format(self.input, self.input_fq2))

            # update ip_dir, input_dir
            self.ip_name = ip_names.pop()
            self.input_name = input_names.pop()
            self.ip_dir = self.outdir + '/' + self.ip_name
            self.input_dir = self.outdir + '/' + self.input_name

        elif isinstance(self.ip_dir, str) and isinstance(self.input_dir, str):
            self.ip_dir = file_abspath(self.ip_dir)
            self.input_dir = file_abspath(self.input_dir)
            
            c_ip = CnRReader(self.ip_dir)
            c_input = CnRReader(self.input_dir)

            if not c_ip.is_hiseq_rn or not c_input.is_hiseq_rn:
                raise ValueError('ip_dir, input_dir failed: {}, {}'.format(
                    self.ip_dir, self.input_dir))

            self.ip_name = c_ip.args.get('smp_name', None)
            self.input_name = c_input.args.get('smp_name', None)

        else:
            raise ValueError('--ip, --input, or --ip-dir, --input-dir failed')

        # update smp_name
        self.smp_name = '{}.vs.{}'.format(self.ip_name, self.input_name)


    def init_files(self, create_dirs=True):
        """
        Prepare directories, files
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = self.outdir + '/' + self.project_name
        self.config_dir = self.project_dir + '/config'

        default_dirs = {
            'bam_dir': 'bam_files',
            'bw_dir': 'bw_files',
            'peak_dir': 'peak',
            'motif_dir': 'motif',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json',
            'ip_bam': self.bam_dir + '/' + self.ip_name + '.bam',
            'ip_bw': self.bw_dir + '/' + self.ip_name + '.bigWig',
            'input_bam': self.bam_dir + '/' + self.input_name + '.bam',
            'input_bw': self.bw_dir + '/' + self.input_name + '.bigWig',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',

            'ip_over_input_bw': self.bw_dir + '/' + self.ip_name + '.ip_over_input.bigWig',
            'bw': self.bw_dir + '/' + self.ip_name + '.bigWig',
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_fingerprint': self.qc_dir + '/09.fingerprint.png'
        }
        self = update_obj(self, default_files, force=True) # key

        if create_dirs:
            check_path([
                self.project_dir, 
                self.config_dir, 
                self.bam_dir, 
                self.bw_dir, 
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir])


class CnRnConfig(object):
    """
    Prepare directories for R1 single replicates

    require: fq1/fq2, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_fq()
        self.init_index()
        self.init_files()


    def init_args(self):
        """
        required arguments for CnR analysis
        """
        args_init = {
            'is_ip': True,
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'smp_name': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'extra_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'keep_tmp': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'hiseq_rn' # 

        # aligner
        if self.aligner is None:
            self.aligner = 'bowtie2'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, 
            self.parallel_jobs)

        self.outdir = file_abspath(self.outdir)


    def init_fq(self):
        """
        Support fastq files
        fq1 (required)
        fq2 (optional)
        """
        # rep_list
        if isinstance(self.rep_list, list):
            pass

        elif isinstance(self.fq1, list):
            self.fq1 = file_abspath(self.fq1)

            # file exists
            if not file_exists(self.fq1):
                raise ValueError('--fq1, file not exists: {}'.format(self.fq1))

            # rep list
            q_names = [fq_name(i, pe_fix=True) for i in self.fq1]
            self.rep_list = [self.outdir + '/' + i for i in q_names]

            # fq2
            if isinstance(self.fq2, list):
                self.fq2 = file_abspath(self.fq2)

                if not file_exists(self.fq2):
                    raise ValueError('--fq2, file not exists: {}'.format(
                        self.fq2))

                # paired
                if not fq_paired(self.fq1, self.fq2):
                    raise ValueError('--fq1, --fq2, file not paired')
        else:
            raise ValueError('rep_list, fq1,fq2 required')

        # smp_name
        self.smp_name = fq_name_rmrep(self.rep_list).pop()


    def init_index(self):
        """
        alignment index
        genome, extra_index, genome_index

        output: genome_index, spikein_index
        """
        # check genome size
        if isinstance(self.extra_index, str):
            if self.genome_size < 1:
                ai = AlignIndex(index=self.extra_index)
                self.genome_size = ai.index_size()
                self.gsize_file = ai.index_size(return_file=True)
                self.genome_index = self.extra_index
        elif isinstance(self.genome, str):
            self.gsize_file = Genome(genome=self.genome).get_fasize()
            gsize = self.gsize_file
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))
            self.genome_index = AlignIndex(aligner=self.aligner).search(
                genome=self.genome, group='genome')
        elif isinstance(self.genome_index, str):
            pass
        else:
            raise ValueError('index failed, extra_index, genome, genome_index')

        # check spikein
        if isinstance(self.spikein_index, str):
            ai = AlignIndex(index=self.spikein_index, aligner=self.aligner)
            if not ai.is_index():
                raise ValueError('spikein_index failed, {}'.format(
                    self.spikein_index))
        elif isinstance(self.spikein, str):
            self.spikein_index = AlignIndex(aligner=self.aligner).search(
                genome=self.spikein, group='genome')
        else:
            self.spikein_index = None


    def init_files(self, create_dirs=True):
        """
        Prepare directories, files
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')

        default_dirs = {
            'align_dir': 'align',
            'bam_dir': 'bam_files',
            'bw_dir': 'bw_files',
            'bg_dir': 'bg_files',
            'peak_dir': 'peak',
            'motif_dir': 'motif',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json',
            'bam_rmdup': self.bam_dir + '/' + self.smp_name + '.rmdup.bam',
            'bam_proper_pair': self.bam_dir + '/' + self.smp_name + '.proper_pair.bam',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bed': self.bam_dir + '/' + self.project_name + '.bed',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',            
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'align_scale_txt': self.align_dir + '/' + 'scale.txt',
            'align_flagstat': self.align_dir + '/' + self.smp_name + '.flagstat',
            'align_stat': self.align_dir + '/' + self.smp_name + '.bowtie2.stat',
            'align_json': self.align_dir + '/' + self.smp_name + '.bowtie2.json',

            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.pdf',
            'frip_txt': self.qc_dir + '/03.FRiP.txt',
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_npz': self.qc_dir + '/06.bam_cor.npz',
            'bam_cor_counts': self.qc_dir + '/06.bam_cor.counts.tab',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png',
            'peak_idr_png': self.qc_dir + '/07.peak_idr.png',
            'peak_idr_txt': self.qc_dir + '/07.peak_idr.txt',
            'peak_overlap_png': self.qc_dir + '/08.peak_overlap.png',
            'peak_overlap_tiff': self.qc_dir + '/08.peak_overlap.tiff',
            'bam_fingerprint': self.qc_dir + '/09.fingerprint.png'
        }
        self = update_obj(self, default_files, force=True) # key

        if create_dirs:
            check_path([
                self.config_dir,
                self.align_dir,
                self.bam_dir, 
                self.bw_dir, 
                self.bg_dir,
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir])


class CnR1Config(object):
    """
    Prepare directories for R1 single replicates

    require: fq1/fq2, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_fq()
        self.init_index()
        self.init_files()


    def init_args(self):
        """
        required arguments for CnR analysis
        """
        args_init = {
            'is_ip': True,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': None,
            'smp_name': None,
            'aligner': 'bowtie2',
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'extra_index': None,
            'threads': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'gsize_file': None,
            'gene_bed': None,
            'keep_tmp': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'hiseq_r1' # 

        # aligner
        if self.aligner is None:
            self.aligner = 'bowtie2'

        # smp_name
        if not isinstance(self.smp_name, str):
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        self.outdir = file_abspath(self.outdir)


    def init_fq(self):
        """
        Support fastq files
        fq1 (required)
        fq2 (optional)
        """
        # fq1
        if not isinstance(self.fq1, str):
            raise ValueError('--fq1, str expected, got {}'.format(
                type(self.fq1).__name__))

        # abs
        self.fq1 = file_abspath(self.fq1)

        # file exists
        if not file_exists(self.fq1):
            raise ValueError('--fq1, file not exists: {}'.format(self.fq1))

        # fq2
        if not self.fq2 is None:
            if not isinstance(self.fq2, str):
                raise ValueError('--fq2, None or str expected, got {}'.format(
                    type(self.fq2).__name__))

            # abs
            self.fq2 = file_abspath(self.fq2)

            if not file_exists(self.fq2):
                raise ValueError('--fq2, file not exists: {}'.format(self.fq2))

            # paired
            if not fq_paired(self.fq1, self.fq2):
                raise ValueError('--fq1, --fq2, file not paired')


    def init_index(self):
        """
        alignment index
        genome, extra_index, genome_index

        output: genome_index, spikein_index
        """
        # check genome size
        if isinstance(self.extra_index, str):
            if self.genome_size < 1:
                ai = AlignIndex(index=self.extra_index)
                self.genome_size = ai.index_size()
                self.gsize_file = ai.index_size(return_file=True)
                self.genome_index = self.extra_index
        elif isinstance(self.genome, str):
            self.gsize_file = Genome(genome=self.genome).get_fasize()
            gsize = self.gsize_file
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))
            self.genome_index = AlignIndex(aligner=self.aligner).search(
                genome=self.genome, group='genome')
        elif isinstance(self.genome_index, str):
            pass
        else:
            raise ValueError('index failed, extra_index, genome, genome_index')

        # check spikein
        if isinstance(self.spikein_index, str):
            ai = AlignIndex(index=self.spikein_index, aligner=self.aligner)
            if not ai.is_index():
                raise ValueError('spikein_index failed, {}'.format(
                    self.spikein_index))
        elif isinstance(self.spikein, str):
            self.spikein_index = AlignIndex(aligner=self.aligner).search(
                genome=self.spikein, group='genome')
        else:
            self.spikein_index = None


    def init_files(self, create_dirs=True):
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
            'bam_dir': 'bam_files',
            'bg_dir': 'bg_files',
            'bw_dir': 'bw_files',
            'peak_dir': 'peak',
            'motif_dir': 'motif',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json',
            'bam_rmdup': self.bam_dir + '/' + self.smp_name + '.rmdup.bam',
            'bam_proper_pair': self.bam_dir + '/' + self.smp_name + '.proper_pair.bam',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bed': self.bam_dir + '/' + self.project_name + '.bed',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',

            'trim_stat_txt': self.clean_dir + '/' + self.project_name + '.qc.stat',
            'align_scale_txt': self.align_dir + '/' + 'scale.txt',
            'align_flagstat': self.align_dir + '/' + self.smp_name + '.flagstat',
            'align_stat': self.align_dir + '/' + self.smp_name + '.bowtie2.stat',
            'align_json': self.align_dir + '/' + self.smp_name + '.bowtie2.json',
            'spikein_bam': self.spikein_dir + '/' + 'spikein.bam',
            'spikein_scale_txt': self.spikein_dir + '/' + 'scale.txt',
            'spikein_flagstat': self.spikein_dir + '/' + 'spikein.flagstat',
            'spikein_stat': self.spikein_dir + '/' + 'spikein.align.stat',
            'spikein_json': self.spikein_dir + '/' + 'spikein.align.json',

            'trim_summary_json': self.qc_dir + '/00.trim_summary.json',
            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.pdf',
            'frip_txt': self.qc_dir + '/03.FRiP.txt',
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png',
            'peak_idr_png': self.qc_dir + '/07.peak_idr.png',
            'peak_overlap_png': self.qc_dir + '/08.peak_overlap.png',
            'bam_fingerprint': self.qc_dir + '/09.fingerprint.png'
        }
        self = update_obj(self, default_files, force=True) # key

        # raw data
        self.raw_fq_list = [self.raw_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            self.raw_fq_list.append(self.raw_dir + '/' + \
                os.path.basename(self.fq2))

        ## clean data
        self.clean_fq_list = [self.clean_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            self.clean_fq_list.append(self.clean_dir + '/' + \
                os.path.basename(self.fq2))

        if create_dirs:
            check_path([
                self.config_dir, 
                self.raw_dir,
                self.clean_dir,
                self.align_dir,
                self.spikein_dir,
                self.bam_dir, 
                self.bg_dir,
                self.bw_dir, 
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir])


class CnRtConfig(object):
    pass


class CnR(object):
    """
    Main port for CnR analysis

    if ip,ip-fq2, input,input-fq2: 
        (CnRn + CnRn) -> CnRx

    if ip_dir, input_dir:
        CnRx

    if fq1, f2:
        CnRn or CnR1
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = CnRConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def run_CnR_design(self):
        """
        Create design
        """
        CnRDesign(**self.__dict__).save()


    def run_CnR_r1(self):
        """
        Run single
        """
        CnR1(**self.__dict__).run()


    def run_CnR_rn(self):
        """
        Run multiple samples, rep_list
        """
        CnRRn(**self.__dict__).run()



    def run_CnR_from_design(self):
        """
        Run ChIPseq, multiple group

        multiple projects
        """
        design = Json(self.design).reader()

        for k, v in design.items():
            project_name = k
            args_local = v
            self = update_obj(self, args_local, force=True)
            CnRx(**self.__dict__).run()


    def run_CnR_rx(self):
        """
        main:
        Run whole pipeline

        required args:
        input_fq, ip_fq
        """
        CnRx(**self.__dict__).run()

        # # for pipeline
        # rx_args = self.__dict__
        # rx_local = {
        #     'ip_dir': ip.project_dir,
        #     'input_dir': input.project_dir,
        #     'fq1': None,
        #     'fq2': None,
        #     'ip': None,
        #     'ip_fq2': None,
        #     'input': None,
        #     'input_fq2': None,
        #     'spikein': self.spikein,
        #     'spikein_index': self.spikein_index,
        #     'extra_index': self.extra_index,
        #     'gene_bed': self.gene_bed
        # }
        # rx_args.update(rx_local)
        # rx = CnRx(**rx_args)
        # rx.run()


    def run(self):
        """
        Run all
        """
        if self.hiseq_type == 'build_design':
            self.run_CnR_design()
        elif self.hiseq_type == 'hiseq_from_design':
            self.run_CnR_from_design()
        elif self.hiseq_type == 'hiseq_r1':
            self.run_CnR_r1()
        elif self.hiseq_type == 'hiseq_rn':
            self.run_CnR_rn()
        elif self.hiseq_type == 'hiseq_rx':
            self.run_CnR_rx()
        else:
            raise ValueError('unknown hiseqtype: {}'.format(self.hiseq_type))


class CnRx(object):
    """
    Run ChIPseq for ip and input
    input: ip_dir, input_dir
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = CnRxConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def run_CnRn(self):
        """
        Run CnR for ip/input fastq files
        or copy file from CnRn
        """
        c_ip = CnRReader(self.ip_dir).is_hiseq_rn
        c_input = CnRReader(self.input_dir).is_hiseq_rn

        if c_ip and c_input:
            log.info('run CnRx for ip_dir and input_dir')
        elif isinstance(self.ip, list) and isinstance(self.input, list):
            # remove args
            # rep_list, ip_dir, input_dir, ip, ip_fq2, input, input_fq2
            args_local = self.__dict__.copy()
            rm_args = ['rep_list', 'ip_dir', 'input_dir', 'ip', 'ip_fq2',
                'input', 'input_fq2']
            [args_local.pop(i, None) for i in rm_args] # remove

            # for ip
            log.info('run CnRn for ip')                
            # for ip fastq files
            ip_args = args_local.copy()
            ip_local = {
                'build_design': False,
                'is_ip': True,
                'fq1': self.ip,
                'fq2': self.ip_fq2,
                'spikein': self.spikein,
                'spikein_index': self.spikein_index,
                'extra_index': self.extra_index,
                'genome_size': self.genome_size,
                'gene_bed': self.gene_bed,
                'hiseq_type': 'hiseq_rn'
                }
            ip_args.update(ip_local)
            CnRn(**ip_args).run()

            # for input fastq
            input_args = args_local.copy()
            input_local = {
                'build_design': False,
                'is_ip': False,
                'fq1': self.input,
                'fq2': self.input_fq2,
                'spikein': self.spikein,
                'spikein_index': self.spikein_index,
                'extra_index': self.extra_index,
                'genome_size': self.genome_size,
                'gene_bed': self.gene_bed,
                'hiseq_type': 'hiseq_rn'
                }
            input_args.update(input_local)
            CnRn(**input_args).run()
        else:
            log.error('CnRx() failed, unknown args')

        ## update 
        self.ip_args = CnRReader(self.ip_dir).args
        self.input_args = CnRReader(self.input_dir).args
        self.ip_bam_from = self.ip_args.get('bam', None)
        self.ip_bw_from = self.ip_args.get('bw', None)
        self.input_bam_from = self.input_args.get('bam', None)
        self.input_bw_from = self.input_args.get('bw', None)


    def copy_bam_files(self):
        """
        Copy bam files
        """
        file_symlink(self.ip_bam_from, self.ip_bam)
        file_symlink(self.input_bam_from, self.input_bam)
        # file_copy(self.ip_bam_from, self.ip_bam)
        # file_copy(self.input_from, self.input_bam)


    def copy_bw_files(self):
        """
        Copy bw files
        """
        file_symlink(self.ip_bw_from, self.ip_bw)
        file_symlink(self.input_bw_from, self.input_bw)
        # file_copy(self.ip_bw_from, self.ip_bw)
        # file_copy(self.input_bw_from, self.input_bw)


    def get_ip_over_input_bw(self):
        """
        Copy bw files
        """
        # ip over input, subtract
        try:
            bwCompare(self.ip_bw, self.input_bw, self.ip_over_input_bw, 'subtract',
                threads=self.threads, binsize=10)
        except:
            log.error('get_ip_over_input_bw() failed, see {}'.format(
                self.bw_dir))


    def call_peak2(self):
        """
        Call peaks using MACS2
        ip, input
        ...
        """
        args_global = self.__dict__.copy()
        args_required = ['genome_size', 'gsize_file']
        args_local = dict((k, args_global[k]) for k in args_required if
            k in args_global)

        peak = Macs2(
            ip=self.ip_args.get('bam', None),
            control=self.input_args.get('bam', None),
            genome=self.genome,
            output=getattr(self, 'peak_dir', None),
            prefix=self.ip_args.get('smp_name', None),            
            # genome_size=self.genome_size,
            # gsize_file=self.gsize_file,
            **args_local)

        ## call peaks
        peak.callpeak()
        # peak.bdgcmp(opt='ppois')
        # peak.bdgcmp(opt='FE')
        # peak.bdgcmp(opt='logLR')

        # ## annotation
        # peak.broadpeak_annotation()


    def call_peak(self):
        """
        Call peaks using MACS2/SEACR
        """
        args_local = {
            'ip': self.ip_bam,
            'input': self.input_bam,
            'outdir': self.peak_dir,
            'genome': self.genome
        }
        CallPeak(**args_local).run()
        CallPeak(method='macs2', **args_local).run()
        CallPeak(method='seacr', **args_local).run()


    def qc_tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        # bw_list = [self.ip_bw, self.input_bw, self.ip_over_input_bw]
        bw_list = [self.ip_bw, self.input_bw]
        bw_list_arg = ' '.join(bw_list)

        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'reference-point',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(bw_list_arg),
            '-o {}'.format(self.tss_enrich_matrix),
            '--referencePoint TSS',
            '-b 2000 -a 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.tss_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.tss_enrich_matrix),
            '-o {}'.format(self.tss_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.tss_enrich_png) and not self.overwrite:
            log.info('qc_tss_enrich() skipped, file exists: {}'.format(
                self.tss_enrich_png))
        else:            
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.tss_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_tss_enrich() failed, see: ')


    def qc_genebody_enrich(self):
        """
        Calculate the TSS enrichment
        """
        # bw_list = [self.ip_bw, self.input_bw, self.ip_over_input_bw]
        bw_list = [self.ip_bw, self.input_bw]
        bw_list_arg = ' '.join(bw_list)

        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'scale-regions',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(bw_list_arg),
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


    def qc_bam_fingerprint(self):
        """
        Calculate fingerprint for bam files
        """
        Bam2fingerprint(
            bam_list=[self.ip_bam, self.input_bam],
            prefix='09.fingerprint',
            outdir=self.qc_dir,
            threads=self.threads).run()


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        hiseq_report_html = os.path.join(
            self.report_dir, 
            'HiSeq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        if file_exists(hiseq_report_html):
            log.info('report() skipped, file exists: {}'.format(
                hiseq_report_html))
        else:
            run_shell_cmd(cmd) 


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # run
        self.run_CnRn()
        self.copy_bam_files()
        self.copy_bw_files()
        self.get_ip_over_input_bw()
        self.call_peak()
        # qc
        self.qc_tss_enrich()
        self.qc_genebody_enrich()
        self.qc_bam_fingerprint()
        self.report()


class CnRn(object):
    """
    Run CnR for multiple replicates, merge replicates
    input: rep_list/fq1
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = CnRnConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def get_bam_list(self):
        """
        get the proper_mapped.bam files
        """
        return [CnRReader(i).args.get('bam', None) for i in self.rep_list]


    def get_bw_list(self):
        """
        get the proper_mapped.bam files
        """
        return [CnRReader(i).args.get('bw', None) for i in self.rep_list]


    def get_peak_list(self):
        """
        get the .narrowPeak files
        """
        return [CnRReader(i).args.get('peak', None) for i in self.rep_list]


    def merge_bam(self):
        """
        Merge replicates, BAM
        """
        self.bam_list = self.get_bam_list()

        cmd = ' '.join([
            'samtools merge -',
            ' '.join(self.bam_list),
            '| samtools sort -o {} -'.format(self.bam),
            '&& samtools index {}'.format(self.bam)])

        if os.path.exists(self.bam):
            log.info('merge_bam() skipped, file exists: {}'.format(
                self.bam))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('merge_bam() failed.')

        ## calculate norm scale
        self.align_scale = self.cal_norm_scale(self.bam)
        with open(self.align_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.align_scale))


    def get_bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        # remove dup
        if rmdup:
            if file_exists(self.bam):
                if file_exists(self.bam_rmdup) and not self.overwrite:
                    log.info('rmdup() skipped, file exists: {}'.format(
                        self.bam_rmdup))
                else:
                    Bam(self.bam).rmdup(self.bam_rmdup)
                    Bam(self.bam_rmdup).index()


    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam,
            'scaleFactor': self.align_scale,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 0,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }

        Bam2bw(**args_local).run()


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
            cmd_txt = self.bg_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
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
            '{} {} {}'.format(self.bg, self.gsize_file, self.bw)
            ])

        if file_exists(self.bw) and not self.overwrite:
            log.info('bg_to_bw() skipped, file exists:{}'.format(self.bw))
        else:
            cmd_txt = self.bw_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bg_to_bw() failed')


    def call_peak2(self):
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
        genome_size = getattr(self, 'genome_size', 0)
        gsize_file = getattr(self, 'gsize_file', None)

        if check_file(self.peak):
            log.info('call_peak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, genome, output, prefix, atac=False,
                genome_size=genome_size, gsize_file=gsize_file).callpeak()


    def call_peak(self):
        """
        Call peaks using MACS2/SEACR
        """
        args_local = {
            'ip': self.bam,
            'input': None,
            'outdir': self.peak_dir,
            'genome': self.genome
        }
        CallPeak(method='macs2', **args_local).run()
        CallPeak(method='seacr', **args_local).run()


    def cal_norm_scale(self, bam, norm=1000000):
        """
        scale factor
        
        Bam().count_reads()
        """
        bam = Bam(bam)
        is_pe = bam.isPaired()
        n_map = bam.getNumberOfAlignments()
        if is_pe:
            n_map = n_map/2

        if n_map > 0:
            n_scale = norm/n_map
        else:
            log.error('no mapped reads detected')
            n_scale = 1

        return n_scale


    def qc_bam_cor(self, window=500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp
        
        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
            --outRawCounts *counts.tab -b bam
        """
        args_local = {
            'bam': self.get_bam_list(),
            'outdir': self.qc_dir,
            'prefix': '06.bam_cor',
            'threads': self.threads,
            'overwrite': self.overwrite,
            'binsize': self.binsize
        }
        Bam2cor(**args_local).run()


    def qc_peak_overlap(self):
        """
        Compute the overlaps between overlaps
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_dir,
            'prefix': '08.peak_overlap',
            'overwrite': self.overwrite
        }   

        BedOverlap(**args_local).run()


    def qc_peak_idr(self):
        """
        Calculate IDR for replicates
        1 vs 1
        peak files
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_dir,
            'prefix': '07.peak_idr',
            'overwrite': self.overwrite
        }

        PeakIDR(**args_local).run()


    def get_peak_frip(self):
        """
        Save all FRiP.txt file to one
        """
        # get list
        self.frip_list = [CnRReader(i).args.get('frip_txt', None) for \
            i in self.rep_list]
        self.frip_list = [i for i in self.frip_list if file_exists(i)]

        if len(self.frip_list) > 0:
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


    def qc_align_txt(self):
        """
        Copy align.txt files to align/
        """
        pass


    def qc_tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        bw_list = self.get_bw_list()
        bw_list.append(self.bw)
        bw_list_arg = ' '.join(bw_list)

        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'reference-point',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(bw_list_arg),
            '-o {}'.format(self.tss_enrich_matrix),
            '--referencePoint TSS',
            '-b 2000 -a 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.tss_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.tss_enrich_matrix),
            '-o {}'.format(self.tss_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.tss_enrich_png) and not self.overwrite:
            log.info('qc_tss_enrich() skipped, file exists: {}'.format(
                self.tss_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.tss_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_tss_enrich() failed, see: {}'.format(
                        self.tss_enrich_matrix_log))


    def qc_genebody_enrich(self):
        """
        Calculate the TSS enrichment
        """
        bw_list = self.get_bw_list()
        bw_list.append(self.bw)
        bw_list_arg = ' '.join(bw_list)

        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'scale-regions',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(bw_list_arg),
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
                    log.error('qc_genebody_enrich() failed, see: {}'.format(
                        self.genebody_enrich_matrix_log))


    def qc_bam_fingerprint(self):
        """
        Calculate fingerprint for bam files
        """
        Bam2fingerprint(
            bam_list=self.get_bam_list(),
            prefix='09.fingerprint',
            outdir=self.qc_dir,
            threads=self.threads).run()


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        hiseq_report_html = os.path.join(
            self.report_dir, 
            'HiSeq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        if file_exists(hiseq_report_html):
            log.info('report() skipped, file exists: {}'.format(
                hiseq_report_html))
        else:
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
        args_required = ['aligner', 'fq1', 'fq2', 'genome', 'gene_bed',
            'genome_size', 'is_ip', 'trimmed', 'outdir', 'overwrite', 
            'parallel_jobs', 'threads', 'spikein', 'spikein_index',
            'extra_index']
        args_local = dict((k, args_tmp[k]) for k in args_required 
            if k in args_tmp)

        args_input = self.pick_fq_samples(i)
        args_init = {
            'build_design': None,
            'design': None,
            'gene_bed': self.gene_bed
        }
        args_local.update(args_input) # update fq1/rep_list/group
        args_local.update(args_init) #

        CnR1(**args_local).run()


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
            # run
            self.merge_bam()
            self.get_bam_rmdup()
            # self.bam_to_bw()
            self.bam_to_bg()
            self.bg_to_bw()
            self.call_peak()
            # qc
            self.qc_bam_cor()
            self.qc_peak_overlap()
            self.qc_peak_idr()
            self.get_peak_frip()
            self.qc_frip()
            self.qc_align_txt()
            self.qc_tss_enrich()
            self.qc_genebody_enrich()
            self.qc_bam_fingerprint()
            self.report()
        elif len(self.rep_list) == 1:
            # new files, symlinks
            log.warning('merge() skipped, Only 1 replicate detected')
            # copy files: bam, bw, peak
            rep = CnRReader(self.rep_list[0]).args
            bam = rep.get('bam', None)
            bg = rep.get('bg', None)
            bw = rep.get('bw', None)
            peak = rep.get('peak', None)
            peak_seacr = rep.get('seacr', None)

            align_scale_txt = rep.get('align_scale_txt', None)
            align_json = rep.get('align_json')

            file_symlink(align_scale_txt, self.align_scale_txt)
            file_symlink(align_json, self.align_json)
            file_symlink(bam, self.bam)
            Bam(self.bam).index()
            file_symlink(bg, self.bg)
            file_symlink(bw, self.bw)
            file_symlink(peak, self.peak)
            file_symlink(peak_seacr, self.peak_seacr)

            self.get_peak_frip()
            self.qc_frip()
            self.qc_tss_enrich()
            self.qc_genebody_enrich()
            self.qc_bam_fingerprint()
            self.report()

        else:
            log.error('merge() failed, no rep detected')
            raise ValueError('merge() failed, no rep detected')


class CnR1(object):
    """
    Run CUT&RUN basic for single fastq file (ip/input)

    align - bam - rmdup - proper_paired - peak - bw - qc - report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = CnR1Config(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    # main pipeline #
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
            file_symlink(self.fq1, raw_fq1, absolute_path=True)
            file_symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=True):
        """
        using bowtie2, --local
        no need to trim adapters

        Trim reads:
        Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        args_trim = self.__dict__.copy()
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
                file_symlink(fq1, clean_fq1)
            else:
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
            ## fq2
            if not fq2 is None:
                if is_gz(fq2):
                    file_symlink(fq2, clean_fq2)
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
                file_symlink(fq1, clean_fq1)
                file_symlink(fq2, clean_fq2)


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
            'genome': self.genome,
            'extra_index': self.extra_index,
            'aligner': self.aligner,
            'keep_tmp': self.keep_tmp
        }
        args_local.update(args_init)

        # output
        bam = args_local.get('bam')
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            Align(**args_local).run()

        # copy files
        g = CnRReader(self.align_dir).args
        g_bam = g.get('bam', None)
        g_align_bam = g.get('bam', None)
        g_align_stat = g.get('align_stat', None)
        g_align_json = g.get('align_json', None)
        g_align_flagstat = g.get('align_flagstat', None)

        # copy files
        file_symlink(g_align_bam, self.bam)
        file_copy(g_align_stat, self.align_stat)
        file_copy(g_align_json, self.align_json)
        file_copy(g_align_flagstat, self.align_flagstat)

        ## calculate norm scale
        self.align_scale = self.cal_norm_scale(self.bam)
        with open(self.align_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.align_scale))


    def align_spikein(self):
        """
        Alignment PE reads, spikein
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_sp = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.spikein_dir,
            'genome': None,
            'index': self.spikein_index,
            'threads': self.threads,
            'keep_tmp': self.keep_tmp
        }

        # output
        bam = args_local.get('spikein_bam')
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            # print('!AAAA-1')
            # print_dict(args_sp)
            # sys.exit()
            Align(**args_sp).run()

        # copy files
        sp = CnRReader(self.spikein_dir).args
        sp_bam = sp.get('bam', None)
        sp_align_bam = sp.get('bam', None)
        sp_align_stat = sp.get('align_stat', None)
        sp_align_json = sp.get('align_json', None)
        sp_align_flagstat = sp.get('align_flagstat', None)

        # copy files
        file_copy(sp_align_bam, self.spikein_bam)
        file_copy(sp_align_stat, self.spikein_stat)
        file_copy(sp_align_json, self.spikein_json)
        file_copy(sp_align_flagstat, self.spikein_flagstat)

        ## calculate norm scale
        self.spikein_scale = self.cal_norm_scale(self.spikein_bam, 10000)
        with open(self.spikein_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.spikein_scale))


    def get_bam(self, dir):
        """
        Get the align bam file
        """
        bam_list = listfile(dir, '*.bam', recursive=True)
        bam_list = sorted(bam_list)
        bam_list = [i for i in bam_list if not i.endswith('.raw.bam')] #

        return bam_list[-1]


    def get_bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        bam_raw = self.get_bam(self.align_dir)

        # symlink bam
        file_symlink(bam_raw, self.bam)

        # remove dup
        if rmdup:
            Bam(self.bam).rmdup(self.bam_rmdup)
            Bam(self.bam_rmdup).index()
        else:
            file_symlink(self.bam, self.bam_rmdup)


    def call_peak(self):
        """
        Call peaks using MACS2/SEACR
        """
        args_local = {
            'ip': self.bam,
            'input': None,
            'outdir': self.peak_dir,
            'genome': self.genome
        }
        CallPeak(method='macs2', **args_local).run()
        CallPeak(method='seacr', **args_local).run()


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
            cmd_txt = self.bg_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
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
            '{} {} {}'.format(self.bg, self.gsize_file, self.bw)
            ])

        if file_exists(self.bw) and not self.overwrite:
            log.info('bg_to_bw() skipped, file exists:{}'.format(self.bw))
        else:
            cmd_txt = self.bw_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bg_to_bw() failed')


    def pipe_genome(self):
        """
        Run for genome
        """
        self.prep_raw()
        self.trim()
        self.align_genome()
        self.get_bam_rmdup()
        self.call_peak()
        self.bam_to_bg()
        self.bg_to_bw()


    def pipe_spikein(self):
        """
        Run for spikein
        """
        if isinstance(self.spikein_index, str):
            self.align_spikein() # run alignment


    # quality control #
    def cal_norm_scale(self, bam, norm=1000000):
        """
        scale factor
        
        Bam().count_reads()
        """
        bam = Bam(bam)
        is_pe = bam.isPaired()
        n_map = bam.getNumberOfAlignments()
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

        #sample input   output  percent
        fq_rep1      2234501 2234276 99.99%

        """
        if file_exists(self.trim_stat_txt):
            with open(self.trim_stat_txt, 'rt') as r:
                lines = r.readlines()
            lines = [i for i in lines if not i.startswith('#')]
            line = lines.pop()
            fqname, n_input, n_output, n_pct = line.strip().split('\t')
            n_pct = eval(n_pct.strip('%'))

            d = {
                'id': fqname,
                'input': n_input,
                'output': n_output,
                'out_pct': n_pct,
                'rm_pct': 100 - n_pct
            }
        else:
            d = {
                'id': self.smp_name,
                'input': 0,
                'output': 0,
                'out_pct': 100,
                'rm_pct': 0
            }
        Json(d).writer(self.trim_summary_json)


    def qc_align_summary(self):
        """
        Save summary to file
        """
        # align to genome
        align = Json(self.align_json).reader()

        # duplicates
        nodup = Bam(self.bam_rmdup)
        is_pe = nodup.isPaired()
        n_nodup = nodup.getNumberOfAlignments()
        if is_pe:
            n_nodup = n_nodup/2
        
        # spikein 
        if file_exists(self.spikein_bam):
            n_spikein = Bam(self.spikein_bam).getNumberOfAlignments()
            if is_pe:
                n_spikein = n_spikein/2
        else:
            n_spikein = 0

        # Json(align).writer()
        align['nodup'] = n_nodup
        align['spikein'] = n_spikein
        align['chrM'] = self.qc_mito()

        # save to file
        Json(align).writer(self.align_summary_json)


    def qc_lendist(self):
        """
        Create length distribution, txt
        """
        if os.path.exists(self.lendist_txt) and self.overwrite is False:
            log.info('lendist() skipped: file exists: {}'.format(
                self.lendist_txt))
        else:
            BamPEFragSize(self.bam).saveas(self.lendist_txt)


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            frip, n, total = peak_FRiP(self.peak, self.bam, threads=1)

            hd = ['FRiP', "peak_reads", "total_reads", "id"]
            n = list(map(str, [frip, n, total]))
            n.append(self.project_name)
            with open(self.frip_txt, 'wt') as w:
                w.write('\t'.join(hd) + '\n')
                w.write('\t'.join(n) + '\n')


    def qc_mito(self):
        """
        Mito reads, percentage
        """
        lines = pysam.idxstats(self.bam).splitlines()
        # extract chrM, MT
        lines = [i for i in lines if re.search('^(chrM|MT)', i)]

        try:
            nreads = sum(
                map(int, [x.split("\t")[2]
                          for x in lines if not x.startswith("#")]))

        except IndexError as msg:
            raise IndexError(
                "can't get number of reads from bamfile, msg=%s, data=%s" %
                (msg, lines))

        return nreads


    def qc_tss(self):
        """
        TSS enrichment

        computeMatrix referencepoint -b -R gene.bed -S in.bw -o mat.gz

        require:
        plotHeatmap()
        plotProfile()
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'reference-point',
            '--referencePoint TSS',
            '-b 2000 -a 2000 --binSize 10 --sortRegions descend --skipZeros',
            '--samplesLabel {}'.format(self.smp_name),
            '-p {}'.format(self.threads),
            '--regionsFileName {}'.format(self.gene_bed),
            '--scoreFileName {}'.format(self.bw),
            '-o {}'.format(self.tss_enrich_matrix),
            '2> {}'.format(self.tss_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.tss_enrich_matrix),
            '-o {}'.format(self.tss_enrich_png),
            '--dpi 300'
            ])

        if file_exists(self.tss_enrich_png) and not self.overwrite:
            log.info('qc_tss() skipped, file exists')
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                cmd_txt = self.qc_dir + '/tss_enrich.sh'
                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('failed to run computeMatrix, see: {}'.format(
                        self.tss_enrich_matrix_log ))


    def qc_tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'reference-point',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(self.bw),
            '-o {}'.format(self.tss_enrich_matrix),
            '--referencePoint TSS',
            '-b 2000 -a 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.tss_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.tss_enrich_matrix),
            '-o {}'.format(self.tss_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.tss_enrich_png) and not self.overwrite:
            log.info('qc_tss_enrich() skipped, file exists: {}'.format(
                self.tss_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.tss_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_tss_enrich() failed, see: ')


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


    def qc_bam_fingerprint(self):
        """
        Calculate fingerprint for bam files
        """
        Bam2fingerprint(
            bam_list=self.get_bam_list(),
            prefix='09.fingerprint',
            outdir=self.qc_dir,
            threads=self.threads).run()


    def pipe_qc(self):
        """
        Quality control
        """
        self.qc_trim_summary()
        self.qc_align_summary()
        self.qc_lendist()
        self.qc_mito()
        self.qc_tss_enrich()
        self.qc_genebody_enrich()


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        hiseq_report_html = os.path.join(
            self.report_dir, 
            'HiSeq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        if file_exists(hiseq_report_html):
            log.info('report() skipped, file exists: {}'.format(
                hiseq_report_html))
        else:
            run_shell_cmd(cmd)


    def run(self):
        """
        Run all steps for ChIPseq pipeline
        """
        # init dir
        args = self.__dict__.copy()

        self.pipe_genome()
        self.pipe_spikein()
        self.pipe_qc()
        self.report()

        # remove clean fastq files
        del_list = self.clean_fq_list
        if not self.keep_tmp:
            file_remove(del_list, ask=False)


class CnRt(object):
    pass


class CnRReader(object):
    """
    Read config.pickle from the local directory
    """
    def __init__(self, x):
        self.x = x
        self.read()
        self.hiseq_type = self.args.get('hiseq_type', None)

        self.is_hiseq_r1 = self.hiseq_type == 'hiseq_r1'
        self.is_hiseq_rn = self.hiseq_type == 'hiseq_rn'
        self.is_hiseq_rx = self.hiseq_type == 'hiseq_rx'
        self.is_hiseq_rt = self.hiseq_type == 'hiseq_rt'
        self.is_hiseq_rd = self.hiseq_type == 'build_design'


        if isinstance(self.hiseq_type, str):
            self.is_hiseq = self.hiseq_type.startswith('hiseq_')
        else:
            self.is_hiseq = False


    def read(self):
        """
        # hiseq
        hiseq
          |-config
          |   |-config.pickle

        # alignment
        align_dir
          |- smp_nmae
          |    |- index
          |    |    |- config.pickle
        """
        p1x = os.path.join(self.x, 'config', 'config.pickle')
        p2x = os.path.join(self.x, '*', '*', 'config.pickle')
        p1 = glob.glob(p1x)
        p2 = glob.glob(p2x)

        # read config
        if len(p1) == 1:
            self.args = pickle2dict(p1[0])
        elif len(p2) == 1:
            self.args = pickle2dict(p2[0])
        else:
            self.args = {}


class CnRDesign(object):
    """
    Prepare ip/input samples for CnR analysis
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
            'trimmed': True,
            'smp_name': None,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': str(pathlib.Path.cwd()),
            'aligner': None,
            'align_to_chrM': True,
            'extra_index': None,
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

        # check genome size
        if isinstance(self.extra_index, list):
            if self.genome_size < 1:
                ai = AlignIndex(index=self.extra_index[-1])
                self.genome_size = ai.index_size()
                self.gsize_file = ai.index_size(return_file=True)
        elif isinstance(self.genome, str):
            self.gsize_file = Genome(genome=self.genome).get_fasize()
            gsize = self.gsize_file
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))
        else:
            raise ValueError('unknown genome, extra_index')

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
            'aligner': None,
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'extra_index': None
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

        # check genome size
        if isinstance(self.extra_index, list):
            if self.genome_size < 1:
                ai = AlignIndex(index=self.extra_index[-1])
                self.genome_size = ai.index_size()
                self.gsize_file = ai.index_size(return_file=True)
        elif isinstance(self.genome, str):
            gsize = Genome(genome=self.genome).get_fasize()
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))
        else:
            raise ValueError('unknown genome, extra_index')

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
            'aligner': None,
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'extra_index': None
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

        # check genome size
        if isinstance(self.extra_index, list):
            if self.genome_size < 1:
                ai = AlignIndex(index=self.extra_index[-1])
                self.genome_size = ai.index_size()
                self.gsize_file = ai.index_size(return_file=True)
        elif isinstance(self.genome, str):
            gsize = Genome(genome=self.genome).get_fasize()
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))
        else:
            raise ValueError('unknown genome, extra_index')

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
            'aligner': None,
            'align_to_chrM': True,
            'extra_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'trimmed': True
        }
        self = update_obj(self, args_init, force=False)

        # # update fq from design
        # if file_exists(self.design) and not self.build_design:
        #     self.design = file_abspath(self.design)
        #     args_design = Json(self.design).reader()
        #     self = update_obj(self, args_design, force=True)

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

        # check genome size
        if isinstance(self.extra_index, list):
            if self.genome_size < 1:
                ai = AlignIndex(index=self.extra_index[-1])
                self.genome_size = ai.index_size()
                self.gsize_file = ai.index_size(return_file=True)
        elif isinstance(self.genome, str):
            self.gsize_file = Genome(genome=self.genome).get_fasize()
            gsize = self.gsize_file # Genome(genome=self.genome).get_fasize()
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))
        elif not self.build_design:
            raise ValueError('unknown genome, extra_index')


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
                if not distance(q1, q2) == 1:
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
            self.chipseq_type = 'chipseq_rx'
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
            file_symlink(self.fq1, raw_fq1)
            file_symlink(self.fq2, raw_fq2)


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
                file_symlink(fq1, clean_fq1)
            else:
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
            ## fq2
            if not fq2 is None:
                if is_gz(fq2):
                    file_symlink(fq2, clean_fq2)
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
                file_symlink(fq1, clean_fq1)
                file_symlink(fq2, clean_fq2)


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
            'extra_index': self.extra_index,
            'extra_para': '--very-sensitive-local', # --no-mixed -X 800
            'aligner': None
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

        # determine the genome size
        genome_size = getattr(self, 'genome_size', 0)
        gsize_file = getattr(self, 'gsize_file', None)

        if check_file(self.peak):
            log.info('call_peak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, self.genome, self.peak_dir, self.project_name,
                atac=False, genome_size=genome_size, 
                gsize_file=gsize_file).callpeak()


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
        chipseq_report_html = os.path.join(
            self.report_dir, 
            'ChIPseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        if file_exists(chipseq_report_html):
            log.info('report() skipped, file exists: {}'.format(chipseq_report_html))
        else:
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
        self.trim(trimmed=trimmed)

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
        genome_size = getattr(self, 'genome_size', 0)
        gsize_file = getattr(self, 'gsize_file', None)

        if check_file(self.peak):
            log.info('call_peak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, genome, output, prefix, atac=False,
                genome_size=genome_size, gsize_file=gsize_file).callpeak()


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
        chipseq_report_html = os.path.join(
            self.report_dir, 
            'ChIPseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        if file_exists(chipseq_report_html):
            log.info('report() skipped, file exists: {}'.format(chipseq_report_html))
        else:
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
            'genome_size', 'is_ip', 'trimmed', 'outdir', 'overwrite', 
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
            file_symlink(rep_dict.get('bam', None), self.bam)
            file_symlink(rep_dict.get('bw', None), self.bw)
            file_symlink(rep_dict.get('peak', None), self.peak)
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
        file_symlink(self.ip_args.get('bam', None), self.ip_bam)
        file_symlink(self.input_args.get('bam', None), self.input_bam)


    def get_bw_files(self):
        """
        Copy bw files
        """
        file_symlink(self.ip_args.get('bw', None), self.ip_bw)
        file_symlink(self.input_args.get('bw', None), self.input_bw)

        # ip over input, subtract
        bwCompare(self.ip_bw, self.input_bw, self.ip_over_input_bw, 'subtract',
            threads=self.threads, binsize=10)


    def call_peak(self):
        """
        Call peaks using MACS2
        ip, input
        ...
        """
        args_global = self.__dict__.copy()
        args_required = ['genome_size', 'gsize_file']
        args_local = dict((k, args_global[k]) for k in args_required if
            k in args_global)

        peak = Macs2(
            ip=self.ip_args.get('bam', None),
            control=self.input_args.get('bam', None),
            genome=self.genome,
            output=getattr(self, 'peak_dir', None),
            prefix=self.ip_args.get('smp_name', None),            
            # genome_size=self.genome_size,
            # gsize_file=self.gsize_file,
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
        chipseq_report_html = os.path.join(
            self.report_dir, 
            'ChIPseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        if file_exists(chipseq_report_html):
            log.info('report() skipped, file exists: {}'.format(chipseq_report_html))
        else:
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
            'rep_list': None,
            'extra_index': self.extra_index,
            'genome_size': self.genome_size}
        ip_args.update(ip_local)
        ip = ChIPseqRn(**ip_args)

        # for input fastq
        input_args = self.__dict__
        input_local = {
            'is_ip': False,
            'fq1': self.input,
            'fq2': self.input_fq2,
            'build_design': None,
            'rep_list': None,
            'extra_index': self.extra_index,
            'genome_size': self.genome_size}
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
            'input_fq2': None,
            'extra_index': self.extra_index
        }
        rx_args.update(rx_local)
        rx = ChIPseqRx(**rx_args)
        rx.run()


    def run(self):
        """
        Run all
        """
        if self.chipseq_type == 'build_design':
            self.run_chipseq_design()
        elif self.chipseq_type == 'chipseq_rx_from_design':
            self.run_chipseq_rx_from_design()
        elif self.chipseq_type == 'chipseq_r1':
            ChIPseqR1(**self.__dict__).run()
        elif self.chipseq_type == 'chipseq_rn':
            ChIPseqRn(**self.__dict__).run()
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

