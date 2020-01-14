# -*- coding:utf-8 -*-

"""
###############
1. unmap files: 
  - bowtie  : _1.fastq, _2.fastq
  - bowtie2 : .1.fastq, .2.fastq

###############


Align fasta/q files to reference sequence

bowtie
bowtie2
STAR
bwa
hisat2
...


arguments:

fq
index
outdir

# optional

unique
mismatch
insertion (bowtie2 -X)

max-multi

# output
sam

unmap

# custom args
arg_str


## standard

input: fq(s), index, args
output: sam, unmap, log

## update: 2020-01-08
1. uniform code style: self ->  dict -> arguments.txt + pickle


"""

import os
import sys
import re
import pathlib
import shutil
import logging
import pandas as pd

from hiseq.utils.args import args_init
# from args import args_init
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions



def sam2bam(sam, bam, sort=True, extra_para=''):
    """
    Convert sam to bam, sort by position
    """
    samtools_exe = shutil.which('samtools')

    cmd = '{} view -Subh {} {} | samtools sort -o {} -'.format(
        samtools_exe,
        extra_para,
        sam,
        bam)

    run_shell_cmd(cmd)
    return bam


def pickFq(x, check_exists=True, is_str=True):
    """
    arguments:
    x, path to the fastq fils, str or list
    check_exists, boolen, 
    is_str, if not, choose the first one
    
    if x is None, warning
    """
    if isinstance(x, list):
        if is_str:
            log.warning('choose the first file: {}'.format(x))
            fq = x[0]
        else:
            fq = x
    elif isinstance(x, str):
        fq = x
    else:
        # log.warning('expect str, list, not {}'.format(type(x)))
        fq = None

    # exists
    if check_exists:
        if isinstance(fq, list):
            tag = all(map(os.path.exists, fq))
        elif isinstance(fq, str):
            tag = os.path.exists(fq)
        else:
            tag = False
    else:
        tag = True

    return fq if tag else None


## to-do
##   - add extra parameters for aligner: '-X 2000' for bowtie2, ...
## level-1
@Logger('INFO')
class Alignment(object):
    def __init__(self, **kwargs):
        """
        The top level of alignment function
        A port for actual usage
        
        args:
          - aligner
          - fq
          - outdir
          - fq2
          - index_list

        # required:
        fq1, outdir, genome
        """
        args = args_init(kwargs, align=True) # global

        if args['fq1'] is None:
            raise Exception('{:>10} : -i, --fq1 no input detected'.format('error'))

        if isinstance(args['fq1'], str):
            args['fq1'] = [args['fq1']]
        if isinstance(args['fq2'], str):
            args['fq2'] = [args['fq2']]

        self.args = args


    def run(self):
        return AlignNFastx(**self.args).run()


##------ Align core: BEGIN ------##

## level-1: 1 fastx, 1 index
## config for all alignment
## port-config
## !!! to-do
##  - aligner_supported()
## for 1 fq, 1 index
class AlignConfig(object):

    def __init__(self, outdir_fixed=False, search_index=False, **kwargs):
        """
        Config for alignment, for one fq, one index
        directories
        files
        ...
        """
        args = args_init(kwargs, align=True)
        # print('!BBBB', args)

        ## update args: aligner, fq, index, outdir, fq2
        args['outdir_fixed'] = outdir_fixed
        args['search_index'] = search_index
        args['current_dir'] = args.get('current_dir', str(pathlib.Path.cwd()))

        ## update spikein
        if args['spikein'] == args['genome']:
            args['spikein'] = None

        ## global: variables for Class
        self.aligner = args['aligner']
        self.fq1 = args['fq1']
        self.fq2 = args['fq2']
        self.index = args.get('index', None)
        self.outdir = args['outdir']
        self.outdir_fixed = outdir_fixed
        self.fx_type = self.format = Fastx(args['fq1']).format # fasta, fastq # pigz ?!
        self.args = args
        self.index_list = self.get_index_list()
        if not search_index:
            ## main configuration:
            tmp = self.config() # update, global variables

            ## update args:
            args['fqname'] = self.fqname
            args['out_prefix'] = self.out_prefix
            args['sam'] = self.sam
            args['bam'] = self.bam
            args['log'] = self.log
            args['unmap'] = self.unmap
            args['unmap1'] = self.unmap1
            args['unmap2'] = self.unmap2
            args['sub_outdir'] = self.sub_outdir
            args['index_list'] = self.index_list

            ## check step: fq, bam, overwrite
            self.check_status = self.check()


    def config(self):
        """
        Create directories, filenames
        *fq*     The fastq file (or fq1 of PE)
        *index*  The aligner index file
        *args*   Align_path, output directory of the alignment    

        return files:
        prefix, bam, bed, log, unmap
        """
        args = self.args.copy()

        ## check args
        # assert aligner_supported(aligner) # !!! to-do
        assert isinstance(self.aligner, str)
        assert isinstance(self.outdir, str) # the top-level output directory
        assert isinstance(self.outdir_fixed, bool)
        assert os.path.exists(self.fq1)
        assert isinstance(self.index, str)
        assert AlignIndex(self.aligner, self.index).is_index()
        
        fqname = file_prefix(self.fq1)[0]
        if not self.fq2 is None:
            fqname = re.sub('[._][rR]?1$', '', fqname)
        ## remove 'unmap' suffix
        fqname = re.sub('.unmap$', '', fqname) # deprecated

        if args['simplify_name']:
            fqname = re.sub('.not_\w+|.map_\w+', '', fqname)

        if args['smp_name']:
            fqname = args['smp_name']

        ## index name
        index_name = args.get('index_name', None)
        if index_name is None:
            index_name = AlignIndex(self.aligner, self.index).name

        ## output (fixed?!)
        if self.outdir_fixed:
            sub_outdir = self.outdir
        else:
            sub_outdir = os.path.join(self.outdir, fqname, index_name)

        ## global variables:
        out_prefix = os.path.join(sub_outdir, fqname)
        self.out_prefix = out_prefix
        self.sam = out_prefix + '.sam'
        self.bam = out_prefix + '.bam'
        self.log = out_prefix + '.log'
        self.unmap = out_prefix + '.unmap.' + self.fx_type # default

        if args['fq2'] is None:
            self.unmap1 = self.unmap2 = None
        else:
            if self.aligner == 'bowtie':
                self.unmap1 = out_prefix + '.unmap_1.' + self.fx_type
                self.unmap2 = out_prefix + '.unmap_2.' + self.fx_type
            else:
                self.unmap1 = out_prefix + '.unmap.1.' + self.fx_type
                self.unmap2 = out_prefix + '.unmap.2.' + self.fx_type


        ## update args
        self.fqname = fqname
        self.index_name = index_name
        self.sub_outdir = sub_outdir
        # self.outdir = sub_outdir

        ## create directory
        is_path(sub_outdir) # !!! create directories

        return fqname, out_prefix, self.sam, self.bam, self.log, self.unmap, self.unmap1, self.unmap2


    def check(self):
        """
        Check arguments, target file exists
        """
        args = self.args.copy()

        ## check args:
        args_pickle = self.out_prefix + '.arguments.pickle' # out_prefix
        args_file = self.out_prefix + '.arguments.txt'

        ## check files:
        chk1 = args_checker(args, args_pickle)
        chk2 = args['overwrite'] is False
        chk3 = os.path.exists(self.bam)

        args_logger(args, args_file, True) # update arguments.txt

        return all([chk1, chk2, chk3]) # arguments.pickle, changed 


    def get_index_list(self):
        """
        Search index for alignment:
        search index, parameters from command/arguments/parameters
            extra_index (priority=1)
            genome
            spikein
            to_rRNA
            to_te
            to_piRNA_cluster
            ...
        : arguments :
        1. genomme (str), spikein (str|None), align_to_rRNA (True|False)
        2. align_to_te (True|False), te_index (str|None) (mapping)
        3. extra_index (str|None)
        Return the list of index
        """
        args = self.args.copy()

        ## extra
        ##   : list of index, str
        if not args['extra_index'] is None:
            if isinstance(args['extra_index'], str):
                args['extra_index'] = [args['extra_index']]
            elif isinstance(args['extra_index'], list):
                pass
            else:
                raise Exception('{:>10} : --extra-index, expect str or list, \
                    get {}'.format('error', type(args['extra_index'])))
            index_list = args['extra_index']

        elif isinstance(args['genome'], str): # Genome(self.genome).check
            index_list = []
            ## extra index (rRNA, chrM, ...)
            if args['align_to_MT_trRNA']:
                index_list.append(
                    AlignIndex(self.aligner).search(args['genome'], 'MT_trRNA'))
            else:
                if args['align_to_rRNA']:
                    index_list.append(
                        AlignIndex(self.aligner).search(args['genome'], 'rRNA'))
                elif args['align_to_chrM']:
                    index_list.append(
                        AlignIndex(self.aligner).search(args['genome'], 'chrM'))
            ## genome
            index_list.append(AlignIndex(self.aligner).search(args['genome'], 'genome'))
            ## spikein
            if isinstance(args['spikein'], str):
                if args['align_to_MT_trRNA']:
                    index_list.append(
                        AlignIndex(self.aligner).search(args['spikein'], 'MT_trRNA'))
                else:
                    if args['align_to_rRNA']:
                        index_list.append(
                            AlignIndex(args['aligner']).search(args['spikein'], 'rRNA'))
                    elif args['align_to_chrM']:
                        index_list.append(
                            AlignIndex(args['aligner']).search(args['spikein'], 'chrM'))
                ## genome/spikein
                index_list.append(
                    AlignIndex(args['aligner']).search(args['spikein'], 'genome'))

        else:
            raise Exception('{:>10} : --genome and --extra-index not valid'.format('error'))

        ## check
        return [i for i in index_list if AlignIndex(args['aligner'], i).check]


## level-2: N fastx -> N index
class AlignNFastx(object):

    def __init__(self, **kwargs):
        """
        Top level for Alignment
        for N fastq and N index
        """
        args = args_init(kwargs, align=True)
        args['fq1'] = pickFq(args['fq1'], is_str=False)
        args['fq2'] = pickFq(args['fq2'], is_str=False)
        assert isinstance(args['fq1'], list)
        if args['fq2'] is None:
            args['fq2'] = [None] * len(args['fq1'])

        ## pair: fq1 == fq2
        if not len(args['fq1']) == len(args['fq2']):
            raise Exception('{:>10} : fq1 and fq2 not paired for PE reads'.format('error'))

        ## index
        args['index_list'] = args.get('index_list', None)
        self.args = args


    def run(self):
        args = self.args.copy()

        bam_list = []
        for fq1, fq2 in zip(args['fq1'], args['fq2']):
            args['fq1'] = fq1 # update
            args['fq2'] = fq2 # update
            fq_bams = AlignNIndex(**args).run()
            bam_list.append(fq_bams)

        return bam_list


## level-3: 1 fastx -> N index
class AlignNIndex(object):
    """
    Align 1 fq to N index

    required:
    aligner=
    fq1=
    outdir=
    index_list=

    optional:
    fq2=
    """
    def __init__(self, **kwargs):
        """
        required arguments:
        fq1, outdir, fq2, index_list
        """
        args = args_init(kwargs, align=True)
        args['fq1'] = pickFq(args['fq1'], is_str=True)
        args['fq2'] = pickFq(args['fq2'], is_str=True)
        assert isinstance(args['fq1'], str)
        index_list = args.get('index_list', None) # input

        ## get index list
        if index_list is None or len(index_list) == 0:
            args['index'] = '' # require searching for index
            args['search_index'] = True # switch on, for index searching
            index_list = AlignConfig(**args).index_list # exception 
            args['search_index'] = False # switch off !!! important

        if len(index_list) == 0:
            raise Exception('{:>10} : index not detected'.format('error'))

        args['index_list'] = index_list
        self.index_parallel = args.get('index_parallel', False) # for multiple index
        self.args = args


    def wrap_log(self, path=None):
        """
        wrap all alignment stat to one file, for each fastx
        outdir: outdir/sample/index/name.
        """
        if path is None:
            path = self.outdir

        ## staf file
        stat_log = path.rstrip('/') + '.align.txt'

        ## version-1: pandas
        # stat_files = listfiles2('*.align.json', path, True, True)
        # frames = [pd.read_json(i, orient='index') for i in stat_files]
        # df = pd.concat(frames, axis=1)
        # df.columns = df.loc['index_name', ].tolist()
        # df2 = df.T.sort_index() # or df1.transpose() # switch column and index
        # df2.to_csv(stat_log, sep='\t', index=False)

        # ## version-2: cat
        # stat_files = listfiles2('*.align.stat', path, True, True)
        # with open(stat_log, 'wt') as w:
        #     for f in stat_files:
        #         iname=$(basename $(dirname ${f})) #
        #         sname=$(basename ${f})
        #         sname=${sname/.align.stat}
        #         with open(f) as r:
        #             for line in r:
        #                 if line.startswith('#'):
        #                     line = '\t'.join([line.strip(), "index", "sample"])
        #                 else:
        #                     line = '\t'.join([line.strip(), iname, sname])
        #                 # w.write(''.join(r.readlines()))
        #                 w.write(line + '\n')

        stat_files = listfiles2('*.align.stat', path, True, True)
        with open(stat_log, 'wt') as w:
            for f in stat_files:
                with open(f) as r:
                    w.write(''.join(r.readlines()))

        return stat_log


    def run(self):
        """
        Run Alignment() one by one
        """
        args = self.args.copy()

        bam_list = []
        for i, index in enumerate(args['index_list']):
            args['index'] = index # update args: index
            args['index_name'] = str(i + 1) + '.' + AlignIndex(args['aligner'], index).name # update name

            ## single align port:
            align = AlignOneIndex(**args)
            bam_list.append(align.run()) ## run alignment

            ## update args: fq, fq2 # pass to next round
            if not self.index_parallel:
                if args['fq2'] is None:
                    args['fq'] = align.config.unmap
                else:
                    args['fq'] = args['fq1'] = align.config.unmap1
                    args['fq2'] = align.config.unmap2
            
            ## record the outdir/fqname
            smp_dir = os.path.join(align.config.outdir, align.config.fqname)

        ## sum all log files
        stat_log = self.wrap_log(smp_dir)

        return bam_list


## port for any aligner
## level-4: 1 fastx -> 1 index
class AlignOneIndex(object):
    """
    Align 1 fq to 1 index
    """
    def __init__(self, **kwargs):
        args = args_init(kwargs, align=True) # init
        args['fq1'] = pickFq(args['fq1'], is_str=True)
        args['fq2'] = pickFq(args['fq2'], is_str=True)
        assert isinstance(args['fq1'], str)

        index = args.get('index', None)
        if isinstance(index, list):
            log.warning('Choose the first index for alignment, {}'.format(index))
            index = index[0]
        assert isinstance(index, str) # required

        ## pick aligner
        self.args = args
        self.aligner_pick = self.get_aligner() # which aligner class()
        self.config = self.aligner_pick(**args).config # config


    def get_aligner(self):
        """
        Choose Aligner() class
        """
        args = self.args.copy()
        aligner = args['aligner'].lower()

        if aligner == 'bowtie':
            port = Bowtie # class
        elif aligner == 'bowtie2':
            port = Bowtie2 # class
        elif aligner == 'star':
            port = Star # class
        elif aligner == 'bwa':
            port = Bwa # class
        elif aligner == 'hisat2':
            port = Hisat2 # class
        else:
            raise Exception('{:>10} : aligner not supported {}'.format(args['aligner']))

        return port


    def run(self):
        return self.aligner_pick(**self.args).run()

##------ Align core: END ------##


## base-level: aligner
## align, report in json
## return bam
## object().config
## add aligners, modify AlignOneIndex().get_aligner()
##
## save all configs to self.object (dict)

class Bowtie(object):
    def __init__(self, **kwargs):
        """
        required arguments:
        fq1, single fastq file, or read1 of PE
        index, path to the index file
        outdir, path to the output directory
        fq2, optional, read2 of PE
        ...        
        """
        args = args_init(kwargs, align=True) # init
        args['fq1'] = pickFq(args['fq1'], is_str=True)
        args['fq2'] = pickFq(args['fq2'], is_str=True)
        assert isinstance(args['index'], str)
        assert isinstance(args['fq1'], str)

        ## update args: aligner, fq, index, outdir, fq2, outdir_fixed
        args['aligner'] = 'bowtie' # fix
        outdir_fixed = False # whether save to: outdir/fqname 

        ## saving common configs for alignment
        self.args = args
        self.config = AlignConfig(**args)
        self.aligner_exe = shutil.which('bowtie')


    def get_cmd(self):
        """
        Generate the command line for bowtie:
        1. bowtie [options] -U <fq> 
        2. bowtie [options] -1 fq1 -2 fq2
        """
        args = self.args.copy()

        ## multi map
        if args['n_map'] < 1:
            args['n_map'] = 1 # bowtie default: -k 1

        ## unique
        if args['unique_only']:
            cmd_unique = '-m 1' # unique
        else:
            cmd_unique = '-v 2 -k {}'.format(args['n_map']) # multi

        ## fastq type
        cmd_fx = '-f' if self.config.format == 'fasta' else '-q'

        ## cmd
        cmd = '{} {} {} {} -p {} --mm --best --sam --no-unal --un {}'.format(
            self.aligner_exe,
            args['index'],
            cmd_fx,
            cmd_unique,
            args['threads'],
            self.config.unmap)

        ## extra align para
        if not args['extra_para'] is None:
            cmd += ' ' + args['extra_para']

        # print('!aaaa', args)

        ## se or pe
        if args['fq2'] is None: 
            cmd += ' {} 1>{} 2>{}'.format(
                args['fq1'],
                self.config.sam,
                self.config.log)
        else:
            cmd += ' -1 {} -2 {} 1> {} 2>{}'.format(
                args['fq1'],
                args['fq2'],
                self.config.sam,
                self.config.log)

        return cmd


    def wrap_log(self):
        """
        Wrapper bowtie log

        Bowtie:
        # reads processed: 10000
        # reads with at least one reported alignment: 3332 (33.32%)
        # reads that failed to align: 457 (4.57%)
        # reads with alignments suppressed due to -m: 6211 (62.11%)

        or:

        # reads processed: 10000
        # reads with at least one reported alignment: 9543 (95.43%)
        # reads that failed to align: 457 (4.57%)

        unique, multiple, unmap, map, total

        skip: Warning, ...
        """
        args = self.args.copy()

        dd = {}
        with open(self.config.log) as fh:
            for line in fh:
                # if not ':' in line or line.startswith('Warning'):
                #     continue
                if not line.startswith('#'):
                    continue
                num = line.strip().split(':')[1]
                value = num.strip().split(' ')[0]
                value = eval(value)
                if 'reads processed' in line:
                    dd['total'] = value
                elif 'at least one reported alignment' in line:
                    dd['map'] = value
                elif 'failed to align' in line:
                    dd['unmap'] = value
                elif 'alignments suppressed due to -m' in line:
                    dd['multiple'] = value
                else:
                    pass
        
        # unique_only
        dd['unique'] = dd['map']
        dd['multiple'] = dd.get('multiple', 0) # default 0
        if args['unique_only']:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']

        # save fqname, indexname, 
        dd['fqname'] = self.config.fqname
        dd['index_name'] = self.config.index_name

        # # sort by keys
        self.log_dict = dd

        # save dict to plaintext file
        self.align_stat = self.config.out_prefix + '.align.stat'
        with open(self.align_stat, 'wt') as w:
            w.write('#') # header line
            w.write('\t'.join(list(map(str, dd.keys()))) + '\n')
            w.write('\t'.join(list(map(str, dd.values()))) + '\n')

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def get_json(self):
        """
        save alignment stat to Json file
        """
        log_json = self.config.out_prefix + '.align.json'
        self.align_stat = self.wrap_log() # update: self.log_dict
        ## save to json file
        Json(self.log_dict).writer(log_json)

        return log_json


    def run(self):
        """
        Run bowtie command line
        """
        args = self.args.copy()
        cmd = self.get_cmd() # output

        ## para, bam, overwrite
        if self.config.check_status:
            log.info('{:>20} : file exists, alignment skipped'.format(self.config.fqname))
        else:
            # try:
            run_shell_cmd(cmd)
            sam2bam(self.config.sam, self.config.bam, sort=True) # convert to bam
            self.get_json() # log_json
            # except:
            #     log.error('Bowtie().run() failed, outdir: {}'.format(args['outdir']))

        return self.config.bam
        

class Bowtie2(object):

    def __init__(self, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...
        """
        args = args_init(kwargs, align=True) # init
        args['fq1'] = pickFq(args['fq1'], is_str=True) # only 1 fq
        args['fq2'] = pickFq(args['fq2'], is_str=True) # only 1 fq
        assert isinstance(args['index'], str)
        assert isinstance(args['fq1'], str)

        ## update args: aligner, fq, index, outdir, fq2, outdir_fixed
        args['aligner'] = 'bowtie2' # fix
        outdir_fixed = False # whether save to: outdir/fqname 

        ## saving common configs for alignment
        self.args = args
        self.config = AlignConfig(**args)
        self.aligner_exe = shutil.which('bowtie2')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. bowtie [options] -U <fq> 
        2. bowtie [options] -1 fq1 -2 fq2
        """
        args = self.args.copy()

        ## multi map
        n_map = args.get('n_map', 1)
        if n_map < 1:
            cmd_multi = '' # default: best hit
        else:
            cmd_multi = '-k {}'.format(n_map)

        ## fx type
        cmd_fx = '-f' if self.config.format == 'fasta' else '-q'

        ## cmd
        cmd = '{} -x {} {} -p {} {} --sensitive-local --mm --no-unal'.format(
            self.aligner_exe,
            args['index'],
            cmd_fx,
            args['threads'],
            cmd_multi)

        ## extra align para
        ## eg: -X 2000,
        if not args['extra_para'] is None:
            cmd += ' ' + args['extra_para']

        ## se or pe
        if args['fq2'] is None: 
            cmd += ' --un {} {} 1>{} 2>{}'.format(
                self.config.unmap,
                args['fq1'],
                self.config.sam,
                self.config.log)
        else:
            cmd += ' --un-conc {} -1 {} -2 {} 1>{} 2>{}'.format(
                self.config.unmap,
                args['fq1'],
                args['fq2'],
                self.config.sam,
                self.config.log)

        return cmd


    def wrap_log(self):
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
        args = self.args.copy()

        dd = {}
        se_tag = 1 #
        with open(self.config.log, 'rt') as ff:
            for line in ff:
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

        if args['unique_only']:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']

        # save fqname, indexname, 
        dd['fqname'] = self.config.fqname
        dd['index_name'] = self.config.index_name

        # save dict
        # sort by keys
        # dd = dict(sorted(dd.items(), key=lambda kv: kv[1], reverse=True))
        self.log_dict = dd

        # save dict to plaintext file
        self.align_stat = self.config.out_prefix + '.align.stat'
        with open(self.align_stat, 'wt') as w:
            # for k, v in sorted(dd.items()):
            #     w.write('\t'.join([self.config.fqname, self.config.index_name, k, str(v)]) + '\n')
            w.write('#') # header line
            w.write('\t'.join(list(map(str, dd.keys()))) + '\n')
            w.write('\t'.join(list(map(str, dd.values()))) + '\n')

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def get_json(self):
        """
        get alignment records as Json
        """
        log_json = self.config.out_prefix + '.json'

        self.wrap_log() # self.log_dict

        ## save to json file
        Json(self.log_dict).writer(log_json)

        return log_json


    def get_unique(self, bam):
        """
        Get the unique mapped reads, using samtools -q 30
        Move bam to bam.tmp
        """
        bam_old = os.path.splitext(bam)[0] + '.raw.bam'
        if os.path.exists(bam):
            shutil.move(bam, bam_old)
        # unique: -q 30
        sam2bam(bam_old, bam, sort=True, extra_para='-q 30')
        return(bam)


    def run(self):
        """
        Run bowtie command line
        """
        args = self.args.copy()
        cmd = self.get_cmd() # output

        ## para, bam, overwrite
        if self.config.check_status:
            log.info('{:>20} : file exists, alignment skipped'.format(
                self.config.fqname))
        else:
            # try:
            run_shell_cmd(cmd)
            if args['unique_only']:
                self.get_unique(self.config.bam)
            else:
                sam2bam(self.config.sam, self.config.bam, sort=True, 
                    extra_para='-F 4')
            self.get_json() # save to json
            # except:
            #     log.error('Bowtie2().run() failed, outdir: {}'.format(
            #         args['outdir']))

        return self.config.bam


class Star(object):
    def __init__(self, **kwargs):
        """
        required arguments:
        fq1, single fastq file, or read1 of PE
        index, path to the index file
        outdir, path to the output directory
        fq2, optional, read2 of PE
        ...        
        """
        args = args_init(kwargs, align=True) # init
        args['fq1'] = pickFq(args['fq1'], is_str=True)
        args['fq2'] = pickFq(args['fq2'], is_str=True)
        assert isinstance(args['index'], str)
        assert isinstance(args['fq1'], str)

        ## update args: aligner, fq, index, outdir, fq2, outdir_fixed
        args['aligner'] = 'STAR' # fix
        outdir_fixed = False # whether save to: outdir/fqname 

        ## saving common configs for alignment
        self.args = args
        self.config = AlignConfig(**args)
        self.aligner_exe = shutil.which('STAR')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. STAR
        """
        args = self.args.copy()

        ## for small genome
        ## mapping takes too long,
        ## 99% of the reads not mapped
        ## change --seedPerWindowNmax 
        ## max number of seeds per window
        ## https://github.com/alexdobin/STAR/issues/329
        ## https://groups.google.com/d/msg/rna-star/hJL_DUtliCY/HtpiePlMBtYJ 
        if args['small_genome']:
            seed_max = 5 # even smaller
        else:
            seed_max = 50 # default

        ## for unique map
        ## --outFilterMultimapNmax
        ## maximum number of loci the read is allowed to map to. [default: 10]
        ##
        ## filt unique reads by: samtools view -q 30
        n_map = args.get('n_map', 10)
        if n_map == 0:
            n_map = 10 # default
        if args['unique_only']:
            n_map = 1 # unique # not working
        cmd_unique = '--outFilterMultimapNmax {} --seedPerWindowNmax {}'.format(
            n_map, seed_max)
        
        ## file type
        cmd_reader = 'zcat' if args['fq1'].endswith('.gz') else '-'

        ## convert None to empty string ''
        if args['fq2'] is None:
            args['fq2'] = ''

        ## output prefix
        out_prefix = self.config.out_prefix
        cmd = '{} --runMode alignReads \
            --genomeDir {} \
            --readFilesIn {} {} \
            --readFilesCommand {} \
            --outFileNamePrefix {} \
            --runThreadN {} \
            --limitOutSAMoneReadBytes 1000000 \
            --genomeLoad NoSharedMemory  \
            --limitBAMsortRAM 10000000000 \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMismatchNoverLmax 0.07 \
            --seedSearchStartLmax 20 \
            --outReadsUnmapped Fastx {} {}'.format(
                self.aligner_exe,
                args['index'],
                args['fq1'],
                args['fq2'],
                cmd_reader,
                self.config.out_prefix,
                args['threads'],
                self.config.unmap,
                cmd_unique)

        ## extra align para
        if not args['extra_para'] is None:
            cmd += ' ' + args['extra_para']

        return cmd


    def update_names(self, keep_old=False):
        """
        Update the filenames of STAR output
        bam: *Aligned.sortedByCoord.out.bam -> *.bam # mv
        log: *Log.final.out -> *.log # copy
        log: *Log.out -> *.out
        unmap: *Unmapped.out.mate1 -> *.unmap.1.fastq
               *Unmapped.out.mate1 -> *.unmap.1.fastq
        """
        args = self.args.copy()

        ## default output of STAR
        bam_from = self.config.out_prefix + 'Aligned.sortedByCoord.out.bam'
        log_from = self.config.out_prefix + 'Log.final.out'
        unmap = unmap1 = self.config.out_prefix + 'Unmapped.out.mate1'
        unmap2 = self.config.out_prefix + 'Unmapped.out.mate2'

        if keep_old:
            return bam_from
        else:
            ## rename files
            ## move BAM, unmap
            ## copy log
            flag = 0
            if os.path.exists(bam_from) and not os.path.exists(self.config.bam):
                # os.symlink(os.path.basename(bam_from), self.config.bam) # symlink
                shutil.move(bam_from, self.config.bam)
                flag += 1
            if os.path.exists(log_from):
                shutil.copy(log_from, self.config.log)
                flag += 1

            ## why self.fq2 changed, from None to '' !!!
            if not args['fq2']:
                if os.path.exists(unmap):
                    shutil.move(unmap, self.config.unmap)
                    flag += 1
            else:
                if os.path.exists(unmap1):
                    shutil.move(unmap1, self.config.unmap1)
                    flag += 1
                if os.path.exists(unmap2):
                    shutil.move(unmap2, self.config.unmap2)

            return self.config.bam


    def wrap_log(self):
        """
        Wrapper

        STAR:
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

        unique, multiple, unmap, map, total
        """
        args = self.args.copy()

        dd = {}
        with open(self.config.log, 'rt') as ff:
            for line in ff:
                value = line.strip().split('|')
                if not len(value) == 2:
                    continue
                value = value[1].strip()
                if 'Number of input reads' in line:
                    dd['total'] = int(value)
                elif 'Uniquely mapped reads number' in line:
                    dd['unique'] = int(value)
                elif 'Number of reads mapped to multiple loci' in line:
                    dd['multiple'] = int(value)
                else:
                    pass

        if args['unique_only']:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']        

        # save fqname, indexname, 
        dd['fqname'] = self.config.fqname
        dd['index_name'] = self.config.index_name

        # sort by keys
        # dd = dict(sorted(d.items(), key=lambda kv: kv[1], reverse=True))
        self.log_dict = dd

        # save dict to plaintext file
        self.align_stat = self.config.out_prefix + '.align.stat'
        with open(self.align_stat, 'wt') as w:
            # for k, v in sorted(dd.items()):
            #     w.write('\t'.join([self.config.fqname, self.config.index_name, k, str(v)]) + '\n')
            w.write('#') # header line
            w.write('\t'.join(list(map(str, dd.keys()))) + '\n')
            w.write('\t'.join(list(map(str, dd.values()))) + '\n')

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def get_json(self):
        """
        get alignment records as Json
        """
        log_json = self.config.out_prefix + '.json'
        self.wrap_log() # self.log_dict
        ## save to json file
        Json(self.log_dict).writer(log_json)
        return log_json


    def get_unique(self, bam):
        """
        Get the unique mapped reads, using samtools -q 30
        Move bam to bam.tmp
        """
        bam_old = os.path.splitext(bam)[0] + '.raw.bam'
        if os.path.exists(bam):
            shutil.move(bam, bam_old)
        # unique: -q 30
        sam2bam(bam_old, bam, sort=True, extra_para='-q 30')
        return(bam)


    def run(self):
        """
        Run STAR command line
        """
        args = self.args.copy()
        cmd = self.get_cmd() # output

        ## para, bam, overwrite
        if self.config.check_status:
            log.info('{:>20} : file exists, alignment skipped'.format(
                self.config.fqname))
        else:
            # try:
            run_shell_cmd(cmd)
            self.update_names(keep_old=False)
            if args['unique_only']:
                self.get_unique(self.config.bam)
            self.get_json() # save to json
            # except:
            #     log.error('Star().run() failed, outdir: {}'.format(
            #         args['outdir']))

        return self.config.bam      


## !!! to-do: wrap_log()
class Hisat2(object):

    def __init__(self, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...

        """
        args = args_init(kwargs, align=True) # init
        args['fq1'] = pickFq(args['fq1'], is_str=True) # only 1 fq
        args['fq2'] = pickFq(args['fq2'], is_str=True) # only 1 fq
        assert isinstance(args['index'], str)
        assert isinstance(args['fq1'], str)

        ## update args: aligner, fq, index, outdir, fq2, outdir_fixed
        args['aligner'] = 'hisat2' # fix
        outdir_fixed = False # whether save to: outdir/fqname 

        ## saving common configs for alignment
        self.args = args
        self.config = AlignConfig(**args)
        self.aligner_exe = shutil.which('hisat2')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. hisat2 [options] -U <fq> 
        2. hisat2 [options] -1 fq1 -2 fq2
        """
        args = self.args.copy()

        # ## unique
        # if args['unique_only']:
        #     arg_unique = '-q 30' # samtools filtering
        # else:
        #     arg_unique = ''

        # self.arg_unique = arg_unique # pass to other func

        ## multi map
        n_map = args.get('n_map', 1)
        if n_map < 1:
            # n_map = 1 # default 1, report 1 hit for each read
            cmd_multi = '' # default: 1
        else:
            cmd_multi = '-k {}'.format(n_map)

        ## fx type
        cmd_fx = '-f' if self.config.format == 'fasta' else '-q'

        ## cmd
        cmd = '{} -x {} {} {} -p {} --no-discordant --mm --new-summary'.format(
            self.aligner_exe,
            args['index'],
            cmd_fx,
            cmd_multi,
            args['threads'])

        ## se or pe
        if args['fq2'] is None: 
            cmd += ' --un {} {} 1>{} 2>{}'.format(
                self.config.unmap,
                args['fq1'],
                self.config.sam,
                self.config.log)
        else:
            cmd += ' --un-conc {} -1 {} -2 {} 1>{} 2>{}'.format(
                self.config.unmap,
                args['fq1'],
                args['fq2'],
                self.config.sam,
                self.config.log)

        return cmd


    def wrap_log(self):
        """
        Wrapper Hisat2 log

        Hisat2

        unique, multiple, unmap, map, total
        """
        self.log_dict = {}
        pass # !!! to-do


    def get_json(self):
        """
        get alignment records as Json
        """
        log_json = self.config.out_prefix + '.json'
        self.wrap_log() # self.log_dict
        ## save to json file
        Json(self.log_dict).writer(log_json)

        return log_json


    def get_unique(self, bam):
        """
        Get the unique mapped reads, using samtools -q 30
        Move bam to bam.tmp
        """
        bam_old = os.path.splitext(bam)[0] + '.raw.bam'
        if os.path.exists(bam):
            shutil.move(bam, bam_old)
        # unique: -q 30
        sam2bam(bam_old, bam, sort=True, extra_para='-q 30')
        return(bam)


    def run(self):
        """
        Run hisat2 command line
        """
        args = self.args.copy()
        cmd = self.get_cmd() # output

        ## para, bam, overwrite
        if self.config.check_status:
            log.info('{:>20} : file exists, alignment skipped'.format(
                self.config.fqname))
        else:
            # try:
            run_shell_cmd(cmd)
            if args['unique_only']:
                self.get_unique(self.config.bam)
            else:
                sam2bam(self.config.sam, self.config.bam, sort=True, 
                    extra_para='-F 4')
            self.get_json() # save to json
            # except:
            #     log.error('Hisat2().run() failed, outdir: {}'.format(
            #         args['outdir']))

        return self.config.bam


## not tested
## !!! to-do
class Bwa(object):

    def __init__(self, fq, index, outdir, fq2=None, **kwargs):
        """
        Align reads to reference
        return: 
        - sam
        - align count
        - unique
        - multi
        - unmap
        - reference
        - ...

        """
        assert os.path.exists(fq)
        assert AlignIndex(index).is_index('bwa')
        is_path(outdir)

        ## init
        args = args_init(kwargs, align=True) # init

        ## output
        fqname = file_prefix(fq)[0]
        if not fq2 is None:
            fqname = re.sub('_[rR]?1$', '', fqname)
        align_prefix = os.path.join(outdir, fqname)
        
        self.fq = fq
        self.fq1 = fq
        self.fq2 = fq2
        self.index = index
        self.outdir = outdir
        self.fqname = fqname
        self.kwargs = args

        self.sam = align_prefix + '.sam'
        self.bam = align_prefix + '.bam'
        self.log = align_prefix + '.log'
        self.unmap = align_prefix + '.unmap.fq' # fa/fq
        if Fastx(fq).format() == 'fasta':
            self.unmap = align_prefix + '.unmap.fa' # fa/fq

        self.aligner_exe = shutil.which('bowtie2')


    def get_cmd(self):
        """
        Create basic alignment command line
        1. bowtie [options] -U <fq> 
        2. bowtie [options] -1 fq1 -2 fq2
        """
        args = self.kwargs.copy()

        n_map = args.get('n_map', 0)

        ## multi map
        if n_map < 1:
            n_map = 1

        ## unique
        if args['unique_only']:
            arg_unique = '-m 1'
        else:
            arg_unique = '-v 2 -k {}'.format(n_map)

        arg_fx = '-f' if self.fq_type == 'fasta' else '-q'

        cmd = '        '

        return cmd


    def align(self):
        """
        Run bowtie command line
        """
        args = self.kwargs.copy()

        cmd = self.get_cmd()

        if os.path.exists(self.bam) and args['overwrite'] is False:
            log.info('{:>20} : file exists, alignment skipped'.format(self.fqname))
        else:
            run_shell_cmd(cmd)

        ## convert sam to bam
        sam2bam(self.sam, self.bam)

        return self.bam


## for one fq, all index
## deprecated: AlignConfig(search_index=True, **args).index_list
class AlignIndexBuilder(object):
    """
    Search index :
        for alignment
        search index, parameters from command/arguments/parameters
        extra_index (priority=1)
        genome
        spikein
        to_rRNA
        to_te
        to_piRNA_cluster
        ...
    : arguments :
    1. genomme (str), spikein (str|None), align_to_rRNA (True|False)
    2. align_to_te (True|False), te_index (str|None) (mapping)
    3. extra_index (str|None)
    """
    def __init__(self, aligner='bowtie', **kwargs):

        ## args init
        args = args_init(kwargs, align=True)

        ## update arguments
        self.genome = args.get('genome', 'dm3')
        self.spikein = args.get('spikein', 'dm3')
        self.align_to_rRNA = args.get('align_to_rRNA', True)
        self.extra_index = args.get('extra_index', None)
        self.aligner = aligner
        self.kwargs = args
        self.index_list = self.builder()


    def builder(self):
        """
        Return the list of index
        extra
        genome (rRNA) 
        spikein (rRNA)
        ...
        """
        args = self.kwargs.copy()

        ## extra
        ##   : list of index, str
        if not self.extra_index is None:
            if isinstance(self.extra_index, str):
                self.extra_index = [self.extra_index]
            elif isinstance(self.extra_index, list):
                pass
            else:
                raise Exception('{:>10} : --extra-index, expect str or list, \
                    get {}'.format('error', type(self.extra_index)))
            index_list = self.extra_index
        elif isinstance(self.genome, str): # Genome(self.genome).check
            index_list = []
            ## genome
            if self.align_to_rRNA:
                index_list.append(
                    AlignIndex(self.aligner).search(self.genome, 'rRNA'))
            index_list.append(AlignIndex(self.aligner).search(self.genome, 'genome'))
            ## spikein
            if isinstance(self.spikein, str):
                if self.align_to_rRNA:
                    index_list.append(
                        AlignIndex(self.aligner).search(self.spikein, 'rRNA'))
                index_list.append(
                    AlignIndex(self.aligner).search(self.spikein, 'genome'))

        else:
            raise Exception('{:>10} : --genome and --extra-index not valid'.format('error'))

        ## check
        return [i for i in index_list if AlignIndex(self.aligner, i).check]

