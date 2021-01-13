

import os
import re
# import hiseq
# import pysam
import shutil
import tempfile
# import pandas as pd
# from multiprocessing import Pool
import copy # copy objects
from hiseq.utils.helper import *


## design, from txt file
class DesignReader(object):
    """
    Parsing tab file for arguments
    return dict
    optional:
    return: json file, dict

    rnaseq:
    names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
    rnaseq:
    names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']

    atacseq:
    names=['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']
    """
    def __init__(self, file):
        self.file = file
        self.df = self.read_design(file)
        self.hiseq_type, self.hiseq_subtype = self.design_type()
        self.args = self.to_dict()


    def read_design(self, x):
        try:
            df = pd.read_csv(x, '\t', header=None, comment='#')
            df = df.fillna('NULL')
        except:
            df = pd.DataFrame()

        return df


    def design_type(self):
        """
        check the 1st column: RNAseq | ATACseq
        check the 7th column: isdir(), True | False
        check ncolumns == 8: True | False

        rnaseq:
        names=['RNAseq', ...]
        atacseq:
        names=['ATACseq', ...]
        """
        if len(self.df) == 0:
            # empty design
            log.warning('--design, no record detected. {}'.format(self.file))
            return None, None

        name_list = set(self.df.iloc[:, 0].to_list()) # 1-st column

        if len(name_list) > 1:
            raise Exception('Multiple seq-type found in 1-st column: {}'.format(name_list))

        hiseq = list(name_list).pop().lower() # rnaseq | atacseq

        # the 7th column: file
        chk1x = [os.path.isfile(i) for i in self.df.iloc[:, 6].to_list()]
        chk1 = all(chk1x)

        # log
        chk1_log = []
        for b, f in zip(chk1x, self.df.iloc[:, 6].to_list()):
            chk1_log.append('{:<6} : {:<s}'.format(str(b), f))

        if hiseq == 'atacseq' and not chk1:
            raise Exception('Design not correct: {}\n{}'.format(self.file, '\n'.join(chk1_log)))

        # ncolumns
        chk2 = len(self.df.columns) == 8

        if chk2: # RNAseq
            hiseq_sub = 'multiple'
        else:
            hiseq_sub = 'deseq' if hiseq == 'rnaseq' and chk1 else None

        return (hiseq, hiseq_sub)


    def to_dict(self):
        if self.hiseq_type == 'rnaseq':
            if self.hiseq_subtype == 'multiple':
                self.df.columns = ['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
                args = {
                    'fq1': self.df['fq1'].to_list(),
                    'fq2': self.df['fq2'].to_list()
                }
                chk0 = all([os.path.exists(i) for i in args['fq1']]) # check
            elif self.hiseq_subtype == 'deseq':
                self.df.columns = ['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']
                args = {
                    'smp_path': self.df['smp_path'].to_list()
                }
                chk0 = len(args['smp_path']) == len(set(args['smp_path'])) # check
            else:
                args = {} # pass
                chk0 = True # check
            # feature: single
            features = self.df['feature'].to_list()
            chk1 = len(set(features)) == 1 # check
        elif self.hiseq_type == 'atacseq':
            self.df.columns = ['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']
            args = {
                    'fq1': self.df['fq1'].to_list(),
                    'fq2': self.df['fq2'].to_list()
                }
            chk0 = all([os.path.exists(i) for i in args['fq1']]) # check
            chk1 = True # check
        else:
            args = {} # pass
            chk0 = chk1 = True # check
            log.warning("""
                Format failed:
                --design: {},
                Expect format:
                rnaseq1:
                ['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
                rnaseq2:
                ['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']
                atacseq:
                ['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']""".format(self.file))
            return {} # empty dict


        ## check the df
        # unique names
        names = self.df['name'].to_list()
        # update name
        names_new = []
        for i, n in enumerate(names):
            if n.lower() in ['null', 'na', 'nan', '']:
                if 'smp_path' in self.df.columns:
                    n_new = os.path.basename(self.df['smp_path'].to_list()[i])
                else:
                    n_new = fq_name(self.df['fq1'].to_list()[i], pe_fix=True)
            else:
                n_new = n
            names_new.append(n_new)

        chk2 = len(names_new) == len(set(names_new))

        # unique genome: single, str
        genomes = self.df['genome'].to_list()
        chk3 = len(set(genomes)) == 1

        # outdir: single, str
        outdirs = self.df['outdir'].to_list()
        chk4 = len(set(outdirs)) == 1

        # update args
        args['smp_name'] = names_new
        args['genome'] = genomes.pop()
        args['outdir'] = outdirs.pop()

        # for check
        if not all([chk0, chk1, chk2, chk3, chk4]):
            raise Exception(""" Design file format not correct:
                fq1 exists, smp_path exists: {}
                             feature unique: {}
                               names unique: {}
                              genome unique: {}
                              outdir unique: {}""".format(chk0, chk1, chk2, chk3, chk4))
        return args


    # deprecated #
    def design_rnaseq(self):
        """
        Add headers to RNAseq design
        names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
        names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']
        """
        ncols = len(self.df.columns)
        if ncols == 8:
            self.df.columns = ['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
            # fq1: str
            fq1_list = self.df['fq1'].to_list()
            fq2_list = self.df['fq2'].to_list()
            chk0 = [os.path.exists(i) for i in fq1_list]

            # construct args
            args = {
                'fq1': fq1_list,
                'fq2': fq2_list
            }

        elif ncols == 7:
            self.df.columns = ['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'count_txt']
            # count_txt path
            smp_path = self.df['smp_path'].to_list()
            chk0 = len(smp_path) == len(set(smp_path))

            # construct args
            args = {
                'smp_path': smp_path
            }

        else:
            raise Exception('unknown design for RNAseq: {}'.format(self.file))

        ## check the df
        # unique names
        names = self.df['name'].to_list()
        chk1 = len(names) == len(set(names))
        chk1 = True # !!!! skip

        # unique genome: single, str
        genomes = self.df['genome'].to_list()
        chk2 = len(set(genomes)) == 1

        # outdir: single, str
        outdirs = self.df['outdir'].to_list()
        chk3 = len(set(outdirs)) == 1

        # feature: single
        features = self.df['feature'].to_list()
        chk4 = len(set(features)) == 1

        # update args
        args['smp_name'] = names
        args['genome'] = genomes.pop()
        args['outdir'] = outdirs.pop()
        args['feature'] = features.pop()

        # for check
        if not all([chk0, chk1, chk2, chk3, chk4]):
            raise Exception('file format not correct: {}'.format(self.file))
            raise Exception(""" Design file format not correct:
                fq1 or smp_list: {}
                          names: {}
                         genome: {}
                         outdir: {}
                        feature: {}""".format(chk0, chk1, chk2, chk3, chk4))



        return args


    # deprecated #
    def design_atacseq(self):
        """
        atacseq:
        names=['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']
        """
        ncols = len(self.df.columns)
        if ncols == 7:
            self.df.columns = ['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']
            # fq1: str
            fq1_list = self.df['fq1'].to_list()
            fq2_list = self.df['fq2'].to_list()
            chk0 = [os.path.exists(i) for i in fq1_list]

            # construct args
            args = {
                'fq1': fq1_list,
                'fq2': fq2_list,
            }
        else:
            raise Exception('unknown design for ATACseq: {}'.format(self.file))

        ## check the df
        # unique names
        names = self.df['name'].to_list()
        chk1 = len(names) == len(set(names))

        # unique genome: single, str
        genomes = self.df['genome'].to_list()
        chk2 = len(set(genomes)) == 1

        # outdir: single, str
        outdirs = self.df['outdir'].to_list()
        chk3 = len(set(outdirs)) == 1

        # update args
        args['smp_name'] = names
        args['genome'] = genomes.pop()
        args['outdir'] = outdirs.pop()

        # for check
        if not all([chk0, chk1, chk2, chk3]):
            raise Exception('file format not correct: {}'.format(self.file))

        return args


    def to_json(self, json_file=None):
        if json_file is None:
            tmp_prefix = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
            json_file = '{}.design.json'.format(tmp_prefix)

        with open(json_file, 'wt') as fo:
            json.dump(self.d, fo, indent=4, sort_keys=True)

        return json_file


## design, to txt file
class DesignBuilder(object):
    """
    Create design.txt for pipeline
    RNAseq:
    fq1, genome, outdir, smp_name, feature, group, count_txt, fq2

    ATACseq:
    fq1, genome, outdir, smp_name, fq2

    rnaseq:
    names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
    rnaseq:
    names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']

    atacseq:
    names=['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.init_args()        

        if len(kwargs) > 0:
            self.hiseq_type, self.hiseq_subtype = self.mission()
            self.status = self.check() # updated


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
        local_config = RNAseqConfig(**self.__dict__) # pre-defined
        self.update(local_config.__dict__)

        # for fq2
        if self.fq2 is None:
            if isinstance(self.fq1, list):
                self.fq2 = [None] * len(self.fq1)
            elif isinstance(self.fq1, str):
                self.fq2 = None
            else:
                self.fq2 = None


    def demo(self):
        """
        Create demo file for example
        """
        lines = ['\t'.join(['#RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2'])]
        lines.append('\t'.join(['#RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']))
        lines.append('\t'.join(['#ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']))
        print('\n'.join(lines))


    def mission(self):
        """
        Check the type of analysis
        RNAseq, or ATACseq, ...
        """
        if self.feature is None:
            tag = 'ATACseq'
            tag_sub = None
        else:
            tag = 'RNAseq'
            if isinstance(self.smp_path, list):
                tag_sub = 'deseq_multiple'
            elif isinstance(self.dirs_ctl, list) and isinstance(self.dirs_exp, list):
                tag_sub = 'deseq_single'
            elif isinstance(self.fq1, list):
                tag_sub = 'rnaseq_multiple'
            elif isinstance(self.fq1, str):
                tag_sub = 'rnaseq_single'
            else:
                tag_sub = None

        return (tag, tag_sub)


    def check(self):
        """
        Check the arguments
        RNAseq, count_txt/fq1,
        ATACseq, fq1
        smp_name ('NULL'|[list])
        group ('NULL' | [list - by fq1/count])
        """
        ## RNAseq
        if self.hiseq_type == 'RNAseq':
            # feature
            chk0 = isinstance(self.feature, str)
            # rnaseq type
            if self.hiseq_subtype == 'deseq_multiple':
                n_samples = len(self.smp_path)
                chk1 = len(self.smp_path) == len(set(self.smp_path))
            elif self.hiseq_subtype == 'deseq_single':
                n_samples = len(self.smp_path)
                self.smp_path = self.dirs_ctl + self.dirs_exp
                chk1 = len(self.smp_path) == len(set(self.smp_path))
            elif self.hiseq_subtype == 'rnaseq_multiple':
                n_samples = len(self.fq1)
                self.fq1 = self.fq1 if isinstance(self.fq1, list) else [self.fq1]
                self.fq2 = self.fq2 if isinstance(self.fq2, list) else [self.fq2]
                chk1 = True
            elif self.hiseq_subtype == 'rnaseq_single':
                n_samples = len(self.fq1)
                if self.fq2 is None:
                    self.fq2 = 'NULL'
                chk1 = True
            else:
                n_samples = 0
                chk1 = True
                pass

        ## ATACseq
        elif self.hiseq_type == 'ATACseq':
            chk0 = chk1 = True  # tmp
            # # check fq
            # self.fq1 = self.fq1 if isinstance(self.fq1, list) else [self.fq1]
            # n_samples = len(self.fq1)
            # if self.fq2 is None:
            #     self.fq2 = ['NULL'] * len(self.fq1)
            # else:
            #     self.fq2 = self.fq2 if isinstance(self.fq2, list) else [self.fq2]

        ## others
        else:
            chk0 = chk1 = True  # tmp
            raise Exception('unkonwn input arguments')

        # genomes
        chk1 = isinstance(self.genome, str)

        # outdir
        chk2 = isinstance(self.outdir, str)

        # groups
        chk3 = isinstance(self.group, list) and len(self.group) == n_samples

        if not all([chk0, chk1, chk2, chk3]):
            raise Exception(""" Design file format not correct:
                                    feature: {}
                                     genome: {}
                                     outdir: {}
                                      group: {}, {}""".format(
                                        chk0, 
                                        self.genome, 
                                        self.outdir, 
                                        n_samples, self.group))


    def to_txt(self):
        """
        Convert to text file
        atacseq:
        names=['ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']

        Add headers to RNAseq design
        names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']
        names=['RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']
        """
        design_lines = []
        if self.hiseq_type == 'RNAseq':
            if self.hiseq_subtype in ['rnaseq_single', 'rnaseq_multiple']:
                design_lines.append('\t'.join([
                    '#RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'fq1', 'fq2']))
                for i, q in enumerate(self.fq1):
                    design_lines.append('\t'.join([
                        'RNAseq',
                        self.group[i],
                        self.smp_name[i],
                        self.feature,
                        self.genome,
                        self.outdir,
                        q,
                        str(self.fq2[i])]))
            elif self.hiseq_subtype in ['deseq_single', 'deseq_multiple']:
                design_lines.append('\t'.join([
                    '#RNAseq', 'group', 'name', 'feature', 'genome', 'outdir', 'smp_path']))
                for i, q in enumerate(self.smp_path):
                    design_lines.append('\t'.join([
                        'RNAseq',
                        self.group[i],
                        self.smp_name[i],
                        self.feature,
                        self.genome,
                        self.outdir,
                        q]))
            else:
                log.warning('unknown RNAseq type: {}'.format(self.hiseq_subtype))

        elif self.hiseq_type == 'ATACseq':
            design_lines.append('\t'.join([
                '#ATACseq', 'group', 'name', 'genome', 'outdir', 'fq1', 'fq2']))
            for i, q in enumerate(self.fq1):
                design_lines.append('\t'.join([
                    'ATACseq',
                    self.group[i],
                    self.smp_name[i],
                    self.genome,
                    self.outdir,
                    q,
                    str(self.fq2[i])]))
        else:
            pass

        return design_lines


    def to_file(self, file=None):
        """
        Organize the arguments to text file
        """
        if file is None:
            tmp_prefix = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
            file = '{}.design.txt'.format(tmp_prefix)

        design_lines = self.to_txt()

        if os.path.exists(file):
            log.info('file exists, overwrited: {}'.format(file))
        # else:
            with open(file, 'wt') as w:
                w.write('\n'.join(design_lines) + '\n')

        return file


    def to_json(self, json_file=None):
        if json_file is None:
            tmp_prefix = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
            json_file = '{}.design.json'.format(tmp_prefix)

        with open(json_file, 'wt') as fo:
            json.dump(self.__dict__, fo, indent=4, sort_keys=True)

        return json_file



class RNAseqBuildDesign(object):
    """
    Create design.txt for experiment
    """
    def __init__(self, **kwargs):
        """
        required arguments:
        fq1, (fq2)
        genome
        outdir
        """
        self.update(kwargs, remove=True) # init self, fresh new
        self.status = self.init_args() # update all variables: *.config, *.args


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
        local_config = RNAseqConfig(**self.__dict__)
        self.update(local_config.__dict__, force=True)

        ## re-initiate args
        self.pickle = None
        self.build_design = False

        # update rnaseq_type
        rt = RNAseqConfig(**self.__dict__)
        self.update(rt.__dict__, force=True)

        ## update design
        self.design = os.path.join(self.config_dir, 'RNAseq_auto_design.txt')

        # check arguments
        chk1 = args_checker(self.__dict__, self.config_pickle, update=True)
        # Json(self.__dict__).writer(self.config_json)
        chk2 = args_logger(self.__dict__, self.config_txt, overwrite=True)
        # return all([chk1, chk2])


    def run(self):
        DesignBuilder(**self.__dict__).to_file(self.design)
        ## log
        lines = '\n' + '#'*104 + '\n'
        lines += '# {:<100s} #\n'.format('Organize aguments for RNAseq anslysis')
        lines += '# {:<100s} #\n'.format('1. arguments saved to plain text:')
        lines += '# {:<100s} #\n'.format(self.config_txt)
        lines += '# {:<100s} #\n'.format('2. arguments saved to pickle: ')
        lines += '# {:<100s} #\n'.format(self.config_pickle)
        lines += '# {:<100s} #\n'.format('')
        lines += '# {:<100s} #\n'.format('Run the following command for RNAseq analysis:')
        lines += '# {:<100s} #\n'.format('')
        lines += 'hiseq rnaseq --pickle {} \n'.format(self.config_pickle)
        # lines += '\n' + '#'*104 + '\n'
        log.info(lines)

