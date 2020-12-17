# compare two RNAseq data sets


import os
import sys
import hiseq
from hiseq.utils.helper import *


class DeseqPair(object):
    """
    Run RNAseq compare
    support for the results of `hiseq rnaseq` output
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True)
        self.init_args() # all args


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
        default arguments
        required:
          - fq1
          - genome/ext-index
          - outdir
        
        optional:
          - aligner

        prepare index
          - align-to-rRNA
          - align-to-chrM
        """
        args_default = {
            'dirA': None,
            'dirB': None,
            'feature': None,
            'outdir': str(pathlib.Path.cwd())}
        self.update(args_default, force=False) # update missing attrs

        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())

        # check path
        self.dirA = self.dirA.rstrip('/')
        self.dirB = self.dirB.rstrip('/')
        subdir = os.path.basename(self.dirA) + '.compare.' + os.path.basename(self.dirB)
        self.outdir = os.path.join(self.outdir, subdir) # add subdir
        
        # abs path
        self.dirA = os.path.abspath(self.dirA)
        self.dirB = os.path.abspath(self.dirB)


    def run(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        run_rnaseq_cmpR = os.path.join(pkg_dir, 'bin', 'run_deseq_pair.R')
        cmd_file = os.path.join(self.outdir, 'cmd.sh')
        cmd = ' '.join(['Rscript', 
            run_rnaseq_cmpR, 
            self.dirA, 
            self.dirB,
            self.outdir])

        # is_path(self.outdir)
        check_path(self.outdir)

        # save cmd
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')
        
        # run
        try:
            # run_shell_cmd(cmd)
            os.system(cmd)
        except:
            log.warning('rnaseq_cmp() failed.')