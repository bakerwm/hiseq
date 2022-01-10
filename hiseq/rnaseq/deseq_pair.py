# compare two RNAseq data sets


import os
import sys
import pathlib
import hiseq
from hiseq.utils.utils import (
    log, update_obj, Config, get_date, init_cpu, read_hiseq, is_supported, 
    print_dict
)
# from hiseq.utils.helper import *


class DeseqPair(object):
    """
    Run RNAseq compare
    support for the results of `hiseq rnaseq` output
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args() # all args


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
        self = update_obj(self, args_default, force=False)
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