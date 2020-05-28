# -*- coding:utf-8 -*-

"""

Run GO, KEGG for specific gene list

"""


import os
import pathlib
import hiseq
import copy # copy objects
import pandas as pd
from hiseq.utils.helper import *


class Go(object):
    """
    Input: deseq_dir
    Input: genes, organism, foldChange
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
            'all': None,
            'input': None,
            'organism': None,
            'outdir': str(pathlib.Path.cwd()),
            'foldchange': None,
            'feature': 'gene',
            'ctl_vs_exp': '1'}
        self.update(args_default, force=False) # update missing attrs

        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        

    def is_enrich_dir(self):
        """
        contains: go_data.rds, go_plots.rds, go_enrich, go_group, go_gsea
        """
        req_dirs = ['go_data.rds', 'go_plots.rds', 'go_group', 'go_enrich']
        chk_dirs = [os.path.join(self.input, i) for i in req_dirs]
        return all(file_exists(chk_dirs))


    def get_sig_genes(self, x):
        """
        Extract sig genes, based on the column: sig
        """
        # read data
        df = pd.read_csv(x, '\t')
        # check
        if not 'sig' in df.columns:
            info.warning('column: sig, not found')
            raise Exception('argumnets failed')

        if not 'log2FoldChange' in df.columns:
            info.warning('column: log2FoldChange, not found')
            raise Exception('argumnets failed')

        # up
        up = df.sig.str.contains('up')
        down = df.sig.str.contains('down')
        up_down = up | down

        return [df[up], df[down], df[up_down]]


    def mission(self):
        """
        check input, and args
        """
        #if self.is_enrich_dir():
        if not self.all is None:
            if os.path.isfile(self.all):
                # for all sig genes
                return 'all'
            else:
                raise Exception('str, a file expected, {} got'.format(self.all))
        elif os.path.isfile(self.input):
            # require: outdir, organism
            if self.outdir is None:
                info.error('-o outdir required')
            if self.genome is None:
                info.error('-g genome required')
            # out
            if self.outdir is None or self.genome is None:
                raise Exception('arguments faild')

            return 'gene_list'
        elif os.path.isdir(self.input):
            return 'deseq_dir'
        else:
            info.warning('-i, str(file|dir) expected, {}  got'.format(self.input))
            raise Exception('-i arguments failed')


    def run(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        run_goR = os.path.join(pkg_dir, 'bin', 'run_go.R')
        cmd_file = os.path.join(self.outdir, 'cmd.sh')
        report_html = os.path.join(self.outdir, 'go_enrich_report.html')

        if self.mission() == 'all':
            # outdir
            is_path(self.outdir)
            # save sig genes
            sig_types = ['up', 'down', 'up_and_down']
            df_up, df_down, df_up_down = self.get_sig_genes(self.all)
            # up
            for i, n in enumerate(self.get_sig_genes(self.all)):
                log.info('Run GO analysis for: ' + sig_types[i])
                dir_sig = os.path.join(self.outdir, sig_types[i])
                file_sig = os.path.join(dir_sig, sig_types[i] + '.txt')
                is_path(dir_sig)
                # save genes to file
                n.to_csv(file_sig, '\t', header=True, index=False)
                # run GO
                Go(input=file_sig, outdir=dir_sig, genome=self.genome,
                    foldChange=self.all).run()
        elif self.mission() == 'gene_list':

            # outdir
            is_path(self.outdir)
            cmd = ' '.join(['Rscript', run_goR, self.input, self.genome, self.outdir])

            # foldChange
            if not self.foldChange is None:
                cmd += ' {}'.format(self.foldChange)

            # save cmd
            with open(cmd_file, 'wt') as w:
                w.write(cmd + '\n')
        else:
            cmd = ' '.join(['Rscript', run_goR, self.input, self.feature, self.ctl_vs_exp])

        if check_file(report_html):
            log.info('GO report() skipped, file exists: {}'.format(report_html))
        else:
            try:
                # run_shell_cmd(cmd)
                os.system(cmd)
            except:
                log.warning('report() failed.')





