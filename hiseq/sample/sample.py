"""
FxSample

Generate sample file from fastx file

randome: seqkit sample

head:

tail:

"""


import os
import re
import glob
from multiprocessing import Pool
from hiseq.utils.helper import *
from hiseq.utils.seq import Fastx


class FxSample(object):
    """
    Get subset of fastx file
    random: true|false
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for CnR analysis
        """
        args_init = {
            'fx': None,
            'number': None,
            'outdir': None,
            'random': False,
            'overwrite': False,
            'parallel_jobs': 1 
        }
        self = update_obj(self, args_init, force=False)

        # file exists
        if isinstance(self.fx, str):
            self.fx = [self.fx]
        elif isinstance(self.fx, list):
            pass
        else:
            raise ValueError('-i expect str or list, got {}'.format(
                type(self.fx).__name__))

        # number
        if isinstance(self.number, int):
            self.number = abs(self.number)
        else:
            raise ValueError('-n expect int, got {}'.format(
                type(self.number).__name__))

        # output
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)


    def sampleR1(self, fx):
        """
        extract sample for single fq file
        """
        fx_out = os.path.join(self.outdir, os.path.basename(fx))
        if not fx_out.endswith('.gz'):
            fx_out += '.gz' # gzip output

        if file_exists(fx_out) and not self.overwrite:
            log.info('sample() skipped, file exists: {}'.format(fx_out))
        else:
            if self.random:
                log.info('extract {} random records from: {}'.format(self.number, fx))
                Fastx(fx).sample_random(fx_out, n=self.number)
            else:
                log.info('extract {} records from: {}'.format(self.number, fx))
                Fastx(fx).sample(self.outdir, sample_size=self.number, gzipped=True)


    def sampleRn(self):
        """
        Run for multiple fastq files
        """
        # run each fq
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.sampleR1, self.fx)





        # input = args.get('input', None)
        # outdir = args.get('outdir', None) # os.getcwd())
        # nsize = args.get('sample_size', 100)
        # if outdir is None:
        #     outdir = os.getcwd()

        # # get list
        # def get_fq(x):
        #     if os.path.isdir(x):
        #         # f_list = list_fq_files(x)
        #         p = re.compile('\.f(ast)?(a|q)(.gz)?')
        #         f_list = [i for i in listfile(x, "*") if p.search(i)]
        #     elif os.path.isfile(x):
        #         f_list = [x]
        #     else:
        #         raise ValueError('str,list expected, {} got'.format(type(input).__name__))

        #     return f_list

        # # output
        # fq_list = []
        # if isinstance(input, list):
        #     for i in input:
        #         fq_list.extend(get_fq(i))
        # elif isinstance(input, str):
        #     fq_list.extend(get_fq(input))
        # else:
        #     raise ValueError('str,list expected, {} got'.format(type(input).__name__))

        # # run 
        # for i in fq_list:
        #     log.info('sample fastq: {} {}'.format(i, nsize))
        #     Fastx(i).sample(outdir, nsize)


    def run(self):
        self.sampleRn()


