


"""
Run RNAseq pipe 


required: design, outdir, genome, threads

genome
transposon
piRNA_cluster
...

"""

import hiseq
import copy # copy objects
from hiseq.utils.helper import *
from hiseq.align.alignment import AlignIndex
from hiseq.rnaseq.rnaseq import RNAseq


class RNAseqPipe(object):
    def __init__(self, **kwargs):
        self.update(kwargs, force=True) # init
        self.args_init()


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


    def args_init(self):
        """
        Default args for RNAseq()
        """
        args_init = {
            'outdir': str(pathlib.Path.cwd()),
            'aligner': 'STAR', # default
            'align_to_rRNA': True,
            'extra_index': None,
            'extra_para': None,
            'genomeLoad': 'LoadAndRemove',
            'read1_only': False,
            'parallel_jobs': 4,
            'overwrite': False,
            'read1_only': True,
            'spikein': None,
            'trimmed': True,
            'unique_only': True,
            'smp_name': None,
            'mode': 'gtp'}
        self.update(args_init, force=False)

        # required args: design, genome, outdir, threads
        genome = getattr(self, 'genome', None)
        if genome is None:
            raise ValueError('-g required')

        # required args: design, genome, outdir, threads
        design = getattr(self, 'design', None)
        if design is None:
            raise ValueError('--design required')


    def run_genome(self):
        """
        require: genome
        remove: index_list
        """
        args_local = self.__dict__.copy()

        ## default
        args_init = {
            'aligner': 'STAR',
            'feature': 'gene',
            'index_list': None  # remove index_list, if exists
        }
        args_local.update(args_init)

        RNAseq(**args_local).run()


    def run_te(self):
        """
        require: index_list, gtf, feature: te
        """
        args_local = self.__dict__.copy()
        genome = args_local.get('genome', None)

        # get index
        te_index = AlignIndex(aligner='STAR').search(genome=genome, group='te')
        te_gtf = Genome(genome).te('gtf')

        ## default
        args_init = {
            'aligner': 'STAR',
            'feature': 'te',
            'unique_only': True,
            'index_list': te_index,
            'gtf': te_gtf
        }
        args_local.update(args_init)

        if te_index is None or not os.path.exists(te_gtf):
            log.warning('te not found, skipped ...')
        else:
            RNAseq(**args_local).run()


    def run_piRNA_cluster(self):
        """
        require: index_list, gtf, feature: piRNAcluster
        """
        args_local = self.__dict__.copy()
        genome = args_local.get('genome', None)

        # get index
        pi_index = AlignIndex(aligner='STAR').search(genome=genome, group='piRNA_cluster')
        pi_gtf = Genome(genome).piRNA_cluster('gtf')

        ## default
        args_init = {
            'aligner': 'STAR',
            'unique_only': True,
            'feature': 'piRNA_cluster', 
            'index_list': pi_index,
            'gtf': pi_gtf
        }
        args_local.update(args_init)

        if pi_index is None or not os.path.exists(pi_gtf):
            log.warning('piRNA_cluster not found, skipped ...')
        else:
            RNAseq(**args_local).run()


    def run(self):
        if('g' in self.mode):
            self.run_genome()
        if('t' in self.mode):
            self.run_te()
        if('p' in self.mode):
            self.run_piRNA_cluster()



