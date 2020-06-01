


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
        print('!AAAA-3', self.outdir)
        self.args_init()
        print('!AAAA-4', self.outdir)


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
            'smp_name': None}
        self.update(args_init, force=False)

        # require design
        design = getattr(self, 'design', None)
        if design is None:
            raise ValueError('--design required')


    def run_genome(self):
        """
        require: genome
        remove: index_list
        """
        self_local = copy.copy(self)
        setattr(self_local, 'aligner', 'STAR')
        genome = getattr(self_local, 'genome', None)

        if genome is None:
            raise ValueError('-g required')

        # remove Index_list
        self_local.index_list = None
        setattr(self_local, 'feature', 'te')
        RNAseq(**self_local.__dict__).run()


    def run_te(self):
        """
        require: index_list, gtf, feature: te
        """
        self_local = copy.copy(self)
        setattr(self_local, 'aligner', 'bowtie2')
        genome = getattr(self_local, 'genome', None)
        aligner = getattr(self_local, 'aligner', 'bowtie2')

        if genome is None:
            raise ValueError('-g required')

        # get index
        te_index = AlignIndex(aligner=aligner).search(genome=genome, group='te')
        te_gtf = Genome(genome).te()

        if te_index is None:
            log.info('te index not found, skipped...')
        else:
            setattr(self_local, 'index_list', te_index)
            setattr(self_local, 'gtf', te_gtf)
            setattr(self_local, 'feature', 'te')
            RNAseq(**self_local.__dict__).run()


    def run_piRNA_cluster(self):
        """
        require: index_list, gtf, feature: piRNAcluster
        """
        self_local = copy.copy(self)
        setattr(self_local, 'aligner', 'bowtie2')
        genome = getattr(self_local, 'genome', None)
        aligner = getattr(self_local, 'aligner', 'bowtie2')

        if genome is None:
            raise ValueError('-g required')

        # get index
        pi_index = AlignIndex(aligner=aligner).search(genome=genome, group='piRNA_cluster')
        pi_gtf = Genome(genome).piRNA_cluster('gtf')

        if pi_index is None:
            log.info('te index not found, skipped...')
        else:
            setattr(self_local, 'index_list', pi_index)
            setattr(self_local, 'gtf', pi_gtf)
            setattr(self_local, 'feature', 'piRNA_cluster')
            RNAseq(**self_local.__dict__).run()


    def run(self):
        self.run_genome()
        self.run_te()
        self.run_piRNA_cluster()



