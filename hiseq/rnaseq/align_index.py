"""
name: align_index.py

Functions:
1. search index
2. validate index


"""

import os
import tempfile
from hiseq.utils.helper import * # all help functions



## for index
class AlignIndex(object):
    """
    arguments:
    index
    aligner
    group
    genome_path
    ...
    """
    def __init__(self, **kwargs):
        """
        Two keywords: index, aligner
        Required args:
          - aligner
          - index (optional)
          - genome
          - group : genome, rRNA, transposon, piRNA_cluster, ...
          - genome_path
        """
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'index': None,
            'aligner': None,
            'genome_path': os.path.join(str(pathlib.Path.home()), 'data', 'genome')
        }
        self = update_obj(self, args_init, force=False)
        # # update: remove `genome` from object
        # if hasattr(self, 'genome'):
        #     delattr(self, 'genome')


    def get_aligner(self, index=None):
        """
        Search the available index for aligner:
        bowtie, [*.[1234].ebwt,  *.rev.[12].ebwt]
        bowtie2, [*.[1234].bt2, *.rev.[12].bt2]
        STAR,
        bwa,
        hisat2,
        """
        # unknown
        if index is None:
            index = self.index

        if index is None: # required
            log.warning('AlignIndex(index=), required for guessing the aligner')
            return None

        # check
        bowtie_files = ['{}.{}.ebwt'.format(index, i) for i in range(1, 5)]
        bowtie2_files = ['{}.{}.bt2'.format(index, i) for i in range(1, 5)]
        hisat2_files = ['{}.{}.ht2'.format(index, i) for i in range(1, 4)]
        bwa_files = ['{}.{}'.format(index, i) for i in ['sa', 'amb', 'ann', 'pac', 'bwt']]
        star_files = [os.path.join(index, i) for i in [
            'SAindex',
            'Genome',
            'SA',
            'chrLength.txt',
            'chrNameLength.txt',
            'chrName.txt',
            'chrStart.txt',
            'genomeParameters.txt']]

        ## check
        chk0 = all(file_exists(bowtie_files))
        chk1 = all(file_exists(bowtie2_files))
        chk2 = all(file_exists(hisat2_files))
        chk3 = all(file_exists(bwa_files))
        chk4 = all(file_exists(star_files))

        ## check file exists
        if chk0:
            aligner = 'bowtie'
        elif chk1:
            aligner = 'bowtie2'
        elif chk2:
            aligner = 'hisat2'
        elif chk3:
            aligner = 'bwa'
        elif chk4:
            aligner = 'star' # STAR
        else:
            aligner = None

        return aligner


    def is_valid(self, index=None, aligner=None):
        """
        The index is valid, match aligner
        """
        if index is None:
            index = self.index

        if aligner is None:
            aligner = self.aligner
        
        # return the aligner, from index
        if aligner is None:
            log.warning('AlignIndex(aligner=), required')
            return False

        return aligner.lower() == self.get_aligner(index=index)


    def search(self):
        """
        Search the index for aligner: 
        STAR, bowtie, bowtie2, bwa, hisat2
        para:

        *genome*    The ucsc name of the genome, dm6, mm9, mm10, hg19, hg38, ...
        *group*      Choose from: genome, rRNA, transposon, piRNA_cluster, ...

        structure of genome_path:
        default: {HOME}/data/genome/{genome_version}/{aligner}/


        ## bowtie/bowtie2/hisat2/...
        path-to-genome/
            |- Bowtie_index /
                |- genome
                |- rRNA
                |- MT_trRNA
                |- transposon
                |- piRNA_cluster

        ## STAR
        path-to-genome/
            |- Bowtie_index /
                |- genome/
                |- rRNA/
                |- MT_trRNA/
                |- transposon/
                |- piRNA_cluster/
        """
        # required arguments: genome
        if not isinstance(self.genome, str):
            log.error('AlignIndex().search(), require, genome=.')
            return None

        # required arguments: genome
        if not isinstance(self.group, str):
            log.error('AlignIndex().search(), require, group=.')
            return None

        # required arguments: aligner
        if not isinstance(self.aligner, str):
            log.error('AlignIndex().search(), require, aligner=.')
            return None

        ## required arguments: group
        group_list = ['genome', 'genome_rm', 'MT_trRNA', 'rRNA', 'chrM', 
                      'structural_RNA', 'transposon', 'te', 'piRNA_cluster', 
                      'miRNA', 'miRNA_hairpin']
        if not self.group in group_list:
            log.error('AlignIndex().search(group={}) unknown, expect {}'.format(
                self.group, group_list))
            return None

        # require: aligner
        aligner_supported = ['bowtie', 'bowtie2', 'STAR', 'hisat2', 'bwa', 
                             'kallisto', 'salmon']
        if not self.aligner in aligner_supported:
            log.error('AlignIndex(aligner=) required, candidate: {}'.format(
                aligner_supported))
            return None

        ## create index path
        p0 = os.path.join(self.genome_path, self.genome, self.aligner + 
            '_index') # [case sensitive] STAR bowtie
        p1 = os.path.join(p0, self.group)

        return p1 if self.is_valid(p1, self.aligner) else None


    def index_name(self, index=None):
        """
        Get the name of index
        basename: bowtie, bowtie2, hisqt2, bwa
        folder: STAR
        """
        if index is None:
            index = self.index
        
        ## check
        if index is None:
            log.warning('AlignIndex(index=), required')
            return None

        if not self.is_valid(index=index):
            log_msg = '\n'.join([
                'index not exists, or not match the aligner:',
                '{:>30s}: {}'.format('Index', index),
                '{:>30s}: {}'.format('Aligner expected', self.get_aligner(index=index)),
                '{:>30s}: {}'.format('Aligner get', self.aligner)])
            log.warning(log_msg)
            return None

        if file_exists(index):
            return os.path.basename(index)
        elif os.path.basename(index) == 'genome':
            # ~/data/genome/dm3/bowtie2_index/genome
            # bowtie, bowtie2, bwa, hisat2
            # iname = os.path.basename(index)
            return os.path.basename(os.path.dirname(os.path.dirname(index)))
        else:
            return os.path.basename(index)


    def _tmp(self, is_dir=False, suffix='.txt'):
        """
        Create a tmp file to save json object
        """
        if is_dir:
            tmp = tempfile.TemporaryDirectory(prefix='tmp')
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
                delete=False)
        return tmp.name


    def index_size(self, index=None, return_file=False):
        """
        chr size of index

        bowtie:  bowtie-inspect -s <index>
        bowtie2: bowtie-inspect -s <index>
        STAR:  chrNameLength.txt
        """
        if index is None:
            index = self.index

        ## check
        if index is None:
            log.warning('AlignIndex(index=) or AlignIndex().index_name(index=) required')
            return None

        if not self.is_valid(index=index):
            log_msg = '\n'.join([
                'index not exists, or not match the aligner:',
                '{:>30s}: {}'.format('Index', index),
                '{:>30s}: {}'.format('Aligner expected', self.get_aligner(index=index)),
                '{:>30s}: {}'.format('Aligner get', self.aligner)])
            log.warning(log_msg)
            return None

        ## aligner
        gsize = self._tmp(suffix='.chrom.sizes')
        chrLength = 0
        aligner = self.get_aligner(index).lower()

        if aligner in ['bowtie', 'bowtie2', 'hisat2', 'star']:
            # get genome size
            if aligner.lower() == 'star':
                gsize = os.path.join(index, 'chrNameLength.txt')
            else:
                if aligner == 'bowtie':
                    x_inspect = shutil.which('bowtie-inspect')
                elif aligner == 'bowtie2':
                    x_inspect = shutil.which('bowtie2-inspect')
                elif aligner == 'hisat2':
                    x_inspect = shutil.which('hisat2-inspect')
                else:
                    pass

                # inspect
                cmd = ' '.join([
                    '{}'.format(x_inspect),
                    '-s {} |'.format(index),
                    'grep ^Sequence |',
                    "sed -E 's/^Sequence-[0-9]+\t//' > {}".format(gsize)])

                # run
                try:
                    os.system(cmd)
                except:
                    log.error('failed to run: {}'.format(x_inspect))

            # read size
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            chrLength = sum(map(int, s))

        else:
            log.error('unknown aligner: {}'.format(aligner))

        ## output
        return gsize if return_file else chrLength








