#!/usr/bin/env python3

"""
Manage files for Genome

files:
  - fasta 
  - fa_size
  - index [bowtie, bowtie2, star, bwa, ...]
  - bed 
  - gtf
  - te
  - piRNA cluster
  - blacklist
  
check:
  - supported (hiseq.utils.is_supported)

"""

import os
import shutil
import pathlib
import pysam
from xopen import xopen
from hiseq.utils.utils import log, download_file, is_supported, update_obj
from hiseq.utils.file import file_exists, check_path
import zipfile


class Genome(object):
    """
    List related information of specific genome
    1. fa(), genome fasta
    2. fasize(), genome fasta size
    3. gene_bed(),
    4. gene_rmsk(),
    5. gene_gtf(), optional, version='ucsc|ensembl|ncbi'
    6. te_gtf(), optional, version='ucsc'
    7. te_consensus(), optional, fruitfly()
    ...

    directory structure of genome should be like this:
    /path-to-data/{genome}/
        |- bigZips  # genome fasta, fasize, chromosome
        |- annotation_and_repeats  # gtf, bed, rRNA, tRNA, annotation
        |- phylop100
        |- ...

    default: $HOME/data/genome/{genome}
    """
    def __init__(self, genome, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.genome = genome
        self.init_args()
        
        
    def init_args(self):
        args_init = {
            'genome_path': None,
            'repeat_masked_genome': False,
            'annotation_dir': 'annotation_and_repeats'
        }
        self = update_obj(self, args_init, force=False)
        if not isinstance(self.genome_path, str):
            self.genome_path = os.path.join(str(pathlib.Path.home()),
                'data', 'genome')
        self.flag = is_supported(self.genome)
        if not self.flag:
            log.error('unknown genome: {}'.format(self.genome))


    def fa(self):
        """
        Get the fasta file of specific genome
        {genome}/bigZips/{genome}.fa
        also check ".gz" file
        """
        if not self.flag:
            return None
        fa = os.path.join(self.genome_path, self.genome, 'bigZips',
            self.genome+'.fa')
        if not file_exists(fa):
            fa = os.path.join(self.genome_path, self.genome, 'fasta',
                self.genome+'.fa')
        fa_gz = fa+'.gz'
        if file_exists(fa):
            out = fa
        else:
            if file_exists(fa_gz):
                with xopen(fa_gz, 'rb') as r:
                    with xopen(fa, 'wb') as w:
                        shutil.copyfileobj(r, w)
            else:
                log.error('fasta file not detected: {}'.format(fa))
                out = None
        return out


    def fasize(self):
        """Get the fasta size file, chromosome size
        optional, fetch chrom size from ucsc
        http://hgdownload.cse.ucsc.edu/goldenPath/<db>/bigZips/<db>.chrom.sizes

        or using UCSC tool: fetchChromSizes
        fetchChromSizes hg39 > hg38.chrom.sizes
        """
        if not self.flag:
            return None
        fa = self.fa()
        fa_size = fa + '.chrom.sizes'
        if not file_exists(fa_size):
            log.warning('file not exists, run samtools faidx to generate it')
            pysam.faidx(fa) # create *.fa.fai
            os.rename(fa + '.fai', fa_size)
        out = fa_size
        return out


    def phylop100(self):
        """Return the phylop100 bigWig file of hg19, only
        for conservation analysis
        """
        if not self.flag:
            return None
        out = os.path.join(
            self.genome_path, self.genome, 'phyloP100way',
            self.genome + '.100way.phyloP100way.bw')
        if not file_exists(out):
            out = None
        return out


    def gene_bed(self, version='ensembl', rmsk=False):
        """
        Return the gene annotation in BED format
        support UCSC, ensembl, gencode
        """
        if not self.flag:
            return None
        version = version.lower() #
        if rmsk:
            suffix = '.rmsk.bed'
        else:
            suffix = '.{}.bed'.format(version)
        out = os.path.join(self.genome_path, self.genome,
            self.annotation_dir, self.genome+suffix)
        if not file_exists(out):
            out = None
        return out


    def gene_gtf(self, version='refseq'):
        """
        Return the gene annotation in GTF format
        support refseq, ensembl, gencode
        """
        if not self.flag:
            return None
        version = version.lower() #
        suffix = '.{}.gtf'.format(version)
        out = os.path.join(self.genome_path, self.genome,
            self.annotation_dir, self.genome+suffix)
        if not file_exists(out):
            out = None
        return out


    def te(self, format='gtf'):
        """
        Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        if not self.flag:
            return None
        out = os.path.join(self.genome_path, self.genome,
            self.genome + '_transposon',
            self.genome + '_transposon.' + format)
        if not file_exists(out):
            out = None
        return out


    def piRNA_cluster(self, format='gtf'):
        """
        Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        if not self.flag:
            return None
        out = os.path.join(self.genome_path, self.genome,
            self.genome + '_piRNA_clusters',
            self.genome + '_piRNA_clusters.' + format)
        if not file_exists(out):
            out = None
        return out
    
    
    def blacklist(self):
        """
        blacklist files were downloaded from github repo:
        https://github.com/Boyle-Lab/Blacklist/
        supported:
        dm3, dm6, ce10, ce11, mm10, hg19, and hg38 
        
        the local path is:
        genome_path/{genome}/annotation_and_repeats/blacklist/{genome}.blacklist.v2.bed
        """
        if not self.flag:
            return None
        # update url
        remote_url = 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists'
        remote_filename = {
            'dm3': 'dm3-blacklist.v2.bed.gz',
            'dm6': 'dm6-blacklist.v2.bed.gz',
            'ce10': 'ce10-blacklist.v2.bed.gz',
            'cd11': 'cd11-blacklist.v2.bed.gz',
            'mm10': 'mm10-blacklist.v2.bed.gz',
            'hg19': 'hg19-blacklist.v2.bed.gz',
            'hg38': 'hg38-blacklist.v2.bed.gz'
        }
        remote_file = {
            k:os.path.join(remote_url, v) for k,v in remote_filename.items()
        }
        local_filename = {k:os.path.splitext(v)[0] for k,v in remote_filename.items()}
        # update path
        local_file = {
            k:os.path.join(self.genome_path, self.genome, self.annotation_dir, 
                           'blacklist', v) for k,v in local_filename.items()
        }
        # check available
        if self.genome in local_file:
            local_bl = local_file.get(self.genome, None)
            local_dir = os.path.dirname(local_bl)
            remote_bl = remote_file.get(self.genome, None)
            if not file_exists(local_bl):
                try:
                    log.info('downloading blacklist file from github [Boyle-Lab/Blacklist]')
                    check_path(local_dir)
                    local_tmp = os.path.join(local_dir, os.path.basename(remote_bl))
                    # download
                    download_file(remote_bl, local_tmp)
                    # gunzip
                    with xopen(local_tmp, 'rb') as r:
                        with xopen(local_bl, 'wb') as w:
                            shutil.copyfileobj(r, w)
                except:
                    msg = '\n'.join([
                        'failed downloading blacklist',
                        'Using the following command in terminal to download file manually:',
                        '$ mkdir -p {}'.format(os.path.dirname(local_bl)),
                        '$ wget -O {} {}'.format(local_bl+'.gz', remote_bl),
                        '$ gunzip {}'.format(local_bl+'.zip')
                    ])
                    log.error(msg)
        else:
            log.error('no blacklist files found for {}'.format(self.genome))
            local_bl = None
        return local_bl