
# -*- coding: utf-8 -*-

"""
Common functions for pipeline construction
file modification, 
...
"""

import os
import sys
import gzip
import shutil
import json
import pickle
import fnmatch
import logging
import functools
import subprocess
import pysam
import pybedtools
import pathlib
import binascii
from .args import args_init


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)


# Decorator:
class Logger(object):
    def __init__(self, level='INFO'):
        logging.basicConfig(
            format='[%(asctime)s %(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            stream=sys.stdout)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(level)

    def __call__(self, fn):
        @functools.wraps(fn)
        def decorated(*args, **kwargs):
            try:
                self.logger.info('{0} - {1} - {2}'.format(
                    fn.__name__, 
                    args, 
                    kwargs))
                result = fn(*args, **kwargs)
                self.logger.info(result)
                # return result
            except Exception as ex:
                self.logger.info('Exception {0}'.format(ex))
                raise ex
            return result
        return decorated


def is_gz(filepath):
    if os.path.exists(filepath):
        with open(filepath, 'rb') as test_f:
            return binascii.hexlify(test_f.read(2)) == b'1f8b'
    else:
        if filepath.endswith('.gz'):
            return True
        else:
            return False

def is_path(path, create = True):
    """
    Check path, whether a directory or not
    if not, create it
    """
    assert isinstance(path, str)
    if os.path.exists(path):
        return True
    else:
        if create:
            try:
                os.makedirs(path)
                return True
            except IOError:
                log.error('failed to create directories: %s' % path)
        else:
            return False


def run_shell_cmd(cmd):
    """This command is from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py
    """
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid, 
        pgid, 
        rc,
        stderr.strip(), 
        stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')


def file_prefix(fn, with_path=False):
    """
    extract the prefix of a file
    remove extensions
    .gz, .fq.gz
    """
    assert isinstance(fn, str)
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.endswith('gz') or px.endswith('.bz2'):
        px = os.path.splitext(p1)[1] + px
        p1 = os.path.splitext(p1)[0]
    if not with_path:
        p1 = os.path.basename(p1)
    return [p1, px]


def args_checker(d, x, update=False):
    """Check if dict and x are consitent
    d is dict
    x is pickle file
    """
    assert isinstance(d, dict)
    flag = None
    if os.path.exists(x):
        # read file to dict
        with open(x, 'rb') as fh:
            d_checker = pickle.load(fh)
        if d == d_checker:
            flag = True
        else:
            if update:
                with open(x, 'wb') as fo:
                    pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    elif isinstance(x, str):
        # save dict to new file
        with open(x, 'wb') as fo:
            pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        log.error('illegal x= argument: %s' % x)

    return flag


def args_logger(d, x, overwrite=False):
    """Format dict, save to file
        key: value
    """
    assert isinstance(d, dict)
    n = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
    if os.path.exists(x) and overwrite is False:
        return True
    else:
        with open(x, 'wt') as fo:
            fo.write('\n'.join(n) + '\n')
        return '\n'.join(n)


def gzip_cmd(src, dest, decompress=True, rm=True):
    """
    Gzip Compress or Decompress files using gzip module in python 
    rm, True/False, whether remove old file

    # check the src file by extension: .gz
    """
    if decompress:
        if is_gz(src):
            with gzip.open(src, 'rb') as r, open(dest, 'wb') as w:
                shutil.copyfileobj(r, w)
        else:
            log.warning('not a gzipped file: {}'.format(src))
            shutil.copy(src, dest)
    else:
        if is_gz(src):
            log.warning('input is gzipped file, no need gzip')
            shutil.copy(src, dest)
        else:
            with open(src, 'rb') as r, gzip.open(dest, 'wb') as w:
                shutil.copyfileobj(r, w)

    # output
    if rm is True:
        os.remove(f)

    return t


def listfiles(path, full_name=True, recursive=False, include_dir=False):
    """
    List all the files within the path
    """
    out = []
    for root, dirs, files in os.walk(path):
        if full_name:
            dirs = [os.path.join(root, d) for d in dirs]
            files = [os.path.join(root, f) for f in files]
        out += files

        if include_dir:
            out += dirs

        if recursive is False:
            break
    return out


def listfiles2(pattern, path='.', full_name=True, recursive=False):
    """
    List all the files in specific directory
    fnmatch.fnmatch()

    pattern:

    *       matches everything
    ?       matches any single character
    [seq]   matches any character in seq
    [!seq]  matches any char not in seq

    An initial period in FILENAME is not special.
    Both FILENAME and PATTERN are first case-normalized
    if the operating system requires it.
    If you don't want this, use fnmatchcase(FILENAME, PATTERN).

    example:
    listfiles('*.fq', './')
    """
    fn_list = listfiles(path, full_name, recursive, include_dir=False)
    fn_list = [f for f in fn_list if fnmatch.fnmatch(f, pattern)]
    return fn_list


class Json(object):

    def __init__(self, x):
        """
        x 
          - dict, save to file
          - json, save to file
          - file, read as dict
        Save dict to json file
        Read from json file as dict
        ...
        """
        self.x = x # input

        if isinstance(x, Json):
            self.dict = x.dict
        elif isinstance(x, dict):
            # input a dict, 
            # save to file
            self.dict = x
        elif os.path.exists(x):
            # a file saving json content
            self.dict = self.reader()
        else:
            raise Exception('unknown objec: {}'.format(x))

    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json', 
            delete=False)
        return tmp.name


    def reader(self):
        """
        Read json file as dict
        """
        if os.path.getsize(self.x) > 0:
            with open(self.x) as r:
                d = json.load(r)
        else:
            d = {}

        return d


    def writer(self, f=None):
        """
        Write d (dict) to file x, in json format
        """
        # save to file
        if f is None:
            f = self._tmp()

        assert isinstance(f, str)
        # assert os.path.isfile(f)

        if isinstance(self.x, dict):
            with open(f, 'wt') as w:
                json.dump(self.x, w, indent=4, sort_keys=True)

        return f


class Genome(object):
    """
    List related information of specific genome
    1. get_fa(), genome fasta
    2. get_fasize(), genome fasta size
    3. bowtie_index(), bowtie index, optional, rRNA=True
    4. bowtie2_index(), bowtie2 index, optional, rRNA=True
    5. star_index(), STAR index, optional, rRNA=True
    6. gene_bed(), 
    7. gene_rmsk(), 
    8. gene_gtf(), optional, version='ucsc|ensembl|ncbi'
    9. te_gtf(), optional, version='ucsc'
    10. te_consensus(), optional, fruitfly()
    ...

    directory structure of genome should be like this:
    /path-to-data/{genome}/
        |- bigZips  # genome fasta, fasize, chromosome 
        |- annotation_and_repeats  # gtf, bed, rRNA, tRNA, annotation
        |- bowtie_index
        |- bowtie2_index
        |- STAR_index
        |- hisat2_index
        |- phylop100
        |- ...

    default: $HOME/data/genome/{genome}

    """
    def __init__(self, genome, genome_path=None, 
        repeat_masked_genome=False, **kwargs):
        assert isinstance(genome, str)
        self.genome = genome
        self.repeat_masked_genome = repeat_masked_genome
        self.kwargs = kwargs

        # path
        if genome_path is None:
            genome_path = os.path.join(str(pathlib.Path.home()), 
                'data', 'genome')
        self.genome_path = genome_path

        # # deprecated
        # if not supportedGenome(genome):
        #     log.error('genome not supported: {}'.foramt(genome))


    def get_fa(self):
        """
        Get the fasta file of specific genome
        {genome}/bigZips/{genome}.fa
        also check ".gz" file
        """
        fa = os.path.join(self.genome_path, self.genome, 'bigZips', 
            self.genome + '.fa')
        if not os.path.exists(fa):
            # gencode version
            fa = os.path.join(self.genome_path, self.genome, 'fasta', 
                self.genome + '.fa')

        fa_gz = fa + '.gz'
        if not os.path.exists(fa):
            if os.path.exists(fa_gz):
                log.error('require to unzip the fasta file: %s' % fa_gz)
            else:
                log.error('fasta file not detected: %s' % fa)
            return None
        else:
            return fa


    def get_fasize(self):
        """Get the fasta size file, chromosome size
        optional, fetch chrom size from ucsc
        http://hgdownload.cse.ucsc.edu/goldenPath/<db>/bigZips/<db>.chrom.sizes

        or using UCSC tool: fetchChromSizes
        fetchChromSizes hg39 > hg38.chrom.sizes
        """
        fa = self.get_fa()
        fa_size = fa + '.chrom.sizes'

        if not os.path.exists(fa_size):
            # log.info('Downloading chrom.sizes from UCSC')
            # url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes' % self.genome
            # download(url, fa_size)
            log.warning('file not exists, run samtools faidx to generate it')
            pysam.faidx(fa) # create *.fa.fai
            os.rename(fa + '.fai', fa_size)

        return fa_size


    def phylop100(self):
        """Return the phylop100 bigWig file of hg19, only
        for conservation analysis
        """
        p = os.path.join(self.genome_path, self.genome, 'phyloP100way',
            self.genome + '.100way.phyloP100way.bw')
        if not os.path.exists(p):
            p = None
        return p


    def gene_bed(self, version='refseq', rmsk=False):
        """Return the gene annotation in BED format
        support UCSC, ensembl, gencode
        """
        if rmsk:
            suffix = '.rmsk.bed'
        else:
            suffix = '.refseq.bed'
        g = os.path.join(self.genome_path, self.genome, 
            'annotation_and_repeats', self.genome + suffix)
        if not os.path.exists(g):
            g = None
        return g


    def gene_gtf(self, version='refseq'):
        """Return the gene annotation in GTF format
        support refseq, ensembl, gencode
        """
        version = version.lower() #

        gtf = os.path.join(
            self.genome_path, 
            self.genome, 
            'annotation_and_repeats',
            self.genome + '.' + version + '.gtf')

        if not os.path.exists(gtf):
            gtf = os.path.join(
            self.genome_path, 
            self.genome, 
            'gtf',
            self.genome + '.' + version + '.gtf')

        if not os.path.exists(gtf):
            gtf = None

        return gtf


    def te_gtf(self, format='gtf'):
        """Return TE annotation of the genome
        or return TE consensus sequence for the genome (dm3)
        """
        # only dm3 supported
        te_gtf = os.path.join(self.genome_path, self.genome, 
            self.genome + '_transposon', 
            self.genome + '_transposon.gtf')
        if not os.path.exists(te_gtf):
            te_gtf = None

        return te_gtf


class Bam(object):
    """
    Manipulate BAM files
    - sort 
    - index 
    - merge
    - count 
    - to_bed
    - rmdup  
    - ...

    Using Pysam, Pybedtools, ...

    code from cgat: 
    """
    def __init__(self, infile, threads=4):
        self.bam = infile
        self.threads = threads
        # self.bed = self.to_bed()


    def index(self):
        """
        Create index for bam
        """
        bai = self.bam + '.bai'
        if not os.path.exists(bai):
            pysam.index(self.bam)

        return os.path.exists(bai)


    def sort(self, outfile=None, by_name=False, overwrite=False):
        """
        Sort bam file by position (default)
        save to *.sorted.bam (or specify the name)
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.sorted.bam'

        if os.path.exists(outfile) and overwrite is False:
            logging.info('file exists: {}'.format(outfile))
        else:
            tmp = pysam.sort('-@', str(self.threads), '-o', outfile, self.bam)

        return outfile


    def merge(self):
        """
        Merge multiple BAM files using samtools
        """
        
        # pysam.merge('')
        pass


    def count(self):
        """Using samtools view -c"""
        return pysam.view('-c', self.bam)


    def to_bed(self, outfile=None):
        """Convert BAM to BED
        pybetools 
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.bed'

        if not os.path.exists(outfile):
            pybedtools.BedTool(self.bam).bam_to_bed().saveas(outfile)

        return outfile


    def rmdup(self, outfile=None, overwrite=False):
        """
        Remove duplicates using picard/sambamba
        sambamba markdup -r --overflow-list-size 800000 raw.bam rmdup.bam
        picard MarkDuplicates  REMOVE_SEQUENCING_DUPLICATES=True I=in.bam O=outfile.bam
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.rmdup.bam'

        sambamba = shutil.which('sambamba')
        log = outfile + '.sambamba.log'
        cmd = '{} markdup -r -t {} --overflow-list-size 800000 \
            --tmpdir="./" {} {} 2> {}'.format(
            sambamba,
            str(self.threads),
            self.bam,
            outfile,
            log)

        if os.path.exists(outfile) and overwrite is False:
            logging.info('file exists: {}'.format(outfile))
        else:
            run_shell_cmd(cmd)

        return outfile


    def proper_pair(self, outfile=None, overwrite=False):
        """
        Extract proper pair
        samtools view -f 2
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.proper_pair.bam'

        if os.path.exists(outfile) and overwrite is False:
            logging.info('file exists: {}'.format(outfile))
        else:
            pysam.view('-f', '2', '-h', '-b', '-@', str(self.threads),
                '-o', outfile, self.bam, catch_stdout=False)

        return outfile


    ##########################################
    ## code from cgat
    ##########################################
    def isPaired(self, topn=1000):
        """
        Check if infile contains paired end reads
        
        go through the topn alignments in file, 

        return: True, any of the alignments are paired
        """
        samfile = pysam.AlignmentFile(self.bam)
        n = 0
        for read in samfile:
            if read.is_paired:
                break
            n += 1
            if n == topn:
                break

        samfile.close()

        return n != topn


    def getNumReads(self):
        """
        Count number of reads in bam file.

        This methods works through pysam.idxstats.

        Arguments
        ---------
        bamfile : string
            Filename of :term:`bam` formatted file. The file needs
            to be indexed.
        Returns
        -------
        nreads : int
            Number of reads
        """

        lines = pysam.idxstats(self.bam).splitlines()

        try:
            nreads = sum(
                map(int, [x.split("\t")[2]
                          for x in lines if not x.startswith("#")]))

        except IndexError as msg:
            raise IndexError(
                "can't get number of reads from bamfile, msg=%s, data=%s" %
                (msg, lines))
        return nreads
      

    def estimateInsertSizeDistribution(self, topn=10000, n=10, 
        method="picard", similarity_threshold=1.0, max_chunks=1000):
        """
        Estimate insert size from a subset of alignments in a bam file.

        Several methods are implemented.

        picard
            The method works analogous to picard by restricting the estimates
            to a core distribution. The core distribution is defined as all
            values that lie within n-times the median absolute deviation of
            the full data set.
        convergence
            The method works similar to ``picard``, but continues reading
            `alignments` until the mean and standard deviation stabilize.
            The values returned are the median mean and median standard
            deviation encountered.

        The method `convergence` is suited to RNA-seq data, as insert sizes
        fluctuate siginificantly depending on the current region
        being looked at.

        Only mapped and proper pairs are considered in the computation.

        Returns
        -------
        mean : float
           Mean of insert sizes.
        stddev : float
           Standard deviation of insert sizes.
        npairs : int
           Number of read pairs used for the estimation
        method : string
           Estimation method
        similarity_threshold : float
           Similarity threshold to apply.
        max_chunks : int
           Maximum number of chunks of size `alignments` to be used
           in the convergence method.

        """

        assert self.isPaired(self.bam), \
            'can only estimate insert size from' \
            'paired bam files'

        samfile = pysam.AlignmentFile(self.bam)

        def get_core_distribution(inserts, n):
            # compute median absolute deviation
            raw_median = numpy.median(inserts)
            raw_median_dev = numpy.median(numpy.absolute(inserts - raw_median))

            # set thresholds
            threshold_min = max(0, raw_median - n * raw_median_dev)
            threshold_max = raw_median + n * raw_median_dev

            # define core distribution
            return inserts[numpy.logical_and(inserts >= threshold_min,
                                             inserts <= threshold_max)]

        if method == "picard":

            # only get first read in pair to avoid double counting
            inserts = numpy.array(
                [read.template_length for read in samfile.head(n=topn)
                 if read.is_proper_pair
                 and not read.is_unmapped
                 and not read.mate_is_unmapped
                 and not read.is_read1
                 and not read.is_duplicate
                 and read.template_length > 0])
            core = get_core_distribution(inserts, n)

            return numpy.mean(core), numpy.std(core), len(inserts)

        elif method == "convergence":

            means, stds, counts = [], [], []
            last_mean = 0
            iteration = 0
            while iteration < max_chunks:

                inserts = numpy.array(
                    [read.template_length for read in samfile.head(
                        n=topn,
                        multiple_iterators=False)
                     if read.is_proper_pair
                     and not read.is_unmapped
                     and not read.mate_is_unmapped
                     and not read.is_read1
                     and not read.is_duplicate
                     and read.template_length > 0])
                core = get_core_distribution(inserts, n)
                means.append(numpy.mean(core))
                stds.append(numpy.std(core))
                counts.append(len(inserts))
                mean_core = get_core_distribution(numpy.array(means), 2)
                mm = numpy.mean(mean_core)
                if abs(mm - last_mean) < similarity_threshold:
                    break
                last_mean = mm

            return numpy.median(means), numpy.median(stds), sum(counts)
        else:
            raise ValueError("unknown method '%s'" % method)

  
    def estimateTagSize(self, topn=10, multiple="error"):
        """
        Estimate tag/read size from first alignments in file.

        Arguments
        ---------
        bamfile : string
           Filename of :term:`bam` formatted file
        alignments : int
           Number of alignments to inspect
        multiple : string
           How to deal if there are multiple tag sizes present.
           ``error`` will raise a warning, ``mean`` will return the
           mean of the read lengths found. ``uniq`` will return a
           unique list of read sizes found. ``all`` will return all
           read sizes encountered.

        Returns
        -------
        size : int
           The read size (actual, mean or list of read sizes)

        Raises
        ------
        ValueError
           If there are multiple tag sizes present and `multiple` is set to
           `error`.

        """
        samfile = pysam.AlignmentFile(self.bam)
        sizes = [read.rlen for read in samfile.head(topn)]
        mi, ma = min(sizes), max(sizes)

        if mi == 0 and ma == 0:
            sizes = [read.inferred_length for read in samfile.head(alignments)]
            # remove 0 sizes (unaligned reads?)
            sizes = [x for x in sizes if x > 0]
            mi, ma = min(sizes), max(sizes)

        if mi != ma:
            if multiple == "error":
                raise ValueError('multiple tag sizes in %s: %s' % (bamfile, sizes))
            elif multiple == "mean":
                mi = int(sum(sizes) / len(sizes))
            elif multiple == "uniq":
                mi = list(sorted(set(sizes)))
            elif multiple == "all":
                return sizes

        return mi


    def getNumberOfAlignments(self):
        """return number of alignments in bamfile.
        """
        samfile = pysam.AlignmentFile(self.bam)
        return samfile.mapped

    ##########################################
    ## code from cgat
    ##########################################



## for index
## to-do
##   - build index (not recommended)
##
class AlignIndex(object):

    def __init__(self, aligner='bowtie', index=None, **kwargs):
        """
        Required args:
          - aligner
          - index (optional)
          - genome
          - group : genome, rRNA, transposon, piRNA_cluster, ...
          - genome_path
        """
        ## init
        args = args_init(kwargs, align=True) # init

        self.aligner = aligner

        if isinstance(index, str):
            # index given
            self.index = index.rstrip('/') # 
            self.name = self.get_name()
            self.aligner_supported = self.get_aligner() # all
            self.check = index if self.is_index() else None
        else:
            # index not defined, search required
            # log.warning('index=, not defined; .search() required')
            pass
        
        self.kwargs = args        


    def get_aligner(self, index=None):
        """
        Search the available index for aligner:
        bowtie, [*.[1234].ebwt,  *.rev.[12].ebwt]
        bowtie2, [*.[1234].bt2, *.rev.[12].bt2]  
        STAR,
        bwa, 
        hisat2, 
        """
        if index is None:
            index = self.index

        bowtie_files = [index + i for i in [
            '.1.ebwt',
            '.2.ebwt',
            '.3.ebwt',
            '.4.ebwt',
            '.rev.1.ebwt',
            '.rev.2.ebwt']]

        bowtie2_files = [index + i for i in [
            '.1.bt2',
            '.2.bt2',
            '.3.bt2',
            '.4.bt2',
            '.rev.1.bt2',
            '.rev.2.bt2']]

        hisat2_files = ['{}.{}.ht2'.format(index, i) for i in range(1, 9)]


        bwa_files = [index + i for i in [
            '.sa',
            '.amb',
            '.ann',
            '.pac',
            '.bwt']]

        STAR_files = [os.path.join(index, i) for i in [
            'SAindex',
            'Genome',
            'SA',
            'chrLength.txt',
            'chrNameLength.txt',
            'chrName.txt',
            'chrStart.txt',
            'genomeParameters.txt']]

        ## check exists
        bowtie_chk = [os.path.exists(i) for i in bowtie_files]
        bowtie2_chk = [os.path.exists(i) for i in bowtie2_files]
        hisat2_chk = [os.path.exists(i) for i in hisat2_files]
        bwa_chk = [os.path.exists(i) for i in bwa_files]
        STAR_chk = [os.path.exists(i) for i in STAR_files]

        ## check file exists
        aligner = []

        if all(bowtie_chk):
            aligner.append('bowtie')
        elif all(bowtie2_chk):
            aligner.append('bowtie2')
        elif all(hisat2_chk):
            aligner.append('hisat2')
        elif all(bwa_chk):
            aligner.append('bwa')
        elif all(STAR_chk):
            aligner.append('STAR')
        else:
            pass

        return aligner


    def is_index(self, index=None):
        """
        Check if index support for aligner
        """
        if index is None:
            index = self.index
        return self.aligner in self.get_aligner(index)

        # return self.aligner in self.aligner_supported if 
        #    index is None else self.aligner in self.get_aligner(index)


    def get_name(self, index=None):
        """
        Get the name of index
        basename: bowtie, bowtie2, hisqt2, bwa
        folder: STAR
        """
        if index is None:
            index = self.index
        if os.path.isdir(index):
            # STAR
            iname = os.path.basename(index)
        else:
            # bowtie, bowtie2, bwa, hisat2
            iname = os.path.basename(index)

        return iname


    def search(self, genome=None, group='genome'):
        """
        Search the index for aligner: STAR, bowtie, bowtie2, bwa, hisat2
        para:

        *genome*    The ucsc name of the genome, dm3, dm6, mm9, mm10, hg19, hg38, ...
        *group*      Choose from: genome, rRNA, transposon, piRNA_cluster, ...

        structure of genome_path:
        default: {HOME}/data/genome/{genome_version}/{aligner}/

        path-to-genome/
            |- Bowtie_index /
                |- genome
                |- rRNA
                |- MT_trRNA
            |- transposon  
            |- piRNA cluster

        """
        args = self.kwargs.copy()

        g_path = args.get('genome_path', './')
        i_prefix = os.path.join(g_path, genome, self.aligner + '_index')

        # all index
        d = {
            'genome': [
                os.path.join(i_prefix, 'genome'),
                os.path.join(i_prefix, genome)],
            'genome_rm': [
                os.path.join(i_prefix, 'genome_rm'),
                os.path.join(i_prefix, genome + '_rm')],
            'MT_trRNA': [
                os.path.join(i_prefix, 'MT_trRNA')],
            'rRNA': [
                os.path.join(i_prefix, 'rRNA')],
            'chrM': [
                os.path.join(i_prefix, 'chrM')],
            'structural_RNA': [
                os.path.join(i_prefix, 'structural_RNA')],
            'te': [
                os.path.join(i_prefix, 'transposon')],
            'piRNA_cluster': [
                os.path.join(i_prefix, 'piRNA_cluster')],
            'miRNA': [
                os.path.join(i_prefix, 'miRNA')],
            'miRNA_hairpin': [
                os.path.join(i_prefix, 'miRNA_hairpin')]}

        # hit
        i_list = d.get(group, ['temp_temp'])
        i_list = [i for i in i_list if self.is_index(i)]

        return i_list[0] if len(i_list) > 0 else None





################################################################################
## functions

def sam_flag_check(query, subject):
    """
    Check two numbers, query (for filtering) in subject or not
    convert to binary mode
    q: 0000010  (2)
    s: 1011011  (91)
    q in s
    range: 0 - 2048 (SAM flag)
    """
    def to_bin(n):
        return '{0:012b}'.format(n)

    # convert to binary mode
    q = to_bin(eval(query))
    s = to_bin(eval(subject))

    # check q, s
    flag = True
    for j, k in zip(q[::-1], s[::-1]):
        if not j == '1':
            continue
        if eval(j) - eval(k) > 0:
            flag = False
            break

    return flag
    
