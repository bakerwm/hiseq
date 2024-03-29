"""Quality control for HiSeq data

- Correlation: 500 window, pearson correlation  
- peaks, overlap  
- IDR, peaks (ChIP-seq, ATAC-seq)  
...

"""

import os
import sys
import re
import shlex
import subprocess
from xopen import xopen
from shutil import which
from PIL import Image # pillow
import pysam
import pybedtools
from collections import OrderedDict
from hiseq.utils.helper import *
from hiseq.utils.seq import Fastx
import deeptools.countReadsPerBin as crpb


def update_obj(obj, d, force=True, remove=False):
    """
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    # fresh start
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)

    return obj


def fasize(x, sum_all=False):
    """
    Count number of bases for each sequence
    """
    x_type = Fastx(x).format # 'fasta', 'fastq'

    k = {}
    with xopen(x) as r:
        for n, s, _ in Fastx(x).readfq(r):
            k[n] = len(s)

    return sum(k.values()) if sum_all else k


def symlink(src, dest, absolute_path=True):
    """
    Create symlinks within output dir
    ../src
    """
    if src is None or dest is None:
        log.warning('symlink skipped: {}, to: {}'.format(src, dest))
    elif file_exists(dest):
        log.warning('symlink skipped, target exists...'.format(dest))
    else:
        if absolute_path:
            # support: ~, $HOME,
            srcname = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        else:
            # only for directories within the same folder
            srcname = os.path.join('..', os.path.basename(src))

        # if not os.path.exists(dest):
        try:
            if not file_exists(dest):
                os.symlink(srcname, dest)
        except:
            log.error('symlink() failed')


def in_dict(d, k):
    """
    Check the keys in dict or not
    """
    assert isinstance(d, dict)
    if isinstance(k, str):
        k_list = [k]
    elif isinstance(k, list):
        k_list = list(map(str, k))
    else:
        log.warning('expect str and list, not {}'.format(type(k)))
        return False

    return all([i in d for i in k_list])


def in_attr(x, a, return_values=True):
    """
    Check a (attributes) in object a or not
    return the values or not
    """
    if isinstance(a, str):
        a_list = [a]
    elif isinstance(a, list):
        a_list = list(map(str, a))
    else:
        log.warning('expect str and list, not {}'.format(type(a)))
        return False

    # status
    status = all([hasattr(x, i) for i in a_list])

    if status and return_values:
        # values
        return [getattr(x, i) for i in a_list]
    else:
        return status


def check_file(x, show_log=False):
    """
    Check if x file, exists
    """
    if isinstance(x, str):
        x_list = [x]
    elif isinstance(x, list):
        x_list = x
    else:
        log.warning('expect str and list, not {}'.format(type(x)))
        return False        

    if show_log:
        for i in x_list:
            flag = 'ok' if os.path.exists(i) else 'fail'
            print('{:6s} : {}'.format(flag, i))

    return all(map(os.path.exists, x_list))


def check_path(x):
    """
    Check if x is path, Create path
    """
    if isinstance(x, str):
        x_list = [x]
    elif isinstance(x, list):
        x_list = list(map(str, x))
    else:
        log.warning('expect str and list, not {}'.format(type(x)))
        return False

    return all(map(is_path, x_list))


def merge_names(x):
    """
    Get the name of replictes
    common in left-most
    """
    assert isinstance(x, list)
    name_list = [os.path.basename(i) for i in x]
    name_list = [re.sub('.rep[0-9].*$', '', i) for i in name_list]
    return list(set(name_list))[0]


def dict_to_log(d, x, overwrite=False):
    """
    Convert dict to log style
        key | value
    """
    assert isinstance(d, dict)
    logout = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
    if overwrite is True or not os.path.exists(x): 
        with open(x, 'wt') as w:
            w.write('\n'.join(logout) + '\n')

    return '\n'.join(logout)


def dict_to_pickle(d, x):
    """
    Convert dict to pickle
    """
    assert isinstance(d, dict)
    with open(x, 'wb') as w:
        pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_to_dict(x):
    """
    Convert pickle file to dict
    """
    with open(x, 'rb') as r:
        return pickle.load(r)


def bed2gtf(infile, outfile):
    """Convert BED to GTF
    chrom chromStart chromEnd name score strand
    """
    with open(infile) as r, open(outfile, 'wt') as w:
        for line in r:
            fields = line.strip().split('\t')
            start = int(fields[1]) + 1
            w.write('\t'.join([
                fields[0],
                'BED_file',
                'gene',
                str(start),
                fields[2],
                '.',
                fields[5],
                '.',
                'gene_id "{}"; gene_name "{}"'.format(fields[3], fields[3])
                ]) + '\n')
    return outfile


def bam2bigwig(bam, genome, outdir, strandness=0, binsize=1, overwrite=False, **kwargs):
    """Convert bam to bigWig using deeptools
    https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    history:
    1. Mappable sequence of a genome, see Table 1 in 
       url: https://www.nature.com/articles/nbt.1518.pdf
    2. effective genome size:
        - non-N bases
        - regions (of some size) uniquely mappable
    3. UCSC
    http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics
    http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics

    !!! strandness
    default: dUTP-based library (read2 is sense strand, read1 is anti-sense strand)
    general RNA library: (NSR, small-RNA-library), read1 is sense, read2 is antisense
    """
    assert os.path.exists(bam)
    assert isinstance(genome, str)
    assert is_path(outdir)
    assert isinstance(strandness, int)
    assert isinstance(binsize, int)
    assert isinstance(overwrite, bool)

    bamcoverage = which('bamCoverage')
    if bamcoverage is None:
        raise ValueError('failed | program not found: bamCoverage')

    # genome size
    effsize = {'dm3': 162367812,
               'dm6': 142573017,
               'mm9': 2620345972,
               'mm10': 2652783500,
               'hg19': 2451960000,
               'hg38': 2913022398,
               'GRCh38': 2913022398}
    genome_size = effsize.get(genome, 0) 

    # create bam index
    Bam(bam).index()

    # prepare output
    prefix = fq_name(bam)
    bw_log = os.path.join(outdir, prefix + '.deeptools.log')
    # log.info('create bigWig for: %s' % prefix)

    
    def run_run(cmd, bw, bw_log):
        if os.path.exists(bw) and not overwrite:
            log.info('file exists, bigWig skipped ...')
        else:
            with open(bw_log, 'wt') as w:
                p1 = subprocess.run(shlex.split(cmd), stdout=w, stderr=w)
        if not os.path.exists(bw):
            raise ValueError('output file is missing, check log file: %s' % bw_log)


    def run_fwd():
        bw_fwd = os.path.join(outdir, prefix + '.fwd.bigWig')
        bw_log = os.path.join(outdir, prefix + '.fwd.deeptools.log')
        cmd = ' '.join([
                '{} -b {} -o {}'.format(bamcoverage, bam, bw_fwd),
                '--filterRNAstrand forward',
                '--binSize {} --effectiveGenomeSize {}'.format(binsize, genome_size)])
        run_run(cmd, bw_fwd, bw_log)


    def run_rev():
        bw_fwd = os.path.join(outdir, prefix + '.rev.bigWig')
        bw_log = os.path.join(outdir, prefix + '.rev.deeptools.log')
        cmd = ' '.join([
                '{} -b {} -o {}'.format(bamcoverage, bam, bw_rev),
                '--filterRNAstrand reverse',
                '--binSize {} --effectiveGenomeSize {}'.format(binsize, genome_size)])
        run_run(cmd, bw_rev, bw_log)


    def run_non():
        # non strandness
        bw = os.path.join(outdir, prefix + '.bigWig')
        bw_log = os.path.join(outdir, prefix + '.deeptools.log')
        cmd = ' '.join([
                '{} -b {} -o {}'.format(bamcoverage, bam, bw_rev),
                '--binSize {} --effectiveGenomeSize {}'.format(binsize, genome_size)])
        run_run(cmd, bw, bw_log)


    if strandness == 0:
        # non-strandness
        run_non()
    elif strandness == 12:
        # fwd + rev
        run_fwd()
        run_rev()
    elif strandness == 1:
        # fwd
        run_fwd()
    elif strandness == 1:
        # rev
        run_rev()
    else:
        log.error('unkown strandness: [0|1|2|12], {} got'.format(strandness))


def peak_FRiP(inbed, inbam, threads=4):
    """Calculate FRiP for ChIP-seq/ATAC-seq peaks
    Fraction of reads in called peak regions
    input: bam, peak (bed)

    using: bedtools intersect 
    count reads in peaks

    Using: Deeptools
    """
    # init data
    frip = 0
    rip = 0
    total = 0

    if file_exists(inbed):
        if file_row_counter(inbed) > 10:
            # reads in peak
            cr = crpb.CountReadsPerBin([inbam], bedFile=inbed, numberOfProcessors=threads)
            crn = cr.run()
            rip = crn.sum(axis=0)[0]
            # fraction
            bam = pysam.AlignmentFile(inbam)
            total = bam.mapped
            frip = float(rip) / total

    # output
    return (frip, rip, total)


def frag_length(infile, outfile=None):
    """
    Calculate Fragment length from BAM file
    BAM: column-9
    """
    d = {}
    for read in pysam.AlignmentFile(infile):
        if not read.is_proper_pair:
            continue
        if not read.is_read1:
            continue
        size = abs(read.tlen)
        d[size] = d.get(size, 0) + 1
    d = sorted(d.items(), key=lambda x: x[0])

    if not outfile is None:
        with open(outfile, 'wt') as w:
            for k, v in d:
                w.write(str(k) + '\t' + str(v) + '\n')

    return d


def peak_overlap(peak_list, outdir):
    """
    pybedtools.contrib.venn_maker.venn_maker(beds, names, figure_filename, script_filename, run=True)

    make venn plot
    """
    assert is_path(outdir)

    names = [os.path.splitext(i)[0] for i in peak_list]

    # name
    pname = os.path.basename(os.path.splitext(peak_list[0])[0])
    pname = pname.rstrip('_peaks')
    pname = pname.rstrip('_rep1')

    # output
    prefix = os.path.join(outdir, pname)
    vennR = prefix + '.venn.R'
    tiffout = prefix + '.venn.tiff'

    # func
    plt = pybedtools.contrib.venn_maker.venn_maker
    plt(peak_list, names, figure_filename=tiffout, script_filename=vennR, run=True)


def peak_overlap2(peak_list, outdir):
    """
    Calculate the overlap between bedA and bedB
    nameA, nameB, countA, countB, overlap, pctA, pctB
    make venn plot
    """
    if not len(peak_list) == 2:
        raise ValueError('expect two BED files, but {} received.'.format(len(peak_list)))

    assert is_path(outdir)

    # name
    pname = os.path.splitext(peak_list[0])[0]
    pname = pname.rstrip('_peaks')
    pname = pname.rstrip('_rep1')

    # output
    prefix = os.path.join(outdir, pname)
    tiffout = prefix + '.overlap.png'

    # number of overlaps
    a = pybedtools.BedTool(peak_list[0])
    b = pybedtools.BedTool(peak_list[1])

    a_in_b = (a+b).count()
    a_not_b = (a-b).count()
    b_not_a = (b-a).count()

    names = [os.path.splitext(i)[0] for i in peak_list]

    # out
    rpt = names + [a.count(), 
    	b.count(), 
    	a_in_b, 
        '{:.2f}'.format(a_in_b/a.count()*100),
        '{:.2f}'.format(a_in_b/b.count()*100)]
    rpt = map(str, rpt) # convert to string

    peak_overlap(peak_list, outdir)

    return list(rpt)


def peak_idr(peak_list, outdir):
    """
    Evaluate the IDR for peaks of MACS2
    narrowPeak files, require sorted by -log10(p-value) (column-8)
    ...
    """
    # filenames
    pname = os.path.splitext(peak_list[0])[0]
    pname = pname.rstrip('_peaks')
    pname = pname.rstrip('_rep1')

    # output
    prefix = os.path.join(outdir, pname)
    txt = prefix + '.idr.txt'
    png = prefix + '.idr.txt.png'
    log = prefix + '.idr.log'
    is_path(outdir) # create dir

    # sort narrowPeak by pvalue
    def sortPval(infile, outfile):
        cmd1 = 'sort -k8,8nr -o {} {}'.format(
            outfile, infile)
        run_shell_cmd(cmd1)

    peakA = os.path.join(outdir, os.path.basename(peak_list[0]))
    peakB = os.path.join(outdir, os.path.basename(peak_list[1]))
    sortPval(peak_list[0], peakA)
    sortPval(peak_list[1], peakB)

    # cmd
    idr = which('idr')
    if idr is None:
        raise Exception('Command not found: idr')

    cmd = '{} --samples {} {} --input-file-type narrowPeak \
        --rank p.value --output-file {} --plot \
        --log-output-file {}'.format(
            idr,
            peakA,
            peakB,
            txt,
            log)
    if os.path.exists(txt) and not overwrite is True:
        log('file exists, skip idr : {}'.format(pname))
    else:
        run_shell_cmd(cmd)

    return txt


class Bam2bw(object):
    """
    Convert bam to bigWig
    version-1: bedtools genomecov -bdg , bedGraphToBigWig, ...
    version-2: deeptools, genomeCoverage

    Convert bam to bigWig using deeptools
    https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
    history:
    1. Mappable sequence of a genome, see Table 1 in 
       url: https://www.nature.com/articles/nbt.1518.pdf
    2. effective genome size:
        - non-N bases
        - regions (of some size) uniquely mappable
    3. UCSC
    http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics
    http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics

    !!! strandness
    default: dUTP-based library (read2 is sense strand, read1 is anti-sense strand)
    general RNA library: (NSR, small-RNA-library), read1 is sense, read2 is antisense
    """
    def __init__(self, **kwargs):
        """
        bam, file to bam file
        outdir, str, path to save all files
        genome, dm3,dm6,...
        binsize, int, default: 50
        strandness, 0, [0|1|2|12], 1=fwd, 2=rev, 12=rev,rev, 0=non-strandness
        reference, str, fasta file
        overwrite, bool, False

        Required:
        bam, outdir, genome, binsize, strandness, overwrite, reference
        """
        # update self (obj)
        for k, v in kwargs.items():
            setattr(self, k, v)

        # update, defaults
        self.init_args()
        args_logger(self.__dict__, self.config_txt)


    def init_args(self):
        # genome size
        effsize = {
            'dm3': 162367812,
            'dm6': 142573017,
            'mm9': 2620345972,
            'mm10': 2652783500,
            'hg19': 2451960000,
            'hg38': 2913022398,
            'GRCh38': 2913022398}
        # gsize = effsize.get(genome, 0) 

        args_init = {
            'bam': None,
            'outdir': str(pathlib.Path.cwd()),
            'binsize': 50,
            'strandness': 0,
            'overwrite': False,
            'reference': None,
            'genome': None,
            'genome_size': None,
            'scaleFactor': 1.0,
            'normalizeUsing': 'RPKM',
            'threads': 4,
            'config_txt': os.path.join(self.outdir, 'arguments.txt'),
            'flag': True # whether run/not
        }
        # update
        for k, v in args_init.items():
            if hasattr(self, k):
                continue
            setattr(self, k, v)

        # cmd
        self.bamcoverage = which('bamCoverage')

        # bam
        if not file_exists(self.bam):
            self.flag = False
        Bam(self.bam).index()

        # name
        self.prefix = os.path.splitext(os.path.basename(self.bam))[0]
        self.bw_fwd = os.path.join(self.outdir, self.prefix + '.fwd.bigWig')
        self.bw_rev = os.path.join(self.outdir, self.prefix + '.rev.bigWig')
        self.bw = os.path.join(self.outdir, self.prefix + '.bigWig')
        self.bw_fwd_log = os.path.join(self.outdir, self.prefix + '.fwd.deeptools.log')
        self.bw_rev_log = os.path.join(self.outdir, self.prefix + '.rev.deeptools.log')
        self.bw_log = os.path.join(self.outdir, self.prefix + '.deeptools.log')

        # outdir
        check_path(self.outdir)

        # binsize
        if not isinstance(self.binsize, int):
            log.warning('binsize: int expected. {} got, auto-set: 50'.format(self.binsize))
            self.binsize = 50

        # strand
        if not self.strandness in [0, 1, 2, 12]:
            log.warning('strandness: [0|1|2|12] expected, {} got, auto-set: 0'.format(self.strandness))
            self.strandness = 0

        # overwrite
        if not isinstance(self.overwrite, bool):
            log.warning('overwrite: bool expected, {} got, auto-set: False'.format(self.overwrite))
            self.overwrite = False

        # genome/genome_size/reference
        if self.genome in effsize:
            self.genome_size = effsize.get(self.genome, None)
        elif isinstance(self.genome_size, int):
            pass
        elif not self.reference is None:
            self.genome_size = fasize(self.reference, sum_all=True)
        else:
            log.error('genome, genome_size, reference; one of arg is requred')
            self.flag = False # do not run


    def run_cmd(self, cmd, bw, bw_log):
        if os.path.exists(bw) and not self.overwrite:
            log.info('file exists, bigWig skipped ..., {}'.format(bw))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(bw_log, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')
                # p1 = subprocess.run(shlex.split(cmd), stdout=w, stderr=w)
        if not os.path.exists(bw):
            log.error('Bam2bw() failed, bw not found: {}'.format(self.bw))
            # raise ValueError('output file is missing, check log file: %s' % bw_log)


    def run_fwd(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_fwd),
            '--binSize {}'.format(self.binsize), 
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--filterRNAstrand forward',
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing {}'.format(self.normalizeUsing),
            '2> {}'.format(self.bw_log)
            ])
        cmd_txt = os.path.join(self.outdir, 'cmd_fwd.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        self.run_cmd(cmd, self.bw_fwd, self.bw_fwd_log)


    def run_rev(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_rev),
            '--binSize {}'.format(self.binsize), 
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--filterRNAstrand reverse',
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing {}'.format(self.normalizeUsing)])
        cmd_txt = os.path.join(self.outdir, 'cmd_rev.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + 'n')

        self.run_cmd(cmd, self.bw_rev, self.bw_rev_log)
        

    def run_non(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw),
            '--binSize {}'.format(self.binsize),
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing  {}'.format(self.normalizeUsing)])

        cmd_txt = os.path.join(self.outdir, 'cmd.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        self.run_cmd(cmd, self.bw, self.bw_log)
        

    def run(self):
        """
        determined by: strandness
        exit by: flag
        """
        if self.flag is False:
            log.error('bam2bw skipped, ...')
        elif self.strandness == 1:
            self.run_fwd()
        elif self.strandness == 2:
            self.run_rev()
        elif self.strandness == 12:
            self.run_fwd()
            self.run_rev()
        elif self.strandness == 0:
            self.run_non()
        else:
            log.error('unknown strandness')


class Bam2bw2(object):
    """
    Convert bam to bigWig
    version-2: bedtools genomecov -bdg , bedGraphToBigWig, ...

    !!! strandness
    default: dUTP-based library (read2 is sense strand, read1 is anti-sense strand)
    general RNA library: (NSR, small-RNA-library), read1 is sense, read2 is antisense
    """
    def __init__(self, **kwargs):
        """
        bam, file to bam file
        outdir, str, path to save all files
        genome, dm3,dm6,...
        binsize, int, default: 50
        strandness, 0, [0|1|2|12], 1=fwd, 2=rev, 12=rev,rev, 0=non-strandness
        reference, str, fasta file
        overwrite, bool, False

        Required:
        bam, outdir, genome, binsize, strandness, overwrite, reference
        """
        # update self (obj)
        for k, v in kwargs.items():
            setattr(self, k, v)

        # update, defaults
        self.init_args()
        args_logger(self.__dict__, self.config_txt)


    def init_args(self):
        args_init = {
            'bam': None,
            'outdir': str(pathlib.Path.cwd()),
            'binsize': 50,
            'strandness': 0,
            'overwrite': False,
            'reference': None,
            'genome': None,
            'genome_size': None,
            'norm_size': 1000000, # 1 million
            'fragment': True,
            'threads': 4,
            'config_txt': os.path.join(self.outdir, 'arguments.txt'),
            'flag': True # whether run/not
        }
        # update
        for k, v in args_init.items():
            if hasattr(self, k):
                continue
            setattr(self, k, v)

        # cmd
        self.bamcoverage = which('bamCoverage')

        # bam
        if not file_exists(self.bam):
            self.flag = False
        Bam(self.bam).index()

        # name
        self.prefix = os.path.splitext(os.path.basename(self.bam))[0]
        self.bw_fwd = os.path.join(self.outdir, self.prefix + '.fwd.bigWig')
        self.bw_rev = os.path.join(self.outdir, self.prefix + '.rev.bigWig')
        self.bw = os.path.join(self.outdir, self.prefix + '.bigWig')
        self.bw_fwd_log = os.path.join(self.outdir, self.prefix + '.fwd.deeptools.log')
        self.bw_rev_log = os.path.join(self.outdir, self.prefix + '.rev.deeptools.log')
        self.bw_log = os.path.join(self.outdir, self.prefix + '.deeptools.log')

        # outdir
        check_path(self.outdir)

        # binsize
        if not isinstance(self.binsize, int):
            log.warning('binsize: int expected. {} got, auto-set: 50'.format(self.binsize))
            self.binsize = 50

        # norm size
        # default: 1 million
        if isinstance(self.norm_size, int):
            if self.norm_size < 1:
                self.norm_size = 1000000 # 1 million
        else:
            self.norm_size = 1000000 # 1 million

        # strand
        if not self.strandness in [0, 1, 2, 12]:
            log.warning('strandness: [0|1|2|12] expected, {} got, auto-set: 0'.format(self.strandness))
            self.strandness = 0

        # overwrite
        if not isinstance(self.overwrite, bool):
            log.warning('overwrite: bool expected, {} got, auto-set: False'.format(self.overwrite))
            self.overwrite = False

        # genome/genome_size/reference
        if self.genome in effsize:
            self.genome_size = effsize.get(self.genome, None)
        elif isinstance(self.genome_size, int):
            pass
        elif not self.reference is None:
            self.genome_size = fasize(self.reference, sum_all=True)
        else:
            log.error('genome, genome_size, reference; one of arg is requred')
            self.flag = False # do not run


    def get_reads(self, fragment=False):
        """
        fragment, pair of PE reads
        normalize to 1 million mapped reads
        RPM, CPM
        """
        b = Bam(self.bam)
        b_is_paired = b.isPaired()
        b_reads = b.getNumReads()
        if fragment and b_is_paired:
            b_frag = b_reads / 2
        else:
            b_frag = b_reads

        return b_frag


    def get_scale(self):
        """
        Norm to 1 million reads
        """
        r = self.get_reads(fragment=self.fragment)

        # norm by reads/fragment
        return r/self.norm_size


    def run_cmd(self, cmd, bw, bw_log):
        if os.path.exists(bw) and not self.overwrite:
            log.info('file exists, bigWig skipped ..., {}'.format(bw))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(bw_log, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')
                # p1 = subprocess.run(shlex.split(cmd), stdout=w, stderr=w)
        if not os.path.exists(bw):
            log.error('Bam2bw() failed, bw not found: {}'.format(self.bw))
            # raise ValueError('output file is missing, check log file: %s' % bw_log)


    def run_fwd(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_fwd),
            '--binSize {}'.format(self.binsize), 
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--filterRNAstrand forward',
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing {}'.format(self.normalizeUsing)])
        cmd_txt = os.path.join(self.outdir, 'cmd_fwd.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        self.run_cmd(cmd, self.bw_fwd, self.bw_fwd_log)


    def run_rev(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_rev),
            '--binSize {}'.format(self.binsize), 
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--filterRNAstrand reverse',
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing {}'.format(self.normalizeUsing)])
        cmd_txt = os.path.join(self.outdir, 'cmd_rev.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + 'n')

        self.run_cmd(cmd, self.bw_rev, self.bw_rev_log)
        

    def run_non(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw),
            '--binSize {}'.format(self.binsize),
            '--effectiveGenomeSize {}'.format(self.genome_size),
            '--numberOfProcessors {}'.format(self.threads), 
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing  {}'.format(self.normalizeUsing)])

        cmd_txt = os.path.join(self.outdir, 'cmd.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        self.run_cmd(cmd, self.bw, self.bw_log)
        

    def run(self):
        """
        determined by: strandness
        exit by: flag
        """
        if self.flag is False:
            log.error('bam2bw skipped, ...')
        elif self.strandness == 1:
            self.run_fwd()
        elif self.strandness == 2:
            self.run_rev()
        elif self.strandness == 12:
            self.run_fwd()
            self.run_rev()
        elif self.strandness == 0:
            self.run_non()
        else:
            log.error('unknown strandness')


class Bam2cor(object):
    """
    Compute correlation between replicates: 2 or more

    input: bam
    output: count_matrix
            cor_matrix
            cor_plot
            ...

    window = 500bp
    
    eg:
    multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
        --outRawCounts *counts.tab -b bam
    """
    def __init__(self, **kwargs):
        # update self (obj)
        for k, v in kwargs.items():
            setattr(self, k, v)

        # update, defaults
        self.init_args()
        args_logger(self.__dict__, self.config_txt)


    def init_args(self):
        args_init = {
            'multibamsummary': shutil.which('multiBamSummary'),
            'plotcorrelation': shutil.which('plotCorrelation'),
            'plotpca': shutil.which('plotPCA'),
            'cor_method': 'pearson',
            'make_plot': True,
            'bam': None,
            'bam_list': '',
            'outdir': str(pathlib.Path.cwd()),
            'prefix': None,
            'binsize': 500,
            'threads': 1,
            'overwrite': False,
            'config_txt': os.path.join(self.outdir, 'arguments.txt'),
            'flag': True # whether run/not
        }
        # update
        for k, v in args_init.items():
            if hasattr(self, k):
                continue
            setattr(self, k, v)

        # cmd
        if getattr(self, 'multibamsummary', None) is None:
            self.flag = False # do not run
        if getattr(self, 'plotcorrelation', None) is None:
            self.flag = False # do not run
        if getattr(self, 'plotpca', None) is None:
            self.flag = False # do not run

        # outname
        if self.prefix is None:
            self.prefix = 'multibam'

        # outdir
        check_path(self.outdir)

        # bam
        if self.bam is None:
            self.flag = False
        elif isinstance(self.bam, list):
            if len(self.bam) > 1:
                pass
            else:
                self.flag = False
                log.warning('bam: >1 bam files required')
        else:
            self.flag = False
            log.warning('bam: list expected, {} got'.format(type(self.bam).__name__))

        if not all(file_exists(self.bam)):
            self.flag = False
            log.warning('bam file not exists')

        self.bam_list = ' '.join(self.bam) # string

        # index
        [Bam(b).index() for b in self.bam] # 

        # files
        self.config_txt = os.path.join(self.outdir, 'arguments.txt')
        self.bam_npz = os.path.join(self.outdir, self.prefix + '.npz')
        self.log = os.path.join(self.outdir, self.prefix + '.deeptools.log')
        # heatmap
        self.plot_cor_heatmap_png = os.path.join(self.outdir, self.prefix + '.cor_heatmap.png')
        self.log_heatmap = os.path.join(self.outdir, self.prefix + '.cor_heatmap.log')
        self.cor_counts = os.path.join(self.outdir, self.prefix + '.cor_counts.tab')        
        self.cor_matrix = os.path.join(self.outdir, self.prefix + '.cor.matrix')
        # PCA plot
        self.plot_cor_pca_png = os.path.join(self.outdir, self.prefix + '.cor_PCA.png')
        self.log_pca = os.path.join(self.outdir, self.prefix + '.cor_PCA.log')        


    def bam_summary(self):
        cmd = ' '.join([
            '{} bins --binSize {}'.format(self.multibamsummary, self.binsize),
            '-p {}'.format(self.threads),
            '--smartLabels -o {}'.format(self.bam_npz),
            '--outRawCounts {}'.format(self.cor_counts),
            '--bamfiles {}'.format(self.bam_list)
        ])
        cmd_txt = os.path.join(self.outdir, self.prefix + '.bam_cor.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        if os.path.exists(self.bam_npz) and not self.overwrite:
            log.warning('file exists: {}'.format(self.bam_npz))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(self.log, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')

        # check
        if not os.path.exists(self.bam_npz):
            log.error('Bam2cor() failed, output file not found: {}'.format(self.bam_npz))


    def cor_heatmap_plot(self):
        """
        plotCorrelation()
        """
        cmd = ' '.join([
            '{} -in {}'.format(self.plotcorrelation, self.bam_npz),
            '--corMethod {}'.format(self.cor_method),
            '--plotTitle "{} Correlation"'.format(self.cor_method),
            '--whatToPlot heatmap --colorMap RdYlBu --plotNumbers',
            '-o {}'.format(self.plot_cor_heatmap_png),
            '--outFileCorMatrix {}'.format(self.cor_matrix)
            ])

        cmd_txt = os.path.join(self.outdir, self.prefix + '.heatmap.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        if os.path.exists(self.plot_cor_heatmap_png) and not self.overwrite:
            log.warning('file exists: {}'.format(self.plot_cor_heatmap_png))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(self.log_heatmap, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')

        # check
        if not os.path.exists(self.plot_cor_heatmap_png):
            log.error('Bam2cor() failed, output file not found: {}'.format(self.plot_cor_heatmap_png))


    def cor_pca_plot(self):
        """
        plotPCA:
        Make PCA plot, for bam_summary
        correlation
        """
        cmd = ' '.join([
            '{} -in {}'.format(self.plotpca, self.bam_npz),
            '-o {}'.format(self.plot_cor_pca_png),
            '-T "PCA for Bam files"'
            ])
        cmd_txt = os.path.join(self.outdir, self.prefix + '.pca.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        if os.path.exists(self.plot_cor_pca_png) and not self.overwrite:
            log.warning('file exists: {}'.format(self.bam_npz))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(self.log_pca, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')

        # check
        if not os.path.exists(self.plot_cor_pca_png):
            log.error('Bam2cor() failed, output file not found: {}'.format(self.plot_cor_pca_png))


    def run(self):
        # run bam summary
        self.bam_summary() 

        # plot
        if self.make_plot is True:
            self.cor_heatmap_plot()
            self.cor_pca_plot()


class Bam2fingerprint(object):
    """
    QC for bam files
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        Required arguments for ATACseq analysis
        """
        args_init = {
            'bam_dir': None,
            'bam_list': None,
            'threads': 8,
            'outdir': None,
            'prefix': None,
            'title': 'BAM_fingerprint',
            'labels': None,
            'overlap': False,
            'overwrite': False
            }
        self = update_obj(self, args_init, force=False)

        ## bam_dir
        if isinstance(self.bam_dir, str):
            self.bam_list = listfile(self.bam_dir, '*.bam')

        ## bam_list
        if not isinstance(self.bam_list, list):
            raise ValueError('bam_list, expect list, got {}'.format(type(self.bam_list).__name__))

        if not all(file_exists(self.bam_list)):
            raise ValueError('bam_list, file not exists')

        ## outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        check_path(self.outdir)

        ## labels
        bam_names = [os.path.splitext(os.path.basename(i))[0] for i in self.bam_list]
        if self.labels is None:
            self.labels = bam_names

        if not isinstance(self.labels, list):
            raise ValueError('labels, expect list, got {}'.format(type(self.labels).__name))

        if abs(len(self.bam_list) - len(self.labels)) > 0:
            log.warning('bam_list, labels are not in equal in length')
            self.labels = bam_names

        ## output files
        prefix = re.sub('[^A-Za-z0-9-.]', '_', self.title)
        prefix = re.sub('_+', '_', prefix) # format string
        if isinstance(self.prefix, str):
            prefix = self.prefix
        else:
            self.prefix = prefix

        self.fp_png = os.path.join(self.outdir, prefix + '.png')
        self.fp_tab = os.path.join(self.outdir, prefix + '.tab')


    def run(self):
        # check bam index
        [Bam(i).index() for i in self.bam_list]

        cmd = ' '.join([
            '{}'.format(shutil.which('plotFingerprint')),
            '--minMappingQuality 30 --skipZeros --numberOfSamples 50000',
            '-p {}'.format(self.threads),
            '--plotFile {}'.format(self.fp_png),
            '--outRawCounts {}'.format(self.fp_tab),
            '-b {}'.format(' '.join(self.bam_list)),
            '--labels {}'.format(' '.join(self.labels))
            ])

        # save cmd
        cmd_txt = os.path.join(self.outdir, self.prefix + '_cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        # run
        if file_exists(self.fp_png) and not self.overwrite:
            log.info('Bam2fingerprint() skipped, file exists: {}'.format(self.fp_png))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('Bam2fingerprint() failed, see {}'.format(self.outdir))


class PeakIDR(object):
    """
    peak
    outdir

    Check IDR
    Irreproducibility Discovery Rate
    """
    def __init__(self, **kwargs):
        # update self (obj)
        for k, v in kwargs.items():
            setattr(self, k, v)

        # update, defaults
        self.init_args()
        args_logger(self.__dict__, self.config_txt)


    def init_args(self):
        args_init = {
            'peak': None,
            'idr_cmd': shutil.which('idr'),
            'outdir': str(pathlib.Path.cwd()),
            'prefix': None,
            'input_type': 'narrowPeak', # broadPeak, bed, gff
            'cor_method': 'pearson', # spearman
            'overwrite': False,
            'flag': True # run all
        }
        # update
        for k, v in args_init.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        # file
        self.config_txt = os.path.join(self.outdir, self.prefix + '.config.txt')

        # outdir
        check_path(self.outdir)

        # cmd
        if getattr(self, 'idr_cmd', None) is None:
            self.flag = False # do not run
            log.warning('idr: command not found')

        # peak
        if isinstance(self.peak, list):
            if len(self.peak) > 1:
                pass
            else:
                self.flag = False
                log.warning('peak: at least 2 peak files required')
        else:
            self.flag = False
            log.warning('peak: list expected, {} got'.format(type(self.peak).__name__))

        if not all(file_exists(self.peak)):
            self.flag = False
            log.warning('peak: file not exists')


    def idr(self, peakA=None, peakB=None):
        """
        Compute IDR for two group of peaks
        """
        if not peakA is None or not peakB is None:
            pass
        else:
            if len(self.peak) >= 2:
                peakA, peakB = self.peak[:2]
            else:
                log.warning('at least 2 peaks required for idr')

        # prefix
        pA = os.path.splitext(os.path.basename(peakA))[0] #
        pB = os.path.splitext(os.path.basename(peakB))[0] #
        prefix = '{}.vs.{}.idr'.format(pA, pB)

        if isinstance(self.prefix, str):
            prefix = self.prefix

        # files
        cmd_txt = os.path.join(self.outdir, prefix + '.cmd.sh')
        idr_txt = os.path.join(self.outdir, prefix + '.txt')
        idr_png = idr_txt + '.png' # os.path.join(self.outdir, prefix + '.png')
        idr_log = os.path.join(self.outdir, prefix + '.log')

        cmd = ' '.join([
            '{} --samples {} {}'.format(self.idr_cmd, peakA, peakB),
            '--input-file-type {}'.format(self.input_type),
            '--plot', #.format(idr_png),
            '--output-file {}'.format(idr_txt),
            '--log-output-file {}'.format(idr_log)
            ])

        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        if os.path.exists(idr_png) and not self.overwrite:
            log.warning('file exists: {}'.format(idr_png))
        else:
            # require, peaks > 20
            peakA_nrow = file_row_counter(peakA)
            peakB_nrow = file_row_counter(peakB)
            if peakA_nrow >= 100 and peakB_nrow >= 100:
                _, stdout, stderr = run_shell_cmd(cmd)
                with open(idr_log, 'wt') as w:
                    w.write(stdout + '\n' + stderr + '\n')
            else:
                log.warning('idr: peak files required at least 100 peaks, {} and {} got'.format(peakA_nrow, peakB_nrow))

        # check
        if not os.path.exists(idr_png):
            log.error('Peak().idr() failed, output file not found: {}'.format(idr_png))


    def run(self):
        """
        Run A, B
        """
        if self.flag is False:
            log.error('Peak().idr() argumnets failed')
        else:
            for peakA, peakB in combinations(self.peak, 2):
                self.idr(peakA, peakB)
    

class PeakFRiP(object):
    """
    Calculate the FRiP 
    see ENCODE: https://www.encodeproject.org/data-standards/terms/#enrichment

    1. bedtools intersect
    2. featureCounts
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'peak': None,
            'bam': None,
            'method': 'bedtools', # featureCounts
            'gsize': None
        }
        self = update_obj(self, args_init, force=False)

        if not file_exists(self.peak):
            raise ValueError('peak, file not exists: {}'.format(self.peak))

        if not file_exists(self.bam):
            raise ValueError('bam, file not exists: {}'.format(self.bam))

        if not self.method in ['bedtools', 'featureCounts']:
            self.method = 'bedtools' # default


    def run_bedtools(self):
        """
        see: https://www.biostars.org/p/337872/#338646

        sort -k1 -k2,2n -o peak.bed peak.bed
        bedtools intersect -c -a peak.bed -b file.bam -sorted -g ref | awk '{i+=$n}END{print i}'
        """
        # total reads
        total = Bam(self.bam).count()

        # peak reads
        rip_txt = self._tmp(delete=False)
        cmd = ' '.join([
            'sort -k1,1 -k2,2n -o {} {}'.format(self.peak, self.peak),
            '&& bedtools intersect -c -a {} -b {} -g {}'.format(
                self.peak, self.bam, self.gsize),
            r"| awk '{i+=$NF}END{print i}'",
            ' > {}'.format(rip_txt)
            ])
        run_shell_cmd(cmd)

        with open(rip_txt, 'rt') as r:
            rip = int(r.read().strip())

        # output
        if total > 0:
            frip = round(rip/total, 4) #  '{:.2f}'.format(rip/total*100)
        else:
            frip = 0

        # fragments (PE)
        if Bam(self.bam).isPaired():
            total = int(total / 2.0)
            rip = int(rip / 2.0)

        # remove temp file
        file_remove(rip_txt, ask=False)
        return (total, rip, frip)


    def run_featureCounts(self):
        """
        see: https://www.biostars.org/p/337872/#337890

        -p , fragments
        """
        fc = FeatureCounts(gtf=self.peak, bam_list=self.bam)
        fc.run()
        df = fc.wrap_log()
        # df = pd.DataFrame(stat).T
        # df.columns = ['total', 'map', 'FRiP']
        return df.values.tolist()[0] # list


    def _tmp(self, delete=True):
        """
        Create a tmp filename
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.txt', delete=delete)
        return tmp.name


    def run(self):
        if self.method == 'bedtools':
            return self.run_bedtools()
        elif self.method == 'featureCounts':
            return self.run_featureCounts()
        else:
            log.error('method, unknown, [bedtools|featureCounts], got {}'.format(
                self.method))


class BedOverlap(object):
    """
    pybedtools API, to calculate the overlaps between bed files

    pybedtools.contrib.venn_maker.venn_maker(beds, names, figure_filename, script_filename, run=True)

    make venn plot
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.init_args()


    def init_args(self):
        """
        required: args
        """
        args_init = {
            'peak': None,
            'outdir': str(pathlib.Path.cwd()),
            'flag': False,
            'prefix': None,
            'overwrite': False
        }
        for k, v in args_init.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        # outdir
        check_path(self.outdir)

        # peak
        if isinstance(self.peak, list):
            if len(self.peak) < 2:
                self.flag = False
                log.warning('peak: >1 files required')
        else:
            self.flag = False
            log.warning('peak: list expected, {} got'.format(type(self.peak).__name__))

        if not all(file_exists(self.peak)):
            self.flag = False
            log.error('peak: file not exists')

        # update peak files
        self.valid_peak() #

        # names
        self.peak_names = [os.path.splitext(i)[0] for i in fq_name(self.peak)]
        if self.prefix is None:
            self.prefix = 'peak_overlap'

        # files
        self.tiff = os.path.join(self.outdir, self.prefix + '.tiff')
        self.png = os.path.join(self.outdir, self.prefix + '.png')
        self.venn_R = os.path.join(self.outdir, self.prefix + '.venn.R')


    def valid_peak(self):
        """
        No more than 4
        Peak: >0
        """
        # file rows > 0
        self.peak = [i for i in self.peak if file_row_counter(i) > 0]

        # files < 4
        if len(self.peak) > 4:
            log.warning('peak: support no more than 4 files, subset to 4')
            self.peak = self.peak[:4]


    def overlap(self):
        """
        Overlap between A and B, ...
        """
        plt = pybedtools.contrib.venn_maker.venn_maker
        if os.path.exists(self.tiff) and self.overwrite is False:
            log.info('overlap file exists, skipped ...')
        elif len(self.peak) > 1:
            log.info('Calculating overlaps between BED files')
            plt(self.peak, self.peak_names, figure_filename=self.tiff, 
                script_filename=self.venn_R, run=True)
        else:
            log.warning('peak: files more than 1 required, {} got'.format(len(self.peak)))

        # tiff -> png
        if os.path.exists(self.png) and self.overwrite is False:
            pass
        else:
            # convert to png
            log.info('Coverting Tiff to png')
            if os.path.exists(self.tiff):
                convert_image(self.tiff, 'PNG')
            else:
                log.warning('tiff, file not found: {}'.format(self.tiff))

 
    def run(self):
        self.overlap()


def convert_image(x, out_fmt='PNG'):
    if not out_fmt in ['PNG', 'JPEG', "TIFF"]:
        log.error('out_fmt: [PNG|JPEG|TIFF], {} got'.format(out_fmt))

    out_ext = out_fmt.lower()
    out_img = os.path.splitext(x)[0] + '.' + out_ext

    # read/write
    if os.path.exists(out_img):
        log.warning('file exists, skipping ...: {}'.format(out_img))
    else:
        img = Image.open(x)
        img.save(out_img, out_fmt)


def bwCompare(bw1, bw2, bw_out, operation='log2', **kwargs):
    """
    Compare two bigWig files: ip over input

    example:
    bigwigCompare -b1 bw1 -b2 bw2 --operation 
    {log2, ratio, subtract, add, mean, reciprocal_ratio, 
    first, second}
    -o out.bw
    """
    binsize = kwargs.get('binsize', 50)
    threads = kwargs.get('threads', 1)
    overwrite = kwargs.get('overwrite', False)

    cmd = ' '.join([
        '{}'.format(shutil.which('bigwigCompare')),
        '--bigwig1 {} --bigwig2 {}'.format(bw1, bw2),
        '--operation {}'.format(operation),
        '--skipZeroOverZero',
        '--skipNAs',
        '--binSize {}'.format(binsize),
        '-p {}'.format(threads),
        '-o {}'.format(bw_out)])

    # savd cmd 
    bw_out_dir = os.path.dirname(bw_out)
    check_path(bw_out_dir)
    cmd_txt = os.path.join(bw_out_dir, 'cmd.txt')
    with open(cmd_txt, 'wt') as w:
        w.write(cmd + '\n')

    # run
    if os.path.exists(bw_out) and not overwrite:
        log.info('bwCompare() skipped, file exists: {}'.format(bw_out))
    else:
        run_shell_cmd(cmd)


