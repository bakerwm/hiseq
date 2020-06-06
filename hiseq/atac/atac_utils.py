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
import pysam
import pybedtools
from collections import OrderedDict
from hiseq.utils.helper import *
from hiseq.utils.seq import Fastx


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
    if absolute_path:
        # support: ~, $HOME,
        srcname = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
    else:
        # only for directories within the same folder
        srcname = os.path.join('..', os.path.basename(src))

    if not os.path.exists(dest):
        os.symlink(srcname, dest)


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
            stdout, stderr = run_shell_cmd(cmd)
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
            '--filterRNAstrand forward',
            '--binSize {} --effectiveGenomeSize {}'.format(self.binsize, self.genome_size)])
        cmd_txt = os.path.join(self.outdir, 'cmd_fwd.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        self.run_cmd(cmd, self.bw_fwd, self.bw_fwd_log)


    def run_rev(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw_rev),
            '--filterRNAstrand reverse',
            '--binSize {} --effectiveGenomeSize {}'.format(self.binsize, self.genome_size)])
        cmd_txt = os.path.join(self.outdir, 'cmd_rev.txt')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + 'n')

        self.run_cmd(cmd, self.bw_rev, self.bw_rev_log)
        

    def run_non(self):
        # bedtools cmd
        cmd = ' '.join([
            '{} -b {} -o {}'.format(self.bamcoverage, self.bam, self.bw),
            '--binSize {} --effectiveGenomeSize {}'.format(self.binsize, self.genome_size)])
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


def cal_FRiP(inbed, inbam, genome="dm6"):
    """Calculate FRiP for ChIP-seq/ATAC-seq peaks
    Fraction of reads in called peak regions
    input: bam, peak (bed)

    using: bedtools intersect 
    count reads in peaks
    """
    # total map
    # index
    Bam(inbam).index()
    bam = pysam.AlignmentFile(inbam)
    total = bam.mapped

    gsize = Genome(genome).get_fasize()
    tmp = os.path.basename(inbed) + '.count.tmp'

    try:
        # reads in peak
        cmd = ' '.join([
            'sort -k1,1 -k2,2n {}'.format(inbed),
            '| bedtools intersect -c -a - -b {}'.foramt(inbam),
            "| awk '{s+=$NF}END{print s}' > {}".format(tmp)])
        run_shell_cmd(cmd)

         # "sort -k1,1 -k2,2n {} | \
         #    bedtools intersect -c -a - -b {} | \
         #    {} > {}".format(
         #    inbed, inbam, "awk '{s+=$NF}END{print s}'", tmp)        
        # run_shell_cmd(p)

        with open(tmp) as r:
            n = eval(r.readline().strip())
        frip = '{:.2f}%'.format(n/total*100)

        # clean
        os.remove(tmp)
    except:
        log.warning('cal_FRiP() failed, skip {}'.format(inbed))
        if os.path.exists(tmp):
            os.remove(tmp)
        frip = n = total = 0 # init

    return (frip, n, total)


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


