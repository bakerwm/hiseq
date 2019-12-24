"""Quality control for HiSeq data

- Correlation: 500 window, pearson correlation  
- peaks, overlap  
- IDR, peaks (ChIP-seq, ATAC-seq)  
...

"""

import os
import sys
import pysam
import pybedtools
from collections import OrderedDict
from hiseq.utils.helper import *


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


def cal_FRiP(inbed, inbam, genome="dm6"):
    """Calculate FRiP for ChIP-seq/ATAC-seq peaks
    Fraction of reads in called peak regions
    input: bam, peak (bed)
    """
    # total map
    # index
    Bam(inbam).index()
    bam = pysam.AlignmentFile(inbam)
    total = bam.mapped

    gsize = Genome(genome).get_fasize()
    tmp = os.path.basename(inbed) + '.count.tmp'
    # reads in peak
    p = "sort -k1,1 -k2,2n {} | \
        bedtools intersect -c -a - -b {} | \
        {} > {}".format(
        inbed, inbam, "awk '{s+=$NF}END{print s}'", tmp)
    
    run_shell_cmd(p)

    with open(tmp) as r:
        n = eval(r.readline().strip())

    frip = '{:.2f}%'.format(n/total*100)

    # clean
    os.remove(tmp)

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
    names = [os.path.splitext(i)[0] for i in peak_list]

    # name
    pname = os.path.splitext(peak_list[0])[0]
    pname = pname.rstrip('_peaks')
    pname = pname.rstrip('_rep1')

    # output
    prefix = os.path.join(outdir, pname)
    vennR = prefix + '.venn.R'
    tiffout = prefix + '.venn.png'

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
    idr = shutil.which('idr')
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


