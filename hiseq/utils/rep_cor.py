"""Quality control for HiSeq data

- Correlation: 500 window, pearson correlation  
- peaks, overlap  
- IDR, peaks (ChIP-seq, ATAC-seq)  
...

"""

import os
import sys
import shutil
import pathlib
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


def bam_cor(bam_list, outdir=None, window=500):
    """
    Compute correlation (pearson) between replicates
    window = 500bp
    
    eg:
    multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
        --outRawCounts *counts.tab -b bam
    """
    multiBamSummary = shutil.which('multiBamSummary')

    if outdir is None:
        outdir = str(pathlib.Path.cwd())

    # default files
    check_path(outdir)
    cor_npz = os.path.join(outdir, 'multibam.npz')
    cor_counts = os.path.join(outdir, 'multibam.counts.tab')

    # run
    cmd = ' '.join([
        multiBamSummary,
        'bins --binSize {}'.format(window),
        '-p {}'.format(8),
        '--smartLabels -o {}'.format(cor_npz),
        '--outRawCounts {}'.format(cor_counts),
        '--bamfiles {}'.format(' '.join(bam_list))])

    # save cmd
    cmd_txt = os.path.join(outdir, "cmd.sh")
    with open(cmd_txt, 'wt') as w:
        w.write(cmd + '\n')

    if os.path.exists(cor_counts):
        log.info('bam_cor() skipped, file.exsits: {}'.format(
            cor_counts))
    else:
        try:
            # make a index
            [Bam(i).index() for i in bam_list]
            run_shell_cmd(cmd)
        except:
            log.warning('bam_cor() failed.')

    return cor_npz


def cor_plot(npz):
    npz_dir = os.path.dirname(npz)
    cor_pdf = os.path.join(npz_dir, 'heatmap.pdf')
    cor_tab = os.path.join(npz_dir, 'heatmap.cor.tab')
    pca_pdf = os.path.join(npz_dir, 'PCA.pdf')

    # correlation
    cmd1 = ' '.join([
        'plotCorrelation -in {}'.format(npz),
        '--corMethod pearson --skipZeros',
        '--whatToPlot heatmap --plotNumbers',
        '-o {} --outFileCorMatrix {}'.format(cor_pdf, cor_tab)
        ])

    # save cmd
    cmd1_txt = os.path.join(npz_dir, 'cmd_cor.sh')
    with open(cmd1_txt, 'wt') as w:
        w.write(cmd1 + '\n')

    # run
    if os.path.exists(cor_pdf):
        log.info('file exists, skipped: {}'.format(cor_pdf))
    else:
        try:
            run_shell_cmd(cmd1)
        except:
            log.warning('plotCorrealtion failed')

    # PCA
    cmd2 = ' '.join([
        'plotPCA -in {}'.format(npz),
        '-o {} -T "BAM PCA" '.format(pca_pdf)
        ])

    # save cmd
    cmd2_txt = os.path.join(npz_dir, 'cmd_pca.sh')
    with open(cmd2_txt, 'wt') as w:
        w.write(cmd2 + '\n')

    # run
    if os.path.exists(pca_pdf):
        log.info('file exists, skipped: {}'.format(pca_pdf))
    else:
        try:
            run_shell_cmd(cmd2)
        except:
            log.warning('plotCorrealtion failed')


def main():
    if len(sys.argv) < 3:
        sys.exit('Usage: rep_cor.py <outdir> <bam1> <bam2> ...')

    outdir = sys.argv[1]
    bam_list = sys.argv[2:]

    npz = bam_cor(bam_list, outdir)
    cor_plot(npz)


if __name__ == '__main__':
    main()

