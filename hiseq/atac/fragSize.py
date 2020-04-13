#!/usr/bin/env python3

"""
Calculate the Fragment size (insert size) of BAM file
1. read.tlen
2. ...
"""

import os, sys
from hiseq.atac.atac_utils import frag_length
from hiseq.utils.helper import *


## method-1
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


## method-2
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


def main():
    if len(sys.argv) < 3:
        sys.exit('python fragSize.py <in.bam> <out.txt>')
    bam_in, stat_txt = sys.argv[1:3]

    frag_length(bam_in, stat_txt)


if __name__ == '__main__':
    main()
