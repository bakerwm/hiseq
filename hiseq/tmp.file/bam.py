#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function, processing BAM files

BAM:
  - sort 
  - index 
  - count 
  - rmdup 
  - proper_pair 
  - subset 
  - frag_size * 
  - read_size *
  - to_bed
  - to_bg *
  - to_bw *
  - to_fq *
  - to_xx *
"""

import os
import numpy as np
import pysam
import pybedtools
import logging
from hiseq.utils.helper import *

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
        x = pysam.view('-c', self.bam)
        return int(x.strip())


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


    def subset(self, size=20000, subdir=None):
        """
        Extract N reads from bam list
        """
        if subdir is None:
            subdir = self._tmp(delete=False)
        check_path(subdir)

        # src, dest
        dest = os.path.join(subdir, os.path.basename(self.bam))
        
        # run
        self.index()        
        srcfile = pysam.AlignmentFile(self.bam, 'rb')
        destfile = pysam.AlignmentFile(dest, 'wb', template=srcfile)
        # counter
        i = 0
        for read in srcfile.fetch():
            i +=1
            if i > size:
                break
            destfile.write(read)

        return dest


    def _tmp(self, delete=True):
        """
        Create a tmp filename
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=delete)
        return tmp.name


    ##########################################
    ## code from cgat: BEGIN
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
        from pysam
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
            raw_median = np.median(inserts)
            raw_median_dev = np.median(np.absolute(inserts - raw_median))

            # set thresholds
            threshold_min = max(0, raw_median - n * raw_median_dev)
            threshold_max = raw_median + n * raw_median_dev

            # define core distribution
            return inserts[np.logical_and(inserts >= threshold_min,
                                             inserts <= threshold_max)]

        if method == "picard":

            # only get first read in pair to avoid double counting
            inserts = np.array(
                [read.template_length for read in samfile.head(n=topn)
                 if read.is_proper_pair
                 and not read.is_unmapped
                 and not read.mate_is_unmapped
                 and not read.is_read1
                 and not read.is_duplicate
                 and read.template_length > 0])
            core = get_core_distribution(inserts, n)

            return np.mean(core), np.std(core), len(inserts)

        elif method == "convergence":

            means, stds, counts = [], [], []
            last_mean = 0
            iteration = 0
            while iteration < max_chunks:

                inserts = np.array(
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
                means.append(np.mean(core))
                stds.append(np.std(core))
                counts.append(len(inserts))
                mean_core = get_core_distribution(np.array(means), 2)
                mm = np.mean(mean_core)
                if abs(mm - last_mean) < similarity_threshold:
                    break
                last_mean = mm

            return np.median(means), np.median(stds), sum(counts)
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
        """
        Return number of alignments in bamfile.
        """        
        if not os.path.exists(self.bam + '.bai'):
            pysam.index(self.bam)
        samfile = pysam.AlignmentFile(self.bam)
        return samfile.mapped
    ##########################################
    ## code from cgat: END
    ##########################################