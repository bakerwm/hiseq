#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This is the main script for regular usage,
call sub-commands
...
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2021-05-20'
__version__ = '1.0.1'

# import hiseq
import os
import sys
import argparse
from multiprocessing import Pool

# main args module modules
from hiseq.atac.atac import Atac
from hiseq.atac.atac import get_args as add_atac_args
from hiseq.cnr.cnr import Cnr
from hiseq.cnr.cnr import get_args as add_cnr_args


# to-be-deprecated: replaced by specific get_args() in each command
from hiseq.utils.argsParser import add_sheet_args, add_demx_args, \
    add_demx2_args, add_qc_args, add_p7_args, add_trim_args, add_align_args, \
    add_quant_args, add_peak_args, add_motif_args, add_rnaseq_args, \
    add_rnaseq_args2, add_chipseq_args, \
    add_trackhub_args, add_deseq_pair_args, add_go_args, add_fragsize_args, \
    add_bam2bw_args, add_bam2cor_args, add_peak2idr_args, add_bed2overlap_args, \
    add_sample_args

from hiseq.demx.sample_sheet import SampleSheet
from hiseq.demx.demx import Demx, Demx2
from hiseq.qc.fastqc import Fastqc
from hiseq.qc.parse_i7 import HiSeqP7
from hiseq.trim.trimmer import Trim
from hiseq.align.align import Align
from hiseq.rnaseq.rnaseq import RNAseq
from hiseq.rnaseq.deseq_pair import DeseqPair
from hiseq.chipseq.chipseq import ChIPseq
from hiseq.go.go import Go
from hiseq.utils import download as dl # download.main()
from hiseq.utils.fastx import Fastx
from hiseq.bam2bw.bam2bw import Bam2bw
from hiseq.utils.bam import Bam2cor
from hiseq.bam2bw.bam2bw import Bam2bw
from hiseq.utils.bed import PeakIDR, BedOverlap
from hiseq.fragsize.fragsize import BamPEFragSize2
from hiseq.get_trackhub.get_trackhub import TrackHub
from hiseq.peak.call_peak import Macs2
from hiseq.sample.sample import FxSample
from hiseq.utils.utils import log, get_date


class Hiseq(object):
    """The 1st-level of command, choose which sub-command to use
    qc, trim, align, quant, peak, motif, report, ... (to be continued)

    """
    def __init__(self):
        usage = """ hiseq <command> [<args>]

    subcommands:

        atac         ATACseq pipeline
        rnaseq       RNAseq pipeline
        rnaseq2      RNAseq pipeline, simplify version
        chipseq      ChIPseq pipeline
        cnr          CUN&RUN pipeline

        sheet        Preparing sample_sheet.csv for demx/demx2
        demx         Demultiplexing reads (P7, barcode) from single Lane
        demx2        Demultiplexing multi barcode files
        qc           quality control, fastqc
        p7           Check the P7 of HiSeq library
        
        trim         trim adapters, low-quality bases, ...
        align        Align fastq/a files to reference genome
        quant        Count genes/features
        peak         Call peaks using MACS2
        motif        Check motifs from a BED/fasta file
        report       Create a report to the above commands
        go           Run GO analysis on geneset
        deseq_pair   Run RNAseq compare

        run_trackhub Generate the trackhub urls
        fragsize     Fragment size of PE alignments
        bam2cor      Correlation between bam files
        bam2bw       Convert bam to bigWig
        peak2idr     Calculate IDR for multiple Peaks
        bed2overlap  Calculate the overlap between bed intervals
        sample       Sample fastq file
        download     Download files
    """
        parser = argparse.ArgumentParser(
            prog = 'hiseq',
            description = 'A collection of tools for HiSeq data',
            epilog = '',
            usage = usage,
        )
        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            log.error('unknown command: {}'.format(args.command))
            parser.print_help()
            sys.exit(1)
        getattr(self, args.command)()
        
    
    def init_args(self, p):
        """
        Parameters
        ----------
        p:
            argparse.ArgumentPaser
        """
        if isinstance(p, argparse.ArgumentParser):
            out = vars(p.parse_args(sys.argv[2:]))
        else:
            out = None
        return out


    def atac(self):
        """
        ATACseq pipeline
        """
        args = self.init_args(add_atac_args())
        Atac(**args).run()

    
################################################################################
## to-be-updated
## 2021-05-20
    
    
    def sheet(self):
        """
        Prepare the sample sheet for demx
        1. sample_name,i7,i5,barcode (Demx, whole-lane)
        2. i7_name,i7,reads (Demx2, part)
        """
        args = self.init_args(add_sheet_args())
        SampleSheet(**args).run()


    def demx(self):
        """
        Demx-1: one large to small files
        """
        args = self.init_args(add_demx_args())
        Demx(**args).run()


    def demx2(self):
        """
        Dexm-2: rename files, demx p7, ...
        """
        args = self.init_args(add_demx2_args())
        Demx2(**args).run()


    def qc(self):
        """
        Fastq quality control
        """
        args = self.init_args(add_qc_args())
        Fastqc(**args).run()

        
    def p7(self):
        """
        Check library structure: P7, barcode
        """
        args = self.init_args(add_p7_args())
        HiSeqP7(**args).run()
        

    def trim(self):
        """
        Trim adapters
        """
        args = self.init_args(add_trim_args())
        TrimRn(**args).run()


    def align(self):
        """
        Align reads
        """
        args = self.init_args(add_align_args())
        Align(**args).run()


    def quant(self):
        """
        quantify hiseq reads
        """
        args = self.init_args(add_quant_args())
        print('hiseq quant')


    def peak(self):
        """
        Call peaks
        """
        args = self.init_args(add_peak_args())
        CallPeak(**args).run()


    def motif(self):
        """
        quantify hiseq reads
        ...
        """        
        args = self.init_args(add_motif_args())
        print('hiseq motif')


    def go(self):
        """
        Run GO analysis on geneset
        """
        args = self.init_args(add_go_args())
        if args['all'] is None and args['input'] is None:
            parser.parse_args(['-h'])
            sys.exit('arguments failed, either --all or --input required')
        Go(**args).run()


    def deseq_pair(self):
        """
        Run RNAseq cmp
        """
        args = self.init_args(add_deseq_pair_args())
        DeseqPair(**args).run()


    def rnaseq(self):
        """
        RNA-seq pipeline
        """
        args = self.init_args(add_atac_args())
        RNAseq(**args).run()


    def rnaseq2(self):
        """
        RNA-seq pipeline, simplify version
        """
        args = self.init_args(add_rnaseq_args2)
        RNAseqPipe(**args).run()


    def chipseq(self):
        """
        ChIPseq pipeline
        """
        args = self.init_args(add_chipseq_args())
        ChIPseq(**args).run()


    def cnr(self):
        """
        CUN&RUN pipeline
        """
        args = self.init_args(add_cnr_args())
        Cnr(**args).run()


    def fragsize(self):
        """
        Calculate the fragment size of PE alignment
        """
        args = self.init_args(add_fragsize_args())
        BamPEFragSize2(**args).run()


    def bam2bw(self):
        """
        Convert bam to bw files
        """
        args = self.init_args(add_bam2bw_args())
        Bam2bw(**args_local).run()


    def bam2cor(self):
        """
        Calculate bam correlation
        using deeptools
        """
        args = self.init_args(add_bam2cor_args())
        # args['make_plot'] = not args.get('no_plot', False)
        Bam2cor(**args).run()


    def peak2idr(self):
        """
        Calculate IDR for peak files
        using: idr
        """
        args = self.init_args(add_peak2idr_args())
        PeakIDR(**args).run()


    def bed2overlap(self):
        """
        Calculate IDR for peak files
        using: idr
        """
        args = self.init_args(add_bed2overlap_args())
        BedOverlap(**args).run()


    def sample(self):
        """
        Sub-sample fastq files
        head
        """
        args = self.init_args(add_sample_args())
        FxSample(**args).run()


    def download(self):
        """
        Download files
        """
        dl.main(sys.argv[1:])


    def run_trackhub(self):
        """
        Make trackhub
        """
        args = self.init_args(add_trackhub_args())
        TrackHub(**args).run()


def main():
    Hiseq()


if __name__ == '__main__':
    main()

#