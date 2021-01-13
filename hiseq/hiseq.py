#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This is the main script for regular usage,
call sub-commands

...
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08@hotmail.com'
__date__ = '2019-10-01'
__version__ = '0.0.1'

# import hiseq
import os
import sys
import argparse
from multiprocessing import Pool
from .utils.argsParser import *
from .demx.demx import Demx
from .qc.fastqc import Fastqc
from .trim.trimmer import TrimRn
from .align.alignment import Alignment
from .atac.atac import Atac
from .rnaseq.rnaseq import RNAseq
from .rnaseq.rnaseq_pipe import RNAseqPipe
from .rnaseq.deseq_pair import DeseqPair
from .chipseq.chipseq import ChIPseq
from .cnr.cnr import CnR
from .go.go import Go
from .utils import download as dl # download.main()
from .utils.seq import Fastx
from .atac.atac_utils import Bam2bw, Bam2cor, PeakIDR, BedOverlap
from .utils.helper import *
from .fragsize.fragsize import BamPEFragSize2
from .get_trackhub.get_trackhub import TrackHub 
# from .peak.call_peak import Macs2
from hiseq.peak.call_peak import Macs2
from hiseq.sample.sample import FxSample


class Hiseq(object):
    """The 1st-level of command, choose which sub-command to use
    qc, trim, align, quant, peak, motif, report, ... (to be continued)

    """
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog = 'hiseq',
            description = 'A collection of tools for HiSeq data',
            epilog = '',
            usage = """ hiseq <command> [<args>]

    The most commonly used sub-commands are:

        atac         ATACseq pipeline
        rnaseq       RNAseq pipeline
        rnaseq2      RNAseq pipeline, simplify version
        chipseq      ChIPseq pipeline
        cnr          CUN&RUN pipeline

        demx         Demultiplexing reads (P7, barcode)
        qc           quality control, fastqc
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
        )
        parser.add_argument('command', help='Subcommand to run')

        # parse_args defaults to [1:] for args
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            sys.exit(1)

        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()


    def demx(self):
        """
        Demultiplexing reads: P7, barcode
        """
        parser = add_demx_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
        Demx(**args).run()


    def qc(self):
        """
        Quality control
        fastqc    : Create fastqc report the raw and clean files
        """
        parser = add_qc_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args) # convert to dict

        Fastqc(**args).run()


    def trim(self):
        """
        Trim adapters
        """
        parser = add_trim_args()
        args = vars(parser.parse_args(sys.argv[2:]))
        TrimRn(**args).run()


    # def trim(self):
    #     """Quality control,
    #     cutadapt  : Trimming adapters from reads
    #     trim_ends : Python scripts to trim ends of reads, usually, the barcode sequences
    #     fastqc    : Create fastqc report the raw and clean files
    #     """
    #     parser = add_trim_args()
    #     args = parser.parse_args(sys.argv[2:])
    #     args = vars(args) # convert to dict

    #     # custom args
    #     fq1_list = args.get('fq1', None) # list
    #     outdir = args.get('outdir', None)

    #     ## iterate all fq1
    #     kwargs = args.copy()
    #     if kwargs['fq2'] is None:
    #         # SE mode
    #         for fq1 in fq1_list:
    #             kwargs['fq1'] = fq1
    #             kwargs['fq2'] = None
    #             kwargs['outdir'] = args.get('outdir', None)
    #             Trimmer(**kwargs).run()
    #     else:
    #         # PE mode
    #         if not len(fq1_list) == len(args['fq2']):
    #             log.error('-i, --fq2 not in the same length')

    #         for fq1, fq2 in zip(fq1_list, args['fq2']):
    #             kwargs['fq1'] = fq1
    #             kwargs['fq2'] = fq2
    #             kwargs['outdir'] = args.get('outdir', None)
    #             Trimmer(**kwargs).run()


    def align(self):
        """Alignment
        Align fastq/a files to reference
        fq: SE or PE
        aligner: bowtie, bowtie2, STAR, bwa, hisat2, ...
        output: directory
        args: unique, multi, x-size, ...
        input
        output
        arguments.txt (pickle)
        """
        parser = add_align_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args) # convert to dict
        # print('Running hiseq align, aligner={}, fq1={}, out={}'.format(
        #     args.aligner, args.fq1, args.out))
        Alignment(**args).run()


    def quant(self):
        """
        quantify hiseq reads
        ...
        """
        parser = add_quant_args()
        args = parser.parse_args(sys.argv[2:])
        print('hiseq quant')


    def peak(self):
        """
        quantify hiseq reads
        ...
        """
        parser = add_peak_args()
        args = vars(parser.parse_args(sys.argv[2:]))
        genome = args.get('genome', 'dm6')
        outdir = args.get('outdir', str(pathlib.Path.cwd()))
        bam_list = args.get('bam', None)
        control_list = args.get('control', None)
        name_list = args.get('name', None)

        if bam_list is None:
            log.warning('-i not detect bam files')
        else:
            for i, bam in enumerate(bam_list):
                # control
                if not control_list is None:
                    control = control_list[i]
                else:
                    control = None

                # prefix
                if not name_list is None:
                    name = name_list(i)
                else:
                    name = os.path.basename(os.path.splitext(bam)[0])

                # output
                subdir = os.path.join(outdir, name)
                peak = os.path.join(subdir, name + '_peaks.narrowPeak')

                if os.path.exists(peak):
                    log.info('file exists: {}'.format(peak))
                else:
                    Macs2(
                        ip=bam,
                        genome=args.get('genome', 'dm6'),
                        output=subdir,
                        prefix=name,
                        control=control,
                        atac=args.get('is_atac', False), 
                        gsize=args.get('genome_size', 0),
                        overwrite=args.get('overwrite', False)).callpeak()
        

    def motif(self):
        """
        quantify hiseq reads
        ...
        """
        parser = add_motif_args()
        args = parser.parse_args(sys.argv[2:])
        print('hiseq motif')


    def go(self):
        """
        Run GO analysis on geneset
        """
        parser = add_go_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args) # convert to dict
        # print('Running hiseq align, aligner={}, fq1={}, out={}'.format(
        #     args.aligner, args.fq1, args.out))
        if args['all'] is None and args['input'] is None:
            parser.parse_args(['-h'])
            sys.exit('arguments failed, either --all or --input required')
        Go(**args).run()


    def deseq_pair(self):
        """
        Run RNAseq cmp
        """
        parser = add_deseq_pair_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args) # convert to dict
        DeseqPair(**args).run()


    def atac(self):
        """
        ATACseq pipeline
        """
        parser = add_atac_args()
        args = parser.parse_args(sys.argv[2:])
        # help
        if len(sys.argv) < 3:
            parser.parse_args(['-h'])
        # main
        args = vars(args) # convert to dict
        # check config or --fq1,--fq2,--genome,--outdir
        config = args.get('config', None)
        design = args.get('design', None)
        fq1 = args.get('fq1', None)
        fq2 = args.get('fq2', None)
        genome = args.get('genome', None)
        outdir = args.get('outdir', None)
        chk1 = config is None
        chk2 = design is None
        chk3 = [i is None for i in [fq1, fq2, genome, outdir]]

        # if config is None and not all(chk2):
        if all([chk1, chk2, chk3]):
            sys.exit('required: --config, or --design, or --fq1, --fq2, --genome, --outdir')

        # a = AtacBatch(**args).run()
        Atac(**args).run()


    def rnaseq(self):
        """
        RNA-seq pipeline
        """
        parser = add_rnaseq_args()
        args = parser.parse_args(sys.argv[2:])
        # help
        if len(sys.argv) < 3:
            parser.parse_args(['-h'])
        args = vars(args) # convert to dict
        # check design or --fq1,--fq2,--genome,--outdir or smp_path/dirs_ctl/dirs_exp
        # design = args.get('design', None)
        # fq1 = args.get('fq1', None)
        # fq2 = args.get('fq2', None)
        # genome = args.get('genome', None)
        # outdir = args.get('outdir', None)
        # smp_path = args.get('smp_path', None)
        # dirs_ctl = args.get('dirs_ctl', None)
        # dirs_exp = args.get('dirs_exp', None)
        # chk1 = design is None
        # chk2 = smp_path is None
        # chk3 = dirs_ctl is None and dirs_exp is None
        # chk4 = all([i is None for i in [fq1, fq2, genome, outdir]])
        # # if config is None and not all(chk2):
        # if all([chk1, chk2, chk3, chk4]):
        #     sys.exit('required: --design, or --fq1, --fq2, --genome, \
        #         --outdir or --smp-path, --dirs-ctl, --dirs-exp')
        RNAseq(**args).run()


    def rnaseq2(self):
        """
        RNA-seq pipeline, simplify version
        """
        parser = add_rnaseq_args2()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args) # convert to dict
        RNAseqPipe(**args).run()


    def chipseq(self):
        """
        ChIPseq pipeline
        """
        parser = add_chipseq_args()
        args = parser.parse_args(sys.argv[2:])
        # help
        if len(sys.argv) < 3:
            parser.parse_args(['-h'])
        
        # main
        args = vars(args) # convert to dict
        # # check config or --fq1,--fq2,--genome,--outdir
        # config = args.get('config', None)
        # design = args.get('design', None)
        # fq1 = args.get('fq1', None)
        # fq2 = args.get('fq2', None)
        # genome = args.get('genome', None)
        # outdir = args.get('outdir', None)
        # chk1 = config is None
        # chk2 = design is None
        # chk3 = [i is None for i in [fq1, fq2, genome, outdir]]
        # 
        # # if config is None and not all(chk2):
        # if all([chk1, chk2, chk3]):
        #     sys.exit('required: --config, or --design, or --fq1, --fq2, --genome, --outdir')

        ChIPseq(**args).run()


    def cnr(self):
        """
        CUN&RUN pipeline
        """
        parser = add_cnr_args()
        args = parser.parse_args(sys.argv[2:])
        if len(sys.argv) < 3:
            parser.parse_args(['-h'])

        args = vars(args)
        CnR(**args).run()


    def fragsize(self):
        """
        Calculate the fragment size of PE alignment
        """
        parser = add_fragsize_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
        BamPEFragSize2(**args).run()


    def bam2bw(self):
        """
        Convert bam to bw files
        using: deeptools
        """
        parser = add_bam2bw_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
        # update scaleFactors
        bam_list = args.get('bam', None)
        sf = args.get('scaleFactor', 1.0)
        if len(sf) == 1:
            sf_list = sf * len(bam_list)
        elif not len(sf) == len(bam_list):
            raise ValueError('--scaleFactor not match bam files')
        else:
            sf_list = sf

        # for bam in args.get('bam', None):
        for i, bam in enumerate(bam_list):
            args_local = args.copy()
            args_local['bam'] = bam
            args_local['scaleFactor'] = sf_list[i]
            Bam2bw(**args_local).run()


    def bam2cor(self):
        """
        Calculate bam correlation
        using deeptools
        """
        parser = add_bam2cor_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
        args['make_plot'] = not args.get('no_plot', False)
        Bam2cor(**args).run()


    def peak2idr(self):
        """
        Calculate IDR for peak files
        using: idr
        """
        parser = add_peak2idr_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
        PeakIDR(**args).run()


    def bed2overlap(self):
        """
        Calculate IDR for peak files
        using: idr
        """
        parser = add_bed2overlap_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
        BedOverlap(**args).run()


    def sample(self):
        """
        Sub-sample fastq files
        head
        """
        parser = add_sample_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
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
        parser = add_trackhub_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args)
        TrackHub(**args).run()


def main():
    Hiseq()


if __name__ == '__main__':
    main()

