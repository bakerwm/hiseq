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
from .utils.argsParser import *
from .qc.trimmer import Trimmer
from .align.alignment import Alignment
from .atac.atac2 import Atac

class Hiseq(object):
    """The 1st-level of command, choose which sub-command to use
    qc, align, quant, peak, motif, report, ... (to be continued)

    """
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog = 'hiseq',
            description = 'A collection of tools for HiSeq data',
            epilog = '',
            usage = """ hiseq <command> [<args>]

    The most commonly used sub-commands are:

        atac      ATACseq pipeline

        qc        Basic quality control for fastq files, trim adapters, low-quality bases, ...
        align     Align fastq/a files to reference genome  
        quant     Count genes/features 
        peak      Call peaks using MACS2
        motif     Check motifs from a BED/fasta file    
        report    Create a report to the above commands
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


    def qc(self):
        """Quality control, 
        cutadapt  : Trimming adapters from reads
        trim_ends : Python scripts to trim ends of reads, usually, the barcode sequences
        fastqc    : Create fastqc report the raw and clean files
        """
        parser = add_qc_args()
        args = parser.parse_args(sys.argv[2:])
        args = vars(args) # convert to dict
        
        # custom args
        fq1_list = args.get('fq1', None) # list
        outdir = args.get('outdir', None)

        ## iterate all fq1
        kwargs = args.copy()
        if kwargs['fq2'] is None:
            # SE mode
            for fq1 in fq1_list:
                kwargs['fq1'] = fq1
                kwargs['fq2'] = None
                kwargs['outdir'] = args.get('outdir', None)
                Trimmer(**kwargs).run()
        else:
            # PE mode
            if not len(fq1_list) == len(args['fq2']):
                log.error('-i, --fq2 not in the same length')

            for fq1, fq2 in zip(fq1_list, args['fq2']):
                kwargs['fq1'] = fq1
                kwargs['fq2'] = fq2
                kwargs['outdir'] = args.get('outdir', None)
                Trimmer(**kwargs).run()


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
        args = parser.parse_args(sys.argv[2:])
        print('hiseq peak')


    def motif(self):
        """
        quantify hiseq reads
        ...
        """
        parser = add_motif_args()
        args = parser.parse_args(sys.argv[2:])
        print('hiseq motif')


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

        RNAseq(**args).run()        


def main():
    Hiseq()

        
if __name__ == '__main__':
    main()

