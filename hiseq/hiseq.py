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


class Hiseq(object):
    """The 1st-level of command, choose which sub-command to use
    qc, align, quant, peak, motif, report, ... (to be continued)

    """


    def __init__(self):
        parser = argparse.ArgumentParser(
            prog = 'hiseq',
            description = 'A collection of tools for HiSeq data',
            epilog = '',
            usage = ''' hiseq <command> [<args>]

    The most commonly used sub-commands are:

        qc        Basic quality control for fastq files, trim adapters, low-quality bases, ...
        align     Align fastq/a files to reference genome  
        quant     Count genes/features 
        peak      Call peaks using MACS2
        motif     Check motifs from a BED/fasta file    
        report    Create a report to the above commands
    '''
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
        '''Quality control, 
        cutadapt  : Trimming adapters from reads
        trim_ends : Python scripts to trim ends of reads, usually, the barcode sequences
        fastqc    : Create fastqc report the raw and clean files
        '''
        parser = argparse.ArgumentParser(
            description='Quality control for fastq files'
            )
        parser.add_argument('-i', '--fq1', 
            help='The fastq file(s) of SE or read1 file(s) of PE')
        # now that, we're inside a subcommand, ignore the first
        # TWO argvs, hiseq and qc
        args = parser.parse_args(sys.argv[2:])
        print('Running hiseq qc, fq1={}'.format(args.fq1))


    def align(self):
        '''Alignment
        Align fastq/a files to reference
        fq: SE or PE
        aligner: bowtie, bowtie2, STAR, bwa, hisat2, ...
        output: directory
        args: unique, multi, x-size, ...

        input
        output
        arguments.txt (pickle)
        '''
        parser = argparse.ArgumentParser(
            description='Align short reads to reference sequence')
        parser.add_argument('-n', '--aligner', 
            help='The aligner, [bowtie, bowtie2, STAR, bwa, hisat2]')
        parser.add_argument('-i', '--fq1',
            help='The fastq file(s) of SE or read1 file(s) of PE')
        parser.add_argument('-o', '--out',
            help='The directory to save the results')
        parser.add_argument('-I', '--fq2',
            help='The read2 of PE files')
        # now that, we're inside a subcommand, ignore the first
        # TWO argvs, hiseq and qc        
        args = parser.parse_args(sys.argv[2:])
        print('Running hiseq align, aligner={}, fq1={}, out={}'.format(
            args.aligner, args.fq1, args.out))

def main():
    # Hiseq()
    hiseq.qc.Trimmer()

        
if __name__ == '__main__':
    main()

