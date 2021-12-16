#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Using featureCounts calculate counts on intervals (genes)
"""

import os
import sys
import pathlib
from hiseq.utils.utils import log
from hiseq.utils.featurecounts import LibStrand, FeatureCounts


def fc_quant(gtf, bam, **kwargs):
    """
    Parameters
    ---------
    x: str
        The path to RnaseqR1() dir, hiseq_type=rnaseq_r1

    Run FeatureCounts for the bam file
    """
    args = {        
        'gtf': gtf,
        'bam_list': bam,
        'strandness': 'sens', # sens, anti, both
        'outdir': str(pathlib.Path.cwd()),
        'threads': 1,
        'overwrite': False
    }
    args.update(kwargs)
    if not all(file_exists([gtf, bam])): 
        msg = 'quant() skipped, file not exists, bam: {}, gtf: {}'.format(
            file_exists(bam), file_exists(gtf)
        )
        log.error(msg)
        return None
    # strandness
    if args['strandness'] == 'both':
        args['strandness'] = 0
    else:
        # fc_no = FeatureCounts(**args).run()
        # guess sense strand
        s = LibStrand(bam=bam, gtf=gtf, clean_up=True) # sense_strand, anti_strand
        if args['strandness'] == 'sens':
            args['strandness'] = s.sens_strand
        elif args['strandness'] == 'anti':
            args['strandness'] = s.anti_strand
        else:
            args['strandness'] = 0 #
    fc = FeatureCounts(**args_sense).run()
    

def get_args():
    example = '\n'.join([
        'Examples:',
        '1. support fastq input',
        '$ python quant.py -a gene.gtf -b [bam1, ...] -p 8 -s sens'
    ])    
    parser = argparse.ArgumentParser(
        prog='quant',
        description='python quant.py',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser = argparse.ArgumentParser(
        description='featureCounts')
    parser.add_argument('-a', '--gtf', default=None, required=True,
        help='The gtf file for quantification')
    parser.add_argument('-b', '--bam', nargs='+', required=True,
        help='bam files')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads for each job, default [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    return parser


def main():
    args = vars(get_args().parse_args())
    bam = args.pop('bam', None)
    gtf = args.pop('gtf', None)
    fc_quant(gtf, bam, **args).run()


if __name__ == '__main__':
    main()

#

