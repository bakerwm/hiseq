#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare bam files

correlation: PCA, cor
"""

import os
import pathlib
import argparse
from shutil import which
from hiseq.utils.file import file_exists, check_path, file_prefix
from hiseq.utils.utils import log, Config, run_shell_cmd, get_date, update_obj
from hiseq.utils.bam import Bam2cor, Bw2cor


def get_args():
    parser = argparse.ArgumentParser(description='hiseq bam2cor -i bam -o outdir')
    parser.add_argument('-i', '--bam-list', dest='bam_list', nargs='+', required=True,
        help='BAM files or bigWig files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-m', '--cor-method', default='pearson',
        choices=['pearson', 'spearman'],
        help='method to calculate correlation, default: [pearson]')
    parser.add_argument('-np', '--no-plot', dest='no_plot', action='store_false',
        help='do not make plots')
    parser.add_argument('-n', '--prefix', default=None,
        help='set the prefix for output files, default: [multibam]')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='number of threads, default: [1]')
    parser.add_argument('-b', '--binsize', default=500, type=int,
        help='set binSize for bigWig, default: [500]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Whether overwrite exists files')
    return parser


def main():
    args = vars(get_args().parse_args())
    if all([i.endswith('.bam') for i in args['bam_list']]):
        print('Bam2cor')
        Bam2cor(**args).run()
    elif all([i.endswith('.bigWig') for i in args['bam_list']]):
        print('Bw2cor')
        args['bw_list'] = args['bam_list'] # update
        Bw2cor(**args).run()
    else:
        log.error('no bam/bigWig files found')


if __name__ == '__main__':
    main()

#