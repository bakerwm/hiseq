#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Check overlap between bed files
"""

import os
import pathlib
import argparse
from shutil import which
from hiseq.utils.file import file_exists, check_path, file_prefix
from hiseq.utils.utils import log, Config, run_shell_cmd, get_date, update_obj
from hiseq.utils.bed import BedOverlap


def get_args():
    """
    required:
            'peak': None,
            'outdir': str(pathlib.Path.cwd()),
            'flag': False,
            'prefix': None,
            'overwrite': False
    """
    parser = argparse.ArgumentParser(description='hiseq bed2overlap -i peak -o outdir')
    parser.add_argument('-i', '--peak', nargs='+', required=True,
        help='peak files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-n', '--prefix', default=None,
        help='set the prefix for output files, default: [multibam]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Whether overwrite exists files')
    return parser

def main():
    args = vars(get_args().parse_args())
    BedOverlap(**args).run()


if __name__ == '__main__':
    main()

#