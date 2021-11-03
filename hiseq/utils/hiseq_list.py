#!/usr/bin/env python3
# -*- encoding: utf8 -*- 

"""
list hiseq dir, file

hiseq list -n name -t _rx [dir1, ...]
"""

import os
import argparse
from hiseq.utils.utils import list_hiseq_dir, list_hiseq_file, is_hiseq_dir
from hiseq.utils.file import list_dir, file_abspath


def list_dirs(x):
    # current, 
    # subdir
    if isinstance(x, str):
        if os.path.isdir(x):
            out = [i for i in list_dir(x, include_dir=True) if os.path.isdir(i)]
            out.append(x)
    elif isinstance(x, list):
        out = [j for i in x for j in list_dirs(i)]
    else:
        out = []
    return out
        

def hiseq_list(dirs, **kwargs):
    """
    list hiseq files, dirs
    """
    args = {
        'name': None,
        'hiseq_type': 'auto',
        'recrusive': False,
        'abspath': False,
    }
    args.update(kwargs)
    if isinstance(dirs, str):
        if not os.path.isdir(dirs):
            out = []
        elif isinstance(args['name'], str):
            out = list_hiseq_file(dirs, args['name'], args['hiseq_type'])
        else:
            try:
                out = list_hiseq_dir(dirs, args['hiseq_type'])
            except:
                out = []
                print('!A-3', dirs)
        # convert to list
        if isinstance(out, str):
            out = [out]
    elif isinstance(dirs, list):
        out = [j for i in dirs for j in hiseq_list(i, **kwargs)]
    else:
        out = []
    # fix None
    if out is None:
        out = []
    if args['abspath']:
        out = file_abspath(out)
    # return
    return sorted(list(set(out)))


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. list the hiseq dirs',
        '$ hiseq hiseq_list path',
        ' ',
        '2. list the hiseq files',
        '$ hiseq hiseq_list -n smp_name path'        
    ])
    parser = argparse.ArgumentParser(
        prog='hiseq_list',
        description='list hiseq dir/file',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('dirs', nargs='+',
        help='list of directories')
    parser.add_argument('-n', '--name', default=None, 
        help='name of the hiseq file, if not return the directory')
    parser.add_argument('-t', '--hiseq-type', dest='hiseq_type', default='auto',
        help='the hiseq type, default: [auto]')
    parser.add_argument('-r', '--recursive', action='store_true',
        help='search the directories recursively')
    parser.add_argument('-a', '--abspath', action='store_true',
        help='return the absolute path')
    return parser


def main():
    args = vars(get_args().parse_args())
    dirs = args.pop('dirs', None)
    dirs = list_dirs(dirs) # update, root/subdir
    out = hiseq_list(dirs, **args)
    print('\n'.join(out))


if __name__ == '__main__':
    main()


#