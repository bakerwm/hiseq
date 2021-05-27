#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
ATAC-seq pipeline: level-2 (generate report)
Generate report
"""

import os
import shutil
import hiseq
import argparse
from hiseq.utils.file import file_exists, check_path
from hiseq.utils.utils import log, update_obj, Config, get_date, run_shell_cmd, \
    read_hiseq


class AtacRp(object):
    """
    Parameters
    ----------
    x:  str
        Path to the hiseq project_dir

    >>> AtacRp(project_dir).run()
    >>> AtacRp(project_dir).report()

    $ python atac_rp.py -i prj_dir
    """
    def __init__(self, project_dir, **kwargs):
        self.project_dir = project_dir
        self.overwrite = kwargs.get('overwrite', False)
        self.init_args()


    def init_args(self):
        # check exists
        if not file_exists(self.project_dir):
            raise ValueError('project_dir not valid: {}'.format(
                self.project_dir))
        # check hiseq_dir
        a = read_hiseq(self.project_dir)
        if not a.is_hiseq:
            raise ValueError('project_dir not hiseq_dir: {}'.format(self.x))
        # check default files
        self.atacseq_type = 'atacseq_rt'
        self.report_dir = os.path.join(self.project_dir, 'report')
        self.config_toml = os.path.join(self.report_dir, 'config.toml')
        self.report_html = os.path.join(self.report_dir, 'HiSeq_report.html')
        self.report_stdout = os.path.join(self.report_dir, 'report.stdout')
        self.report_stderr = os.path.join(self.report_dir, 'report.stderr')
        check_path(self.report_dir)
        # save config
        Config().dump(self.__dict__, self.config_toml)


    def report(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        hiseq_report_R = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir,
            'HiSeq_report.html')
        cmd = ' '.join([
            '{}'.format(shutil.which('Rscript')),
            hiseq_report_R,
            self.project_dir,
            self.report_dir,
            '1>{}'.format(self.report_stdout),
            '2>{}'.format(self.report_stderr),
            ])
        # save command
        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # report_html
        if file_exists(atac_report_html) and not self.overwrite:
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            run_shell_cmd(cmd)


    def run(self):
        self.report()


def get_args():
    example = '\n'.join([
        'Examples:',
        '$ python atac_rp.py -i project_dir',
    ])
    parser = argparse.ArgumentParser(
        prog='atac_rp',
        description='Generate hiseq report',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--project-dir', dest='project_dir',
        required=True,
        help='Directory of hiseq project')
    parser.add_argument('-r', '--overwrite', action='store_true',
        help='Overwrite the exists file.')
    return parser


def main():
    args = vars(get_args().parse_args())
    AtacRp(**args).run()


if __name__ == '__main__':
    main()

#