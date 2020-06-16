#!/usr/bin/env python3
#-*- coding: utf8 -*-

"""
Download files from website:
1. download fastq files from Novogene.com, using linuxnc

    lnd login -u <user> -p <passwd>
    lnd status
    lnd list <remote/path>
    lnd cp -d <remote/path/> <local/path/>

2. download SRA files from NCBI/SRA, (GSE/SRA/SRP/...)

3. download supp data from NCBI/GEO, (GSE/GSM)
"""

import os
import sys
import re
import logging
import argparse
import shutil
import subprocess
# from hiseq.utils.helper import run_shell_cmd

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)

log.setLevel('INFO')


def run_shell_cmd(cmd):
    """This command is from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py

    save log to file
    """
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid,
        pgid,
        rc,
        stderr.strip(),
        stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return (stdout.strip('\n'), stderr.strip('\n'))


def get_lnd_args():
    parser = argparse.ArgumentParser(description='download files')
    parser.add_argument('-i', '--config', required=True,
        help='config file for download, save id/passwd/path/...')
    parser.add_argument('-o', '--local-path', dest='local_path', required=True,
        help='path to save the files')
    parser.add_argument('-y', '--yes', action='store_true',
        help='Do not ask for confirmation')
    args = parser.parse_args(sys.argv[2:])
    return vars(args)


class Linuxnd(object):
    """
    command: lnd (linux)
    download sequencing data
     
    for linux: http://data-deliver.novogene.com/linuxnd.zip
    for Windows: http://data-deliver.novogene.com/DataDeliver-win32-x64.zip
    for MacOS: http://data-deliver.novogene.com/macnd.zip

    lnd login -u <user> -p <passwd>
    lnd status
    lnd list <remote/path>
    lnd cp -d <remote/path/> <local/path/>
    
    lnd -i {} -o {} -t {}
    """
    def __init__(self, **kwargs):
        """
        Required args
        user
        passwd
        remote_dir
        local_dir
        """
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.init_args()
        self.load_config()


    def init_args(self):
        self.yes = getattr(self, 'yes', False)
        self.local_path = getattr(self, 'local_path', os.path.join(os.getcwd(), 'from_illumina'))

        # cmd
        self.lnd = shutil.which('lnd')
        if self.lnd is None:
            raise ValueError('''
                lnd, command not found in PATH, try to download from: 
                http://data-deliver.novogene.com/linuxnd.zip, 
                add it to your PATH, and re-run this command''')


    def load_config(self):
        # config
        if not hasattr(self, 'config'):
            raise ValueError('''
                config.txt, file not found, expected, 
                id: ***
                passwd: ***
                path: ***''')

        with open(self.config) as r:
            s = r.read()

        # user/id
        self.user = None
        self.password = None
        self.remote_path = None
        p1 = re.compile('id:\s?([\w|\-|\.]+)') # user
        p2 = re.compile('password:\s?([\w|\!|\@|\#|\$|\%]+)') # pasword
        p3 = re.compile('path:\s?([\w|\-|\/|\:]+)') # path
        if p1.search(s):
            self.user = p1.search(s).group(1)
        if p2.search(s):
            self.password = p2.search(s).group(1)
        if p3.search(s):
            self.remote_path = p3.search(s).group(1)

        # check
        msg = '\n'.join([
            '{:>10}: {}'.format('user', self.user),
            '{:>10}: {}'.format('password', self.password),
            '{:>10}: {}'.format('path', self.remote_path)])

        # status
        if not all([self.user, self.password, self.remote_path]):
            print(msg)
            raise ValueError('checkout the config: {}'.format(self.config))

        # add prefix: oss://
        self.remote_path = 'oss://' + self.remote_path #


    def login(self):
        cmd = ' '.join([self.lnd, 'login', '-u', self.user, '-p', self.password])
        stdout, stderr = run_shell_cmd(cmd)


    def status(self):
        cmd = ' '.join([self.lnd, 'status'])
        stdout, stderr = run_shell_cmd(cmd)
        if stderr.encode('utf-8') == '当前没有登录'.encode('utf-8'):
            log.info('login required.')
            # self.login()


    def listfile(self):
        cmd = ' '.join([self.lnd, 'list', self.remote_path])
        return run_shell_cmd(cmd)


    def download(self):
        # check status
        self.status()
        self.login()
        self.listfile()

        # ask user, whether to download files
        log.info('Saving data to: {}'.format(self.local_path))
        if self.yes:
            chk = 'Y'
        else:
            chk = input("Go on downloading files, [Y|n]: ")
        # check
        if chk == 'Y':
            log.info('Start downloading ...')
            cmd = ' '.join([self.lnd, 'cp', '-d', self.remote_path, self.local_path])
            run_shell_cmd(cmd)
        else:
            log.info('Downloading skipped!')


def main():
    parser = argparse.ArgumentParser(
        prog = 'download',
        description = 'A collection of tools for HiSeq data',
        epilog = '',
        usage = """ download <command> <args>

        lnd        Download HiSeq data from Novogene cloud
        sra        Download SRA data
        geo        Download GEO supp data
    """)
    parser.add_argument('command', help='Subcommand to run')
    args = parser.parse_args(sys.argv[1:2])
    
    if args.command == 'lnd':
        args = get_lnd_args()
        Linuxnd(**args).download()
    elif args.command == 'sra':
        pass
    elif args.command == 'geo':
        pass
    else:
        print('Unrecognized command: {}'.format(args.command))
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()