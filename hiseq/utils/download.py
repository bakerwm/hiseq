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
import math
import logging
import argparse
import shutil
import subprocess
import pathlib
import urllib
from .parallel_download import Downloader
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


def get_lnd_args(argv):
    parser = argparse.ArgumentParser(description='download files')
    parser.add_argument('-i', '--config', required=True,
        help='config file for download, save id/passwd/path/...')
    parser.add_argument('-o', '--local-path', dest='local_path', required=True,
        help='path to save the files')
    parser.add_argument('-y', '--yes', action='store_true',
        help='Do not ask for confirmation')
    args = parser.parse_args(argv) # sys.argv[2:]
    return vars(args)


def get_oss_args(argv):
    parser = argparse.ArgumentParser(description='download by OSS url')
    parser.add_argument('-i', '--url', required=True,
        help='url or a file saved urls, one per line')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True,
        help='path to save the files')
    parser.add_argument('-t', '--threads', type=int, default=4,
        help='Number of threads, to download file. default: [4]')
    parser.add_argument('-n', '--dry-run', dest='dry_run', action='store_true',
        help='List files, do not download')
    args = parser.parse_args(argv) # sys.argv[2:]
    return vars(args)
    

class OSSDownload(object):
    """
    Download aliyun OSS object, by URL with signature
    
    Example:    
    args = {
        'url': 'http://*',
        'outdir': 'download',
        'threads': 4,
    }
    
    OSSDownload(**args).run()
    """
    def __init__(self, **kwargs):
        """
        URL:
        http://{bucket}.oss-cn-shenzhen.aliyuncs.com/{object}?Expires={}&OSSAccessKeyId={}&Signature={}
        """
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.init_args()
        self.url_list = self.load_urls()


    def init_args(self):
        outdir = os.path.join(str(pathlib.Path.cwd()), 'from_illumina')
        self.outdir = getattr(self, 'outdir', outdir)
        self.threads = getattr(self, 'threads', 4)
        self.url = getattr(self, 'url', None)
        self.dry_run = getattr(self, 'dry_run', False)
        if self.threads not in range(1, 20):
            print('--threads {} not valid, set [4]'.format(self.threads))
            self.threads = 4
        if not isinstance(self.url, str):
            raise ValueError('--url not valid, {}'.format(self.url))
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)


    def load_urls(self):
        """
        urls in file (one per line) or url
        """
        url_list = []
        if os.path.isfile(self.url):
            try:
                with open(self.url) as r:
                    for line in r:
                        line = line.lstrip().rstrip()
                        if line.startswith('#'):
                            continue
                        if len(line) > 10: # at least
                            url_list.append(line)
            except IOError as e:
                print(e)
        elif isinstance(self.url, str):
            url_list.append(self.url)
        else:
            raise ValueError('--url not valid, {}'.format(self.url))
        if len(url_list) == 0:
            raise ValueError('--url, no urls detected')
        return url_list    
        
                        
    def get_file_name(self, url):
        """Get the file name from URL
        for OSS downloader
        example:
        http://{bucket}.oss-cn-shenzhen.aliyuncs.com/{object}?Expires={}&OSSAccessKeyId={}&Signature={}
        """
        s = re.sub('http:.*./', '', url)
        s = re.sub('kefu\%2F', '', s)
        s = re.sub('\?Expires.*', '', s)
        s = re.sub('[^A-Za-z0-9\.\-\_]', '_', s)
        return s


    def make_a_request(self, url, referer=None):
        """
        create a http request
        :param url:
        :param referer:
        :return: urllib2.Request
        """
        o = urllib.parse.urlparse(url)
        headers = {
            'Accept': '*/*',
            'Accept-Language': 'zh-CN,zh;q=0.8,en;q=0.6,zh-TW;q=0.4',
            'Connection': 'keep-alive',
            'DNT': '1',
            'Host': o.netloc,
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/40.0.2214.45 Safari/537.36'
        }
        if referer:
            headers['Referer'] = referer
        return urllib.request.Request(url, headers=headers)
    
    
    def get_url_length(self, url, referer=None):
        request = self.make_a_request(url, referer)
        # request.get_method = lambda: 'HEAD'
        try:
            response = urllib.request.urlopen(request)
            content_length = response.getheader('Content-Length')
            if content_length:
                return int(content_length)
        except urllib.error.URLError as e:
            print(e)
        return -1


    def get_file_size(self, url, retry_times=1, referer=None):
        for i in range(0, retry_times):
            try:
                return self.get_url_length(url, referer)
            except:
                print("can't get size of {}, try {}".format(url, i))
        raise Exception("can't get size of {}, try {}".format(url, retry_times))

        
    def readable_size(self, s):
        """
        s is the size in bytes
        Convert to human readable size: K/M/G/T
        """
        if isinstance(s, int):
            if s < 1024: # B
                out = s
            elif s < math.pow(1024, 2): # KB
                out = '{:.1f}K'.format(s/1024)
            elif s < math.pow(1024, 3): # MB
                out = '{:.1f}M'.format(s/math.pow(1024, 2))
            elif s < math.pow(1024, 4): # GB
                out = '{:.1f}G'.format(s/math.pow(1024, 3))
            elif s < match.pow(1024, 5): # TB
                out = '{:.1f}T'.format(s/math.pow(1024, 4))
            elif s < match.pow(1024, 6): # PB
                out = '{:.1f}P'.format(s/math.pow(1024, 5))
            else:
                out = s
        else:
            print('illegal s={}, expect int'.format(s))
            out = s
        return out
            
        
    def list_files(self):
        """
        Only get the file size, do not download files
        """
        lines = []
        lines.append('{:>7}\t{:<30}'.format('Size', 'Filename'))
        for url in self.url_list:
            fname = self.get_file_name(url)
            target_file = os.path.join(self.outdir, fname)
            s = self.get_file_size(url)
            sa = self.readable_size(s)
            lines.append('{:>7}\t{:<30}'.format(sa, target_file))
        msg = '\n'.join(lines)
        print(msg)
    
    
    def download(self):
        for url in self.url_list:
            fname = self.get_file_name(url)
            target_file = os.path.join(self.outdir, fname)
            obj = Downloader(url, self.threads, target_file)
            obj.start_download()
            print(obj.get_metadata())
            print(obj.get_remote_crc32c())
            print(obj.get_downloaded_crc32c())
            
    
    def run(self):
        if self.dry_run:
            self.list_files()
        else:
            self.download()
    
    
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


def main(argv):
    parser = argparse.ArgumentParser(
        prog = 'download',
        description = 'A collection of tools for HiSeq data',
        epilog = '',
        usage = """ download <command> <args>
        oss        Download aliyun OSS, by URL with signature
        lnd        Download HiSeq data from Novogene cloud
        sra        Download SRA data
        geo        Download GEO supp data
    """)
    parser.add_argument('command', help='Subcommand to run')
    args = parser.parse_args(argv[1:2])
    
    if args.command == 'lnd':
        args = get_lnd_args(argv[2:])
        Linuxnd(**args).download()
    elif args.command == 'oss':
        args = get_oss_args(argv[2:])
        OSSDownload(**args).run()
    elif args.command == 'sra':
        pass
    elif args.command == 'geo':
        pass
    else:
        print('Unrecognized command: {}'.format(args.command))
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main(sys.argv)