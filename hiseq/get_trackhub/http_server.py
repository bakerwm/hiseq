#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Demo for trackhub, to create CompositeTrack

Stage:
1. Create simple compositeTrack (for few tracks)
2. Create two compositeTracks
3. Add subgroups for tracks (for large number of tracks)
4. Automatic recognize subtrack, subgroups, ...
"""


import os
import sys
import re
from hiseq.utils.utils import log, update_obj, Config, get_date, run_shell_cmd
from hiseq.utils.file import file_abspath


class HttpServer(object):
    """
    Deposite the local dirs/files to HTTP server

    Main purpose: convert path to http_url, on the localhost server.

    if http_server is 'localhost', it means the machine running this code
    is hosting httpd service.

    # to-do
    # 1. support ftp server (local, remote)
    # 2. support http server (remote)
    # require, server_ip, username, password, root_dir, ...
    # transfer files by SSH/ftp/...

    Convert the file_path to the http://path

    Parameters
    -----------
    remote_dir : str
        Absolute path to the directory/file in the http_root_dir

    http_root_dir : str
        Absolute path to the root_dir of http config, could be found in file
        `/etc/apache2/sites-available/000-browser.conf`, or other *.conf
        files, if you changed the Apache2 configuration. Check the Alias
        also.

    http_root_url : str
        The URL correspond to the http_root_dir.

    Example:

    >>> args = {
        's': '/data/public/hub.txt',
        'http_root_dir': '/data/public',
        'http_root_alias': '/public',
        'http_root_url': None,
        'is_https': False
        }
    >>> HttpServer(**args).to_url()


    ## Default Apache2 configuration

    On Ubuntu 18.04, could find these config in file:
    `/etc/apache2/sites-available/000-default.conf`, or your *.conf

    $grep -i 'DocumentRoot' /etc/apache2/sites-available/000-default.conf

        DocumentRoot /var/www/html

    it means, the http_root_dir is "/var/www/html"

    $ hostname -I
    1.2.3.4 # (fake ip, here)
    
    # in case multiple hosts, choose first one
    $ hostname -I 
    1.2.3.4  5.6.7.8

    you can find the http_root_ip is '1.2.3.4'

    ## For VirtualHost

    If VirtualHost has been build on the server, we can host the files on
    other directories beside the `DocumentRoot /var/www/html`

    for example:

    ```
    # file: /etc/apache2/sites-available/000-default.conf # or other .conf

    Alias /public /data/public
    <Directory "/data/public">
            Require all granted
    </Directory>
    ```

    The `http_root_url` for the virtualHost was: `hostname/publi`, and the
    correspond `http_root_path` was: `/data/public`
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        """Required Parameters

        - s : absolute path of the file/dir on the server
        - http_root_dir
        - http_root_url

        # for example:
        Alias /public /data/public

        http_root_alias = '/public'
        http_root_dir = '/data/public'
        http_root_url = http://ip/public
        """
        self.s = getattr(self, 's', '')
        self.http_root_dir = getattr(self, 'http_root_dir', '') # /data/public
        self.http_root_alias = getattr(self, 'http_root_alias', '') # public
        self.http_root_url = getattr(self, 'http_root_url', '')
        self.is_https = getattr(self, 'https', False)


    def init_files(self):
        """Initate the http_root_{dir,alias,url}

        Purpose
        fetch values for: http_root_{dir,alias,url}

        Option-1:
        http_root_dir is a valid URL,
          - parse http_root_alias, string after the {ip|domain},
            http://www.123.com{http_root_alias}

        Option-2:
        http_root_alias if a valid dirname,
          - parse the IP of 'localhost' HTTP server using shell command:
            `hostname -I`
          - determine HTTP server type: 'https' or 'http'
          - construct the http_root_url:
            {http}://{ip}{http_root_alias}

        """
        # absolute path
        self.s = file_abspath(self.s).rstrip('/')
        self.http_root_dir = file_abspath(self.http_root_dir).rstrip('/')

        # update http_root_{alias,url}
        if self.is_url(self.http_root_url):
            self.http_root_url = self.http_root_url.rstrip('/')
            url = re.sub('^((https?)|^(ftp))://', '', self.http_root_url)
            url_ip, url_dir = url.split('/', 1) #
            self.http_root_alias = '/' + url_dir
        elif len(self.http_root_alias) > 0:
            self.http_root_alias = self.http_root_alias.rstrip('/')
            url_ip = self.get_ip()
            url_http = 'https' if self.is_https else 'http'
            self.http_root_url = '{http}://{ip}{alias}'.format(
                http=url_http,
                ip=url_ip,
                alias=self.http_root_alias)
        else:
            log.error('check http_root_alias, http_root_url')
            raise ValueError(self.http_root_alias, self.http_root_url)

        # log message
        self.msg = '\n'.join([
            'Parameters',
            '{:>16s}: {}, {}'.format(
                's',
                isinstance(self.s, str),
                self.s),
            '{:>16s}: {}, {}'.format(
                'http_root_dir',
                isinstance(self.http_root_dir, str),
                self.http_root_dir),
            '{:>16s}: {}, {}'.format(
                'http_root_alias',
                isinstance(self.http_root_alias, str),
                self.http_root_alias),
            '{:>16s}: {}, {}'.format(
                'http_root_url',
                isinstance(self.http_root_url, str),
                self.http_root_url)
        ])

        # check http_root_dir and s
        # s contains http_root_dir
        if not self.s.startswith(self.http_root_dir):
            log.error('s and http_root_dir not match')
            raise ValueError(self.msg)

        # argument type
        if not all([isinstance(i, str) for i in [
            self.s,
            self.http_root_dir,
            self.http_root_url,
            self.http_root_alias]]):
            raise ValueError(self.msg)


    def get_ip(self):
        """Return the public ip address of the server

            - Ubuntu
            $ hostname -I
            1.2.3.4
            
            $ hostname -I
            1.2.3.4   5.6.7.8 
            # in case multiple hosts
        """
#         return run_shell_cmd('hostname -I')[1].strip()
        ip = run_shell_cmd('hostname -I')[1].strip()
        ip_list = ip.split(' ')
        if len(ip_list) > 1:
            ip = ip_list[0] # first
            log.warning('Multiple ip addresses found: {}, choose: [{}]'.format(
                ip_list, ip))
        # validate
        if self.is_ipv4(ip):
            return ip
        else:
            raise ValueError('not valid ip address found: {}'.format(ip))
        
    
    def is_ipv4(self, s):
        # validate ipv4
        def is_valid(i):
            try: 
                return str(int(i)) == i and 0 <= int(i) <= 255
            except: 
                return False
        # check
        return s.count(".") == 3 and all(is_valid(i) for i in s.split("."))


    def is_url(self, s):
        """
        Check if s is url:

        Parameters
        -----------
        s : str
            The url, http://, https://, ftp://

        - http://
        - https://
        - ftp://

        auto fix: www.abc.com -> http://www.abc.com
        127.0.0.1
        www.name.com
        http://127.0.0.1
        http:name.com
        """
        if not isinstance(s, str):
            return False

        # autofix: add 'http://' to url
        # www.abc.com
        # 127.0.0.1
        p1 = re.compile('(^www)|(^\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})',
            flags=re.IGNORECASE)
        if p1.search(s):
            s = 'http://' + s

        # search url
        p2 = re.compile('^((https?)|^ftp)://', flags=re.IGNORECASE)

        return p2.search(s)


    def to_url(self, s=None):
        """Convert s -> s_alias -> s_url

        Generate the URL of the hub_txt file, (validate the URL?)

        Parameters
        ----------
        s : str
            Path to the file/directory
        """
        if not isinstance(s, str):
            s = self.s

        if s.startswith(self.http_root_dir):
            s_alias = self.s[len(self.http_root_dir):]
            return self.http_root_url + s_alias
        else:
            log.error('s and http_root_dir not in the same dir')
            raise ValueError(self.msg)
