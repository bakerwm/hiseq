#!/usr/bin/env python3

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
import signal
import time
import yaml
import logging
import trackhub
from urllib.request import urlopen
from bs4 import BeautifulSoup
from hiseq.utils.helper import *


# logging.basicConfig(
#     format='[%(asctime)s %(levelname)s] %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     stream=sys.stdout)
# log = logging.getLogger(__name__)
# log.setLevel('INFO')

## run
def get_ticks():
    """Returns ticks.
       Mark the time stamp
        - Python3: Use time.perf_counter().
    """
    return getattr(time, 'perf_counter', getattr(time, 'time'))()


def run_shell_cmd(cmd):
    """
    This command is adapted from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_lib_common.py

    Run shell command
    capture the output, code
    """
    # mark time
    t0 = get_ticks() # start_time, to run process
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to save mode bash
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        start_new_session=True) # to make child process
    # process id
    pid = p.pid
    pgid = os.getpgid(pid)

    # log
    msg_log = 'run_shell_cmd(): PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd)
    log.info(msg_log)

    # run
    out, err = p.communicate(cmd)
    rc = p.returncode

    # mark time
    t1 = get_ticks()
    time_dur = '{:.1f}'.format(t1 - t0) # seconds

    # err msg
    msg_out = 'PID={}, PGID={}, RC={}, DURATION_SEC={} \nSTDOUT={} \nSTDERR={}'.format(
        pid, pgid, rc, time_dur, out, err)

    # check status
    if rc:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(msg_out)
    else:
        log.info(msg_out)

    return (rc, out, err)


## dict to yaml
def dict2yaml(d, s):
    """
    convert dict to yaml, save to file
    
    Parameters
    -----------
    d : dict
        A dict saving objects 
    s : str
        path to a file, will saving the dict as YAML format
    
    dimX:
        name: stage,
        label: stage,
        mapping:
            1h: 1h,
            2h: 2h
    dimY:
        name: gene,
        label: gene,
        mapping:
            geneA: geneA,
            geneB: geneB

    Dict format:

    {
        'dimX': {
            'name': 'stage,',
            'label': 'stage,',
            'mapping': {
                '1h': '1h,', 
                '2h': '2h'
            }
        },
        'dimY': {
            'name': 'gene,',
            'label': 'gene,',
            'mapping': {
                'geneA': 'geneA,', 
                'geneB': 'geneB'
            }
        }
    }
    """
    with open(s, 'w') as w:
        try:
            yaml.safe_dump(d, w, default_flow_style=False, allow_unicode=True, 
                encoding='utf-8', indent=4)
        except yaml.YAMLError as exc:
            log.error(exc)


def yaml2dict(s):
    """
    load file (.yaml), to dict

    YAML format:

    dimX:
        name: stage,
        label: stage,
        mapping:
            1h: 1h,
            2h: 2h
    dimY:
        name: gene,
        label: gene,
        mapping:
            geneA: geneA,
            geneB: geneB

    Dict format:

    {
        'dimX': {
            'name': 'stage,',
            'label': 'stage,',
            'mapping': {
                '1h': '1h,', 
                '2h': '2h'
            }
        },
        'dimY': {
            'name': 'gene,',
            'label': 'gene,',
            'mapping': {
                'geneA': 'geneA,', 
                'geneB': 'geneB'
            }
        }
    }
    """
    with open(s, 'r') as r:
        try:
            return yaml.safe_load(r)
        except yaml.YAMLError as exc:
            log.error(exc)


# Stage4. Auto recognize composite tracks, subgroups
#
# 1. if subgroups.json exists, the files in the folder should organize into one
#    compositeTrack()
#
# 2. According to TrackFile, add signal_view or region_view, or both
#
# 3. Create, Init compositeTrack, based on subgroups


################################################################################
## functions for Track Hub ##
class HubUrl(object):
    """For hub.txt URL

    Purpose
    -------
    1. Check, validate hub_url 
    2. convert hub_url to trackhub_url

    Parameters
    ----------
    hub_url : str
        URL direct to the hub.txt file, that is open accessiable

    genome : str or None
        The genome relase of the tracks, if [None], parsing this value from
        hub_url, in this way: hub_url -> hub.txt -> genomes.txt -> genome.
        default: [None]

    mirror : str
        The UCSC genome browser mirror, one of ['usa', 'euro', 'asia'] or 
        the URL of other mirror, default: [usa]

    position : str or None
        The coordinates of the position to show in default in trackhub_url,
        in this format: 'chr1:2000-4000', default: [None]

    validate_url : bool
        Validate the hub_url using hubCheck command line tool. default: [False]

    Example:
    >>> a = HubUrl(hub_url=hub_url, genome='dm6', 
        position='chr2:10,000,000-12,000,000').trackhub_url()
    >>> print(a)

    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        Parameters
        -----------
        hub_url : str
            URL to the hub.txt file, open accessiable
        """
        self.hub_url = getattr(self, 'hub_url', None)
        self.genome = getattr(self, 'genome', None) # specify the genome
        self.mirror = getattr(self, 'mirror', 'usa')
        self.position = getattr(self, 'position', None)
        self.validate_url = getattr(self, 'validate_url', False)

        if not isinstance(self.hub_url, str):
            raise ValueError('hub_url failed, str expected, got {}'.format(
                type(self.hub_url).__name__))

        if not self.is_url(self.hub_url):
            raise ValueError('hub_url failed, not a URL: {}'.format(
                self.hub_url))

        if self.validate_url:
            if self.validate_hub(): # return code: 0=ok
                raise ValueError('hub_url failed, hubCheck failed: {}'.format(
                    self.hub_url))


    def validate_hub(self):
        """
        Using the hubCheck tool, to check the url
        """
        if self.is_url(self.hub_url):
            cmd = ' '.join([
                shutil.which('hubCheck'),
                '-noTracks',
                self.hub_url
                ])

            # run
            try:
                return run_shell_cmd(cmd)[0]
            except:
                log.error('hub_url failed, not a valid hub_url: {}'.format(
                    self.hub_url))
                return 1 # fail


    def is_url(self, s):
        """
        Check if hub_url is url:

        Parameters
        -----------
        s : str or None
            The hub.txt url

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


    def read_url(self, s):
        """
        Read the page of url, save as dict

        Parameters
        -----------
        s : str or None
            The hub.txt url        
        """
        d = {}
        try:
            log.info('Parsing url: {}'.format(s))
            html = urlopen(s).read()
            text = BeautifulSoup(html, features='html.parser')
            x = re.split('\s', text.text)
            xit = iter(x)
            d = dict(zip(xit, xit))
        except:
            log.error('reading url failed, check your internet connection:\
                {}'.format(s))

        return d


    def is_hub_url(self, s):
        """s is a URL of hub.txt file

        Parsing 'hub', 'genomesFile' in the file
        """
        args_hub = self.read_url(s)
        return all([i in args_hub for i in 
            ['hub', 'shortLabel', 'longLabel', 'genomesFile']])


    def is_genome_url(self, s):
        """s is a URL of genomes.txt file
        Parsing 'genome', 'trackDb' attributes
        """
        args_hub = self.read_url(s)
        return all([i in args_hub for i in ['genome', 'trackDb']])


    def get_genome(self):
        """Fetch the genome release

        Option-1:
        self.genome

        Option-2:
        parsing(hub_url) -> genomes.txt -> genome (the 1st one)
        """
        if isinstance(self.genome, str):
            genome = self.genome
        else:
            log.info('Parsing the [genome] from URL')
            # Parsing the hub_url webpage (plain text)
            args_hub = self.read_url(self.hub_url)

            # parsing the genome (hg19, hg38, dm6, ...)
            # construct the genome_url
            genome_fname = args_hub.get('genomesFile', 'genomes.txt')
            genome_url = os.path.join(os.path.dirname(self.hub_url), 
                genome_fname)
            args_genome = self.read_url(genome_url)

            genome = args_genome.get('genome', None)

            if not isinstance(genome, str):
                raise ValueError('genome failed, illegal file: {}'.format(
                    genome_url))

            log.info('Got [genome={}]'.format(genome))

        return genome


    def init_mirror(self, s):
        """Iniate the mirror of UCSC Genome Browser

        Parameters
        -----------
        s : str
            name of the mirrors: ['usa', 'euro', 'asia'], or the url to an 
            mirror

        return the url of mirror
        """
        # common mirrors
        mirrors = {
            'usa': 'http://genome.ucsc.edu',
            'asia': 'http://genome-asia.ucsc.edu',
            'euro': 'http://genome-euro.ucsc.edu'}

        # determine the mirror
        if self.is_url(s):
            mirror = s
        elif s in mirrors:
            mirror = mirrors.get(s, mirrors['usa'])
        else:
            mirror = mirrors['usa']
            log.warning((
                'mirror, illegal value, http://, or'
                '[usa|asia|euro] expect, got {}, set mirror=usa').format(s))

        return mirror


    def format_position(self, position=None):
        """Convert position to url format
        Input: chr1:1-1000
        Output: &position=chr1%3A1-1000

        Parameters
        -----------
        position : str or None
            The coordinates of position, format: chr:start-end

        see: UCSC Documentation
        https://genome.ucsc.edu/goldenPath/help/customTrack.html#optParams
        """
        if not isinstance(position, str):
            posotion = self.position

        fmt_pos = '' # formated position

        if isinstance(position, str):
            position = re.sub('[^\w:-]', '', position) # sanitize chr
            p = re.compile('^(\w+):([\d,]+)-([\d,]+)$') # chr1:200-400
            m = p.search(position)

            # output format: chr%3A{}-{}
            if m:
                fmt_pos = '&position={chr}%3A{start}-{end}'.format(
                    chr=m.group(1),
                    start=m.group(2), 
                    end=m.group(3))

        return fmt_pos


    def trackhub_url(self, position=None):
        """Convert hub_url to trackhub_url

        It is user-friendly to share trackHub using trackhub_url, open the url
        will direct user to the trackHub

        If hub_url only, user need to 
        1. open UCSC genome browser: 'http://genome.ucsc.edu', 
        2. choose 'My Data' -> 'Track Hubs', 
        3. paste hub_url to the box next to 'URL:', under 'My Hubs' Tab, 
        4. click 'Add Hub'
    
        For example:
        hub_url: 
            - http://ip/hub.txt
        trackhub_url: 
            - http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://ip/hub.txt 

        Parameters
        -----------
        mirror : str
        The UCSC genome browser mirror, one of ['usa', 'euro', 'asia'] or 
        the URL of other mirror, default: [usa]

        usa:  http://genome.ucsc.edu
        asia: http://genome-asia.ucsc.edu
        euro: http://genome-euro.ucsc.edu
        custome: ''
        """
        if not isinstance(position, str):
            position = self.position

        # constructure trackhub_url
        mirror = self.init_mirror(self.mirror)
        genome = self.get_genome()
        hub_url = self.hub_url

        fmt_pos = self.format_position(position)

        return (
            '{mirror}/cgi-bin/hgTracks?db={genome}'
            '{position}&hubUrl={hub_url}').format(
            mirror=mirror, genome=genome, position=fmt_pos, hub_url=hub_url)


class HttpServer(object):
    """Deposite the local dirs/files to HTTP server

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
            log.errorr('s and http_root_dir not match')
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
        """
        return run_shell_cmd('hostname -I')[1].strip()


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


class TrackFile(object):
    """Extract information from file

    Parameters
    ----------
    s : str
        The path to a track file, 

    subgroups : dict
        The subgroups structure in dict format, including dimX, dimY, ...

    colorDim : str
        Assign colors to the tracks based on which dimension, dimX, dimY, ...
        default: [dimX]

    colorPal : int or str
        choose the group of colors, 1 to N, or the name of fish
        candidate list: [1, 2, 3, 'Scarus_hoefleri', 'Gramma_loreto',
        'Centropyge_loricula'], default: [1]

    label_rm_list : list
        A list of strings, remove from the 'short_label', to make sure the 
        length of short label less than 17 characters. default: []
    """
    def __init__(self, s, **kwargs):
        """
        Parse attributes for track file

        Parameter
        ---------
        s : str
            path to a file
        """
        for k, v in kwargs.items():
            if not hasattr(self, k):
                setattr(self, k, v)
        self.s = s
        self.init_args()


    def init_args(self):
        """Default values

        subgroups support: 
            dimX: 
            dimY:
            dimA:
            ... (A-Z)
        self.subgroups = self.get_subgroups() 
        """
        # default args
        self.subgroups = getattr(self, 'subgroups', {}) # pass groups
        self.colorDim = getattr(self, 'colorDim', 'dimX') # assign colors
        self.label_rm_list = getattr(self, 'label_rm_list',  []) # rm
        self.colorPal = 1 # color palette, option: 1, 2, 3
        self.fname, self.fext = os.path.splitext(os.path.basename(self.s))
        self.fname = self.sanitize(self.fname)
        self.label = self.sanitize(self.fname, self.label_rm_list)
        self.ftype = self.filetype(self.fext)
        # self.fcolor = self.pick_color(self.fname)        
        self.subgroup = self.get_subgroup() # blank dict
        self.color = self.get_color() # '255,0,0'

        self.track = self.signal_track() if self.ftype in \
            ['bigWig', 'bedGraph'] else self.region_track()


    def is_track_file(self):
        """trackhub supported files: 
        - bigWig
        - bigBed
        - bedGraph
        """
        return any([self.is_bw(), self.is_bb(), self.is_bg()])


    def is_bw(self):
        """bigWig files:
        - bigWig 
        - bigwig
        - bw
        """
        return self.ftype == 'bigWig'


    def is_bb(self):
        """bigBed files:
        - bigBed
        - bigbed
        - bb
        """
        return self.ftype == 'bigBed'


    def is_bg(self):
        """bedGraph files:
        - bedGraph
        - bedgraph
        - bg
        """
        return self.ftype == 'bedGraph'


    def filetype(self, s):
        """Determine the type of file
        
        Parameters
        ----------
        s : str
            path to a file

        Supported format list:
        - bigWig: bigWig, bigwig, bw
        - bigBed: bigBed, bigbed, bb
        - bedGraph: bedGraph, bedgraph, bg
        """
        s = s.lstrip('.').lower()

        # predefined ftypes
        ftypes = {
            'bw': 'bigWig',
            'bigwig': 'bigWig',
            'bb': 'bigBed', 
            'bigbed': 'bigBed',
            'bg': 'bedGraph',
            'bdg': 'bedGraph'
        }

        return ftypes.get(s, None)


    def sanitize(self, s, remove_list=None):
        """
        Sanitize a string, for shortLabel, letters, numbers allowed only [A-Za-z0-9]
        
        Convert spaces to unserscore, remove other chaaracters

        Parameters
        -----------
        s : str
            String to sanitize

        remove_list : list or None
            if list, remove the items from string; if None, skipped
        """
        # remove non-characters:
        s1 = re.sub('[^A-Za-z0-9_ ]', '', s)

        # replace spaces to underscore
        s1 = re.sub('\s+', '_', s1)

        # remove specific strings
        if isinstance(remove_list, list) and len(remove_list) > 0:
            for i in remove_list:
                if not isinstance(i, str):
                    continue
                s1 = re.sub(i, '', s1, flags=re.IGNORECASE)

        return s1


    def get_subgroup(self):
        """Assign subgroups for file

        Searching 'mapping()' of each subgroup in filename
    
        return the 1st match

        Format for subgroups
        {
            'dimX': {
                'name': 'stage',
                'label': 'stage',
                'mapping': {
                    '2h': '2hrs',
                    '4h': '4hrs'
                }
            },
            'dimY': {
                'name': 'gene',
                'label': 'Gene',
                'maping': {
                    'geneA': 'GeneA',
                    'geneB': 'GeneA'
                }
            }
        }
    
        # output format:
        # subgroup = {'stage': '2h', 'gene': 'geneA'}
        """
        subgroup = {}
        if isinstance(self.subgroups, dict) and len(self.subgroups) > 0:
            for d, v in self.subgroups.items():
                # iterate, dimensions
                # dimX, dimY, dimA, ...
                d_name = v.get('name', None)
                d_mapping = v.get('mapping', None)
                for m, n in d_mapping.items():
                    # iterate, groups
                    # 2h, 4h, ... (in dimX['mapping'])
                    if m in self.fname:
                        subgroup[d_name] = m
                        break # stop iteration
                # check mapping
                if subgroup.get(d_name, None) is None:
                    log.error('subgroup_list not found in file: {}, {}'.format(
                        self.fname, list(d_mapping.keys())))
                    raise ValueError('subgroups failed')

        # format:
        # subgroup = {'stage': '2h', 'gene': 'geneA'}
        return subgroup 


    def fish_colors(self, colorPal=1, n=10):
        """Pick colors
        
        Parameters
        ----------
        colorPal : int or str
            choose the group of colors, 1 to N, or the name of fish
            candidate list: [1, 2, 3, 'Scarus_hoefleri', 'Gramma_loreto',
            'Centropyge_loricula'], default: [1]

        n : int
            number of colors return, default: 10

        # colors from fishualize R package
        # https://nschiett.github.io/fishualize/articles/overview_colors.html

        Colors from R package: 
        > fishualize::fish_pal(option = "Scarus_hoefleri")(10)
         [1] "#D2372CFF" "#E25719FF" "#F0780BFF" "#F89814FF"
         [5] "#EFB52BFF" "#C6CB43FF" "#7AD45CFF" "#00CE7DFF"
         [9] "#00BCABFF" "#0499EAFF"
        
        > fishualize::fish_pal(option = "Gramma_loreto")(10)
         [1] "#020122FF" "#1E2085FF" "#4029CBFF"
         [4] "#6628EEFF" "#901CEDFF" "#B804CAFF"
         [7] "#D61693FF" "#E6445DFF" "#EE7A30FF"
        [10] "#F0BF0BFF"

        > fishualize::fish_pal(option = "Centropyge_loricula")(10)
         [1] "#8F1D1EFF" "#B30029FF" "#DF002AFF"
         [4] "#FF7D1AFF" "#FFBD17FF" "#E7BE5AFF"
         [7] "#988591FF" "#0043A0FF" "#001A72FF"
        [10] "#000000FF"
        """
        # pre-defined colors
        c1 = ['#D2372C', '#E25719', '#F0780B', '#F89814', '#EFB52B', 
            '#C6CB43', '#7AD45C', '#00CE7D', '#00BCAB', '#0499EA']
        c2 = ['#020122', '#1E2085', '#4029CB', '#6628EE', '#901CED', 
            '#B804CA', '#D61693', '#E6445D', '#EE7A30', '#F0BF0B']
        c3 = ['#8F1D1E', '#B30029', '#DF002A', '#FF7D1A', '#FFBD17', 
            '#E7BE5A', '#988591', '#0043A0', '#001A72', '#000000']

        # RGB (10 colors, repeat twice, 20 items)
        color_d = {
            'Scarus_hoefleri': list(map(self.hex2rgb, c1*2)),
            'Gramma_loreto': list(map(self.hex2rgb, c1*2)),
            'Centropyge_loricula': list(map(self.hex2rgb, c1*2))
        }

        # get the fish_list
        fish_list = list(color_d.keys())

        # determine the fish (colorBy)
        if isinstance(colorPal, int):
            colorPal = colorPal - 1 # to 0-indexed
            if not colorPal in range(len(fish_list)):
                colorPal = 0
            fish = fish_list[colorPal]

        elif isinstance(colorPal, str):
            fish = colorPal if colorPal in color_d else 'Scarus_hoefleri'

        else:
            fish = 'Scarus_hoefleri'

        # output colors
        return color_d.get(fish, c1)[:n] 


    def hex2rgb(self, h):
        """Convert Hex to RGB code

        see: https://stackoverflow.com/a/29643643

        Parameters
        ----------
        h : str
            String hex color value

        >>> hex2rgb('#D2372C')

        """
        # extract the first 6 characters, exclude '#'
        if isinstance(h, str):
            if h.startswith('#') and len(h) >= 7:
                # h = h.lstrip('#')[:6]
                pass
            else:
                raise ValueError('hex color expected, got {}'.format(h))
        else:
            raise ValueError('hex color, str expected, got {}'.format(
                type(h).__name__))

        return ','.join(map(str, 
            [int(h[1:3], 16),
            int(h[3:5], 16),
            int(h[5:7], 16)]
            ))


    def get_color(self, colorDim=None):
        """Get color, based on sugroups (dimX, dimY, dimA)

        Parameter
        ---------
        colorDim : str
            The dimension, to assign colors, dimX, dimY, dimA, ...
    
        subgroups:
        {
            'dimX': {
                'name': 'stage',
                'label': 'stage',
                'mapping': {
                    '2h': '2hrs',
                    '4h': '4hrs'
                }
            },
            'dimY': {
                'name': 'gene',
                'label': 'Gene',
                'maping': {
                    'geneA': 'GeneA',
                    'geneB': 'GeneA'
                }
            }
        }
        """
        if colorDim is None:
            colorDim = self.colorDim

        if not colorDim in self.subgroups:
            colorDim = 'dimX'

        # self.subgroups is the subgroups structure {dimX: , dimY: , ...}
        # extract the mapping from the colorDim (dimX)
        g_mapping = self.subgroups.get(colorDim, {}).get('mapping', {})

        # assign colors to each groups in mapping
        g_list = list(g_mapping.keys()) # 2h, 4h, ...
        c_list = self.fish_colors(colorPal=self.colorPal, n=len(g_list))

        # build the color dict for the mapping
        # {2h: 'rgb', 4h: 'rgb'}
        c_dict = {k:v for k, v in zip(g_list, c_list)}

        # self.subgroup, is groups for the file
        # format: stage=2h gene=geneA

        # extract the name of colorDim
        # dimX.name == 'stage'
        colorDimName = self.subgroups.get(colorDim, {}).get('name', None)

        # extract the group for colors
        # g = '2h'
        g = self.subgroup.get(colorDimName)

        # assign color
        return c_dict.get(g, '0,0,255')


    def signal_track(self):
        """Create Track for the file: signal
        ready to insert to CompositeTrack.views
        """
        if self.ftype:
            track = trackhub.Track(
                name=self.fname+'_signal',
                short_label=self.label,
                source=self.s,
                visibility='full',
                viewLimits='0:10',
                maxHeightPixels='8:40:128',
                subgroups=self.subgroup,
                color=self.color,
                tracktype=self.ftype)
        else:
            track = None

        return track


    def region_track(self):
        """Create Track for the file: region
        ready to insert to CompositeTrack.views
        """
        if self.ftype:
            track = trackhub.Track(
                name=self.fname+'_region',
                short_label=self.label,
                source=self.s,
                visibility='dense',
                subgroups=self.subgroup,
                color=self.color,
                tracktype=self.ftype)
        else:
            track = None

        return track


    def to_track(self):
        """Choose {signal|region} Track, based on the filetype

        - signal Track, bigWig, bedGraph
        - region Track, bigBed
        """
        if self.ftype in ['bigWig', 'bedGraph', 'bigBed']:
            return self.track
        else:
            raise ValueError('bigWig, bedGraph, bigBed expected, got{}'.format(
                self.s))


class TrackHubConfig(object):
    """Parsing files, subgroups, in data_dir

    1. directory structure:
    data_dir
      |- dirA
      |   |- subgroups.yaml
      |   |- bw, bb files
      |
      |- dirB
      |   |- subgroups.yaml
      |   |- bw, bb files
    
    2. template of subgroups.yaml
    dimX:
      label: stage
      mapping:
        1h: 1h
        2h: 2h
        3h: 3h
        L1: L1
        L2: L2
        L3: L3
      name: stage
    dimY:
      label: gene
      mapping:
        shPiwi: piwi-KD
        shWhite: white-KD
      name: gene

    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_check()
        self.init_files()


    def init_args(self):
        """Required arguments for trackhub
        TrackHub args,
        HttpServer args,
        HubUrl args,
        """
        args_init = {
            'config': None, # yaml, for global args
            'data_dir': None,
            'remote_dir': None,
            'hub_name': None,
            'genome': None,
            'user': 'UCSC',
            'email': 'abc@abc.com',
            'dry_run': False,
            'short_label': None,
            'long_label': None,
            'descriptionUrl': '',
            'colorDim': 'dimX',
            'colorPal': 1,
            'label_rm_list': [],
            'http_root_dir': None,
            'http_root_alias': None,
            'http_root_url': None,
            'is_https': False,
            'mirror': None,
            'position': None,
            'validate_url': False
        }
        self = update_obj(self, args_init, force=False)

        # update from config # !!!! loading args from config.yaml
        if file_exists(self.config):
            args_config = yaml2dict(self.config)
            self = update_obj(self, args_config, force=True)

        # label
        if not isinstance(self.short_label, str):
            self.short_label = self.hub_name

        if not isinstance(self.long_label, str):
            self.long_label = self.short_label

        if not isinstance(self.email, str):
            self.email = 'abc@abc.com'


    def init_check(self):
        """
        Check files, path, ...
        """
        self.msg_args = '\n'.join([
            'Parameters',
            '{:>16s}: {}, {}'.format(
                'data_dir',
                isinstance(self.data_dir, str),
                self.data_dir),
            '{:>16s}: {}, {}'.format(
                'hub_name',
                isinstance(self.hub_name, str),
                self.hub_name),
            '{:>16s}: {}, {}'.format(
                'genome',
                isinstance(self.genome, str),
                self.genome),
            '{:>16s}: {}, {}'.format(
                'user',
                isinstance(self.user, str),
                self.user),
             '{:>16s}: {}, {}'.format(
                'email',
                isinstance(self.email, str),
                self.email)
        ])

        if not all([isinstance(i, str) for i in [self.data_dir, self.hub_name, 
            self.genome, self.user, self.email]]):
            raise ValueError(self.msg_args)


    def init_files(self):
        """Parsing files from "data_dir"

        split files into group (composite) if contains subgroups.yaml file

        directory structure:
        dirA
          |- dirB
          |   |- subgroups.yaml
          |   |- bw, bb files
          |
          |- dirC
          |   |- subgroups.yaml
          |   |- bw, bb files

        files in dirB, into compositeTrack1
        files in dirC, into compositeTrack2
        
        Support only 2-level directory structure:
        """
        # 1st-level: data_dir
        subgrp_yaml = os.path.join(self.data_dir, 'subgroups.yaml')
        g_lv1 = self.parse_group_files(data_dir=self.data_dir, 
            grp_name=self.hub_name)

        d_grp = {}
        # level-1: subgroups in data_dir
        if g_lv1:
            d_grp = {self.hub_name: g_lv1}
        else:
        # 2nd-level: data_dir/subdir
            d_lv1 = [i for i in listdir(self.data_dir, include_dir=True) 
                if os.path.isdir(i)]
            if len(d_lv1) > 0:
                g_lv2 = [self.parse_group_files(data_dir=d) for d in d_lv1]
                g_lv2 = [i for i in g_lv2 if not i is None] # remove None
                if len(g_lv2) > 0:
                    d_grp = {i['name']:i for i in g_lv2}
        # check
        if len(d_grp) == 0:
            msg = '\n'.join([
                'bigWig files or subgroups.yaml not found.',
                'Check files in: [{}]'.format(self.data_dir)
            ])
            log.error(msg)
            sys.exit(1)
            raise ValueError('bigWig files and subgroups.yaml not found: \
                check directory {}'.format(self.data_dir))

        self.subgroups_list = d_grp


    def parse_group_files(self, data_dir, grp_name=None):
        """Make sure 'subgroups.yaml' exists in the dir, "s/"

        Parameter:
        ----------
        s : str
            The directory with bw/bg files

        grp : str
            A YAML file, specify the subgroups for the track files in the dir

        return: dict{
            name: '',
            subgroups: "d/subgroups.yaml",
            bw_files: [...],
            bb_files: [...]

        else:

        return: None

        Directory structure:
        |- data_dir
        |   |- subgroups.yaml
        |   |- bw_files/
        |   |- bb_files/

        # if bw_files, subgroups.yaml exists
        """
        d_grp = None
        # init_dir
        if isinstance(data_dir, str):
            if os.path.isdir(data_dir):
                if not isinstance(grp_name, str):
                    grp_name = os.path.basename(data_dir)
                subgrp_lv1 = os.path.join(data_dir, 'subgroups.yaml')
                if file_exists(subgrp_lv1):
                    # files: level=1
                    f_lv1 = [i for i in listdir(data_dir)
                        if TrackFile(i).is_track_file()]

                    # files: level=2
                    d_lv1 = [i for i in listdir(data_dir, include_dir=True) 
                        if os.path.isdir(i)]
                    f_lv2 = [] # empty
                    if len(d_lv1) > 0:
                        [f_lv2.extend([i for i in listdir(d) 
                            if TrackFile(i).is_track_file()]) for d in d_lv1]
                    f_lv1.extend(f_lv2)

                    # output
                    d_grp = {
                        'name': grp_name,     
                        'subgroups': subgrp_lv1,
                        'bw_files': [i for i in f_lv1 if TrackFile(i).is_bw()],
                        'bb_files': [i for i in f_lv1 if TrackFile(i).is_bb()]
                    }

                    # check
                    if len(d_grp.get('bw_files')) == 0:
                        log.warning('bigWig files not found: {}'.format(
                            data_dir))
                        d_grp = None
                else:
                    log.warning('required file missing: {}'.format(subgrp_lv1))
            else:
                log.warning('not a directory: data_dir={}'.format(data_dir))
        else:
            log.warning('data_dir expect str, got {}'.format(
                type(data_dir).__name__))

        # output
        return d_grp


class TrackHub():
    """Generate multiple CompositeTracks

    Purpose:
        Organize files in directory, as CompositeTrack(s)
        return hub.txt

    Parameters
    ----------
    data_dir : str
        Path to the dir, host track files, bigWig*, bigBed 

    remote_dir : str
        Path to host the trackhub files, open accessiable, on the http server

    hub_name : str
        The name of the trackhub, [A-Za-z0-9]  

    genome : str
        The name of the UCSC genome release name, see the following link:
        https://genome.ucsc.edu/FAQ/FAQreleases.html 

    user : str
        The maintainer of this trackHub, default: [UCSC]

    email : str
        The email address of the user, default: [abc@abc.com] 

    dry_run : bool
        Try to render the tracks, create symlinks for the track files, instead
        copying files to target directory. use this option to check the 
        configuration

    short_label : str or None
        Show on the left of tracks, if [None], use hub_name instead

    long_label : str or None
        Longer labels show in the middle, if None, use short_label instead

    descriptionUrl : str
        Path to the html page, including description of the trackHub,
         default: ['']

    colorDim : str
        Assign colors to the tracks based on which dimension, dimX, dimY, ...
        default: [dimX]

    colorPal : int or str
        choose the group of colors, 1 to N, or the name of fish
        candidate list: [1, 2, 3, 'Scarus_hoefleri', 'Gramma_loreto',
        'Centropyge_loricula'], default: [1]       

    label_rm_list : list
        A list of strings, remove from the 'short_label', to make sure the 
        length of short label less than 17 characters. 

    # for HttpServer

    remote_dir : str
        Path to host the track files, under the http_root_dir, The directory 
        on the HTTP server requires open accessiable 

    http_root_dir : str
        The root_dir of HTTP server. Could be found in file
        `/etc/apache2/sites-available/000-browser.conf`, or other *.conf
        files. If VirtualHost was set, Check the active *.conf file, for the 
        <VirtualHost: ...> section.

    http_root_alias : str or None
        The Alias name of the VirtualHost directory, if [None], parse the dir
        from http_root_url

    http_root_url : str
        The URL correspond to the http_root_dir. Be aware the VirtualHost 
        configuration

    is_https : bool
        The HTTP server is using 'https'. default: [False] 

    # for HubUrl #

    hub_url : str
        URL direct to the hub.txt file, that is open accessiable

    genome : str or None
        The genome relase of the tracks, if [None], parsing this value from
        hub_url, in this way: hub_url -> hub.txt -> genomes.txt -> genome.
        default: [None]

    mirror : str
        The UCSC genome browser mirror, one of ['usa', 'euro', 'asia'] or 
        the URL of other mirror, default: [usa]

    position : str or None
        The coordinates of the position to show in default in trackhub_url,
        in this format: 'chr1:2000-4000', default: [None]

    validate_url : bool
        Validate the hub_url using hubCheck command line tool. default: [False]

    Example:
    Saving all args in config.yaml file
    >>> th = TrackHub(config='config.yaml')
    >>> th.render()
    >>> th.upload()
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        if self.demo or self.config is None:
            # trackhub_config_template()
            pass
        else:
            args_local = TrackHubConfig(**self.__dict__)
            self = update_obj(self, args_local.__dict__, force=True)
            self.init_hub() # main port, hub, trackdb, ...

            # path to hub.txt
            self.remote_hub_dir = os.path.join(self.remote_dir, self.hub_name)
            self.hub_txt = os.path.join(self.remote_hub_dir, self.hub.filename)


    def init_hub(self):
        """Initiate the trackhub
        required args:

        hub - hub.txt content
        genomes_file - genomes.txt contnet
        genome - genomes.txt content
        trackdb - main port for tracks
        """
        self.hub, self.genome_file, genome, self.trackdb = trackhub.default_hub(
            hub_name=self.hub_name,
            short_label=self.short_label,
            long_label=self.long_label,
            genome=self.genome,
            email=self.email)
            # descriptionUrl=self.descriptionUrl)


    def init_composite(self, subgroups_config):
        """Add CompositeTrack, 
        Add views to CompositeTrack

        add subgroups, if exists

        Parameter:
        ----------
        track_files : dict
            A dict saving track files and subgroups.yaml, in the following 
            format:

        track_files format:

        {
            'name': 'RNAseq', 
            'subgroups': 'subgroups.yaml',
            'bw_files': [...],
            'bb_files': [...]

        }

        subgroups.yaml format:

        dimX:
            name: stage,
            label: stage,
            mapping: 
                1h: 1h,
                2h: 2h
        dimY:
            name: gene,
            label: gene,
            mapping:
                geneA: geneA,
                geneB: geneB

        ################
        trackhub -> compositeTrack -> view -> tracks

        return [composite, signal_view, region_view]
        """
        composite_name = subgroups_config.get('name', 'composite')
        composite_subgroup = subgroups_config.get('subgroups', None)
        bw_files = subgroups_config.get('bw_files', [])
        bb_files = subgroups_config.get('bb_files', [])

        # For subgroups:
        #
        # required arguments:
        # 
        # 1. subgroups, YAML file, saving subgroups structure
        # 2. dimensions, dimX, dimY, ...
        # 3. sortOrder, gene=+ stage=+ (default: dimX, dimY, ...) 
        # 4. filterComposite, dimA 

        # load subgroups.yaml to dict
        # dimX, dimY, dimA, ...
        subgroups_dict = yaml2dict(composite_subgroup)

        if isinstance(subgroups_dict, dict):
            # prepare trackhub.SubGroupDefinition()
            # the order: dimX, dimY, dimA, ...
            subgroups_def_list = [trackhub.SubGroupDefinition(
                name=i.get('name', None),
                label=i.get('label', None),
                mapping=i.get('mapping', None)
                ) for i in list(subgroups_dict.values())]

            # prepare: dimensions
            # format: 'dimX=stage dimY=gene', or dict
            dimensions = ' '.join([
                '{}={}'.format(k, v.get('name', None)) for k, v in 
                    subgroups_dict.items()])

            # prepare: sortOder
            # format: stage=+, gene=+
            dim_list = [i.get('name', None) for i in 
                list(subgroups_dict.values()) if i.get('name', None)]
            sortOrder = ' '.join(['{}=+'.format(i) for i in dim_list])

            # prepare: filterComposite
            # format: dimA, ...
            if len(dim_list) > 2:
                filterComposite = ' '.join(list(subgroups_dict.keys())[2:])
            else:
                filterComposite = None
        else:
            subgroups_def_list = []
            dimensions = None
            sortOrder = None,
            filterComposite = None

        # init composite
        composite = trackhub.CompositeTrack(
            name=composite_name,
            short_label=self.sanitize(composite_name),
            long_label=self.sanitize(composite_name),
            dimensions=dimensions,
            sortOrder=sortOrder,
            filterComposite=filterComposite,
            tracktype='bigWig', # only for bigWig tracks
            visibility='full'
        )

        # add subgroups
        composite.add_subgroups(subgroups_def_list)

        # add signal view (bigWig, bigBed, ...)
        signal_view = trackhub.ViewTrack(
            name=composite_name+'_signal',
            view='signal',
            visibility='full',
            tracktype='bigWig',
            short_label='Signal')
        if len(bw_files) > 0:
            composite.add_view(signal_view)
        else:
            signal_view = None

        # add region view (bigBed)
        region_view = trackhub.ViewTrack(
            name=composite_name+'_region',
            view='regions',
            visibility='dense',
            tracktype='bigBed',
            short_label='Region')
        if len(bb_files) > 0:
            composite.add_view(region_view)
        else:
            region_view = None

        # output
        return (composite, subgroups_dict, signal_view, region_view)


    def sanitize(self, s, remove_list=None):
        """Sanitize a string, for shortLabel, letters, numbers allowed only [A-Za-z0-9]
        
        Convert spaces to unserscore, remove other chaaracters

        Parameters
        -----------
        s : str
            String to sanitize

        remove_list : list or None
            if list, remove the items from string; if None, skipped
        """
        # remove non-characters:
        s1 = re.sub('[^A-Za-z0-9_ ]', '', s)

        # replace spaces to underscore
        s1 = re.sub('\s+', '_', s1)

        # remove specific strings
        if isinstance(remove_list, list):
            for i in remove_list:
                if not isinstance(i, str):
                    continue
                s1 = re.sub(i, '', s1, flags=re.IGNORECASE)

        return s1


    def add_CompositeTrack(self):
        """Check: if multiple compositeTrack() required?!  self.subgroups_list

        if = 1: 
            composite_name = hub_nmae
        if > 1:
            composite_name = grp_config.get('name')
        """
        if len(self.subgroups_list) > 0:
            for subgroups_config in list(self.subgroups_list.values()):
                composite, subgroups_dict, signal_view, \
                    region_view = self.init_composite(subgroups_config)

                # add composite to trackdb
                self.trackdb.add_tracks(composite)

                # add bigWig files
                for bw in subgroups_config.get('bw_files', []):
                    tk = TrackFile(bw, 
                        subgroups=subgroups_dict,
                        colorDim=self.colorDim, # dimX, dimY, dimA, ...
                        label_rm_list=self.label_rm_list).track
                    signal_view.add_tracks(tk)

                # add bigBed files
                for bb in subgroups_config.get('bb_files', []):
                    tk = TrackFile(bb, 
                        subgroups=subgroups_dict,
                        colorDim=self.colorDim, # dimX, dimY, dimA, ...
                        label_rm_list=self.label_rm_list).track
                    region_view.add_tracks(tk)
        else:
            raise ValueError('Track() failed')


    def _tmp(self, is_dir=False):
        """Create a tmp file/directory
        """
        if is_dir:
            tmp = tempfile.TemporaryDirectory(prefix='trackhub')
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='trackhub', delete=False)
        return tmp.name


    def render(self):
        """To get ready for uploading, render the hub files, symlink all the
        track source files into a local directory

        The directory is ready for transferring to a host.
        """
        # construct tracks()
        self.add_CompositeTrack()

        # create symlinks
        tmp_dir = self._tmp(is_dir=True)
        log.info('Render track hub to: {}'.format(tmp_dir))
        trackhub.upload.stage_hub(self.hub, staging=tmp_dir)

        # hub.txt file
        return os.path.join(tmp_dir, self.hub.filename)


    def upload(self):
        """Upload the hub files, source files to server
        """
        if isinstance(self.remote_dir, str):
            # dir, to save track files, hub files
            # self.remote_hub_dir = os.path.join(self.remote_dir, self.hub_name)
            check_path(self.remote_hub_dir)

            # construct tracks()
            tmp_dir = self.render()
            path_remove(tmp_dir, ask=False)

            # upload files to remote/server
            trackhub.upload.upload_hub(self.hub, host='localhost',
                remote_dir=self.remote_hub_dir)

            # return hub filename
            return self.hub_txt

        else:
            log.error('remote_dir failed, str expected, got {}'.format(
                type(self.remote_dir).__name__))


    def get_hub_url(self):
        """Convert hub.txt to hub_url

        required parameters
        - http_root_dir
        - http_root_alias
        - http_root_url

        optional parameters
        - genome 
        - mirror
        - validate_url
        - position
        """
        # upload files
        if not file_exists(self.hub_txt):
            self.upload()

        # convert hub_txt to hub_url
        self.hub_url = HttpServer(s=self.hub_txt, **self.__dict__).to_url()

        # save hub_url to hub.yaml
        hub_yaml = self.remote_hub_dir + '/hub.yaml'
        hub_dict = {'hub_url': self.hub_url}
        dict2yaml(hub_dict, hub_yaml)

        return self.hub_url


    def get_trackhub_url(self):
        """Convert hub.txt to trackhub_url

        required parameters
        - http_root_dir
        - http_root_alias
        - http_root_url

        optional parameters
        - genome 
        - mirror
        - validate_url
        - position
        """
        # upload files
        if not file_exists(self.hub_txt):
            self.upload()
        
        # hub_url
        self.get_hub_url() # self.

        # trackhub_url
        self.trackhub_url = HubUrl(**self.__dict__).trackhub_url()
        
        # save hub_url to hub.yaml
        hub_yaml = self.remote_hub_dir + '/hub.yaml'
        hub_dict = {'hub_url': self.hub_url, 'trackhub_url': self.trackhub_url}
        dict2yaml(hub_dict, hub_yaml)

        # target files
        print('>> find hub_txt: {}'.format(hub_yaml))

        return self.trackhub_url


    def run(self):
        """Run for all

        1. --dry-run, self.render() 
        2. self.trackhub_url()

        Parameters
        ----------
        dry_run : bool
            Construct trackhub files, but not copy files

        """
        if self.demo or self.config is None:
            trackhub_config_template()
        elif self.dry_run:
            self.render()
        else:
            self.get_trackhub_url()


def trackhub_config_template():
    """Generate config.yaml template
    set default values

    save args to file: config_template.yaml
    save subgroups to file: subgroups_ytemplate.yaml
    """
    default_args = {
        'hub_name': 'tracks',
        'data_dir': '', 
        'remote_dir': '',
        'label_rm_list': ['RNAseq_', 'ATACseq_', 'CnR_'],
        'colorDim': 'dimX',
        'http_root_dir': '/data/public/upload',
        'http_root_alias': '/upload',
        'http_root_url': None,
        'genome': 'dm6',
        'user': 'user',
        'email': 'abc@abc.com',
        'dry_run': False,
        'short_label': None,
        'long_label': None,
        'descriptionUrl': '',
        'colorPal': 1,
        'is_https': False,
        'mirror': 'usa',
        'position': None,
        'validate_url': False
    }

    default_subgroups = {
        'dimX': {
            'name': 'stage',
            'label': 'Stage',
            'mapping': {
                '1h': '1h',
                '2h': '2h'
            }
        },
        'dimY': {
            'name': 'gene',
            'label': 'gene',
            'mapping': {
                'geneA': 'geneA',
                'geneB': 'geneB'
            }
        }
    }

    # default files
    config_f = 'config_template.yaml'
    subgroups_f = 'subgroups_template.yaml'
    dict2yaml(default_args, config_f)
    dict2yaml(default_subgroups, subgroups_f)

    # show message
    cf = '{:>16s} : {}'.format('config.yaml', config_f),
    gf = '{:>16s} : {}'.format('subgroups.yaml', subgroups_f),

    msg = '\n'.join([
        '{}'.format('#'*80),
        '# {:<77s}#'.format('Mini-tutorial [run_trackhub]'),
        '# {:<77s}#'.format('1. Generating template files:'),
        '# {:<77s}#'.format('$ python run_trackhub --demo'),
        '# {:<77s}#'.format(str(cf)),
        '# {:<77s}#'.format(str(gf)),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('2. Move the YAML files to {data_dir}'),
        '# {:<77s}#'.format('$ mv subgroups.yaml {data_dir}/subgroups.yaml'),
        '# {:<77s}#'.format('$ mv config.yaml {data_dir}/config.yaml'),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('Update the YAML files, according to your data'),
        '# {:<77s}#'.format('Required fields - config.yaml'),
        '# {:<77s}#'.format('  - data_dir        # absolute path'),
        '# {:<77s}#'.format('  - genome          # dm6'),
        '# {:<77s}#'.format('  - label_rm_list   # string, removed from label'),
        '# {:<77s}#'.format('  - hub_name        # RNAseq_piwi'),
        '# {:<77s}#'.format('  - remote_dir      # see http_root_dir'),
        '# {:<77s}#'.format('  - position        # chr2L:1-1000'),
        '# {:<77s}#'.format('  - http_root_alias # /upload, see: /var/apache2/sites-available/'),
        '# {:<77s}#'.format('  - http_root_dir   # /data/public/upload, as above'),
        '# {:<77s}#'.format('  - http_root_url   # null, parse IP of HTTP server'),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('Required fields - subgroups.yaml'),
        '# {:<77s}#'.format('  - mapping'),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('3. Generating trackhub files'),
        '# {:<77s}#'.format('$ hiseq run_trackhub -c {data_dir}/config.yaml'),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('4. Find the hub.txt file'),
        '# {:<77s}#'.format('{remote_dir}/hub_name/{hub_name}_hub.txt'),
        '{}'.format('#'*80)
    ])

    print(msg)


def main():
    trackhub_config_template()


if __name__ == '__main__':
    main()
