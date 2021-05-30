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
import trackhub
from shutil import which
from urllib.request import urlopen
from bs4 import BeautifulSoup
from hiseq.utils.utils import log, update_obj, Config, get_date, run_shell_cmd



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
                which('hubCheck'),
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
