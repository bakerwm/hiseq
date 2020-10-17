"""
Generate a trackHub directory for bigWig, bigBed files

Deatails about Track Hubs: http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html#Intro 

The trackhub python package [trackhub](https://daler.github.io/trackhub/quickstart.html) is designed to generate the files

1. rename samples
2. pick colors
3. 

bigWig + bigBed (bigNarrowPeak)

bigWig


Features:
1. Search files (bigWig, bigBed)
2. Generate trackhub for local files  
3. Make local files available online (server)  


4. Utils: 
  - narrowPeak -> bb
  - change values    
  - choose colors 
  - rename files

"""

import os
import sys
import re
import pathlib
import shutil
import fnmatch
import argparse
import json
import glob
import shlex
import logging
import subprocess
import collections
import trackhub
from hiseq.utils.helper import Json
from urllib.request import urlopen
from bs4 import BeautifulSoup

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)

log.setLevel('INFO')


def print_dict(d):
    d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog='get_trackhub',
        description='Generate trackhub for bigWig and bigBed files',
        epilog='Example: \n\
               python get_trackhub.py -i bigWig -n ChIPseq -g dm6')
    parser.add_argument('-i', '--data-dir', required=True, dest='data_dir',
        help='The directory of bigWig and bigBed files')
    parser.add_argument('-o', '--remote-dir', required=True, dest='remote_dir',
        help='The directory to save the track files')
    parser.add_argument('-r', '--recursive', action='store_true',
        help='search files in data_dir Recursively')
    parser.add_argument('-n', '--hub-name', metavar='hub_name', required=False,
        default=None, help='hub name')
    parser.add_argument('-g', '--genome', metavar='GENOME', required=False,
        default='dm6', help='genome for the trackhub, UCSC genome build, \
        [hg19, hg38, mm9, mm10, dm3, dm6]')
    parser.add_argument('-l', '--short-label', default=None, dest='short_label',
        help='short label for the hub, default: [--hub-name]')
    parser.add_argument('-L', '--long-label', default=None, dest='long_label',
        help='long label for the hub, default: [--hub]')
    parser.add_argument('-u', '--user', default='UCSC',
        help='Who maintain the trackhub')
    parser.add_argument('-e', '--email', default='abc@abc.com',
        help='email of the maintainer')
    parser.add_argument('-m', '--mirror', default='usa',
        help='The mirror of UCSC, [usa|asia|euro], or custome mirror, input \
        url of your UCSC_mirror: default: [usa]')
    parser.add_argument('-s', '--subgroups-config', dest='subgroups_config',
        default=None,
        help='The config for subgroups, default: [None]')
    parser.add_argument('-t', '--http-config', dest='http_config',
        default=None,
        help='The config for http, open access, including host, root_dir, \
        default [None]')
    parser.add_argument('--http-host', dest='http_host', default=None,
        help='The http server host url, example: http://abc.com/upload')
    parser.add_argument('--http-root-dir', dest='http_root_dir', default=None,
        help='The http server, root_dir, example: /data/upload')
    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
        help='Do not copy the files')
    args=parser.parse_args()
    return args


def update_obj(obj, d, force=True, remove=False):
    """
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    # fresh start
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)

    return obj


## utils
def listdir(path, full_name=True, recursive=False, include_dir=False):
    """
    List all the files within the path
    """
    out = []
    for root, dirs, files in os.walk(path):
        if full_name:
            dirs = [os.path.join(root, d) for d in dirs]
            files = [os.path.join(root, f) for f in files]
        out += files

        if include_dir:
            out += dirs

        if recursive is False:
            break

    return sorted(out)


def listfile(path='.', pattern='*', full_name=True, recursive=False):
    """
    Search files by the pattern, within directory
    fnmatch.fnmatch()

    pattern:

    *       matches everything
    ?       matches any single character
    [seq]   matches any character in seq
    [!seq]  matches any char not in seq

    An initial period in FILENAME is not special.
    Both FILENAME and PATTERN are first case-normalized
    if the operating system requires it.
    If you don't want this, use fnmatchcase(FILENAME, PATTERN).

    example:
    listfile('./', '*.fq')
    """
    fn_list = listdir(path, full_name, recursive, include_dir=False)
    fn_list = [f for f in fn_list if fnmatch.fnmatch(os.path.basename(f), pattern)]
    return sorted(fn_list)


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


## functions for Track Hub ##
class HubUrl(object):
    """
    For Huburl

    1. hub.txt -> hub_url 
    2. hub.txt -> hub_url -> trackhub_url

    Option2:
    read: hub.txt, genomes.txt from string (trackhub.default_hub() output)
    
    hub, genome_file, _, trackdb = trackhub.default_hub()


    Example:
    HubUrl(
        hub='/data/tracks/hub.txt', 
        http_host='http://abc.com/tracks',
        http_root_dir='/data/tracks')

    HubUrl(
        hub='http://abc.com/tracks/hub.txt', 
        http_host='http://abc.com/tracks',
        http_root_dir='/data/tracks')
    
    HubUrl(
        hub='/data/tracks/hub.txt',
        http_config='/data/http_config.json'
    )

    ## work with trackhub package

    hub, genome_file, _, trackdb = trackhub.default_hub(
        hub_name='hub_name',
        short_label='short_label',
        long_label='long_label',
        genome='genome',
        email='email',
        descriptionUrl='descriptionUrl')

    HubUrl(
        hub=str(hub),
        genome_file=str(genome_file),
        remote_dir='/data/tracks',
        http_config='/data/http_config.json'
    )

    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.config()


    def config(self):
        """        
        trackhub.default_hub()
        """
        args_init = {
            'hub': None,
            'genome_file': None,
            'remote_dir': None,
            'http_config': None,
            'http_host': None,
            'http_root_dir': None}
        self = update_obj(self, args_init, force=False)

        # url: str
        if not isinstance(self.hub, str):
            raise ValueError('url, expect str, got {}'.format(
                type(self.hub).__name__))

        # update hub_url
        self.init_hub() #

        # update genome
        args_genome = self.read_genome_txt()
        self.genome = args_genome.get('genome', None)

        # http server
        self.init_http()


    def init_hub(self):
        """
        Construct url from: args, {hub,genome_file,remote_dir}

        generate url: http://hub.txt
        
        option1:
        hub: 
            - path to hub.txt
            - content of hub.txt
            - path to http://hub.txt

        genome_file:
            - {hub.txt exists} hub -> genome.txt
            - {hub.txt not exists} hub string -> genome.txt

        remote_dir:
            - required if hub is content of hub.txt
        """
        if self.is_hub_txt(self.hub):
            if not isinstance(self.remote_dir, str):
                raise ValueError('remote_dir required, if hub is text')
            args_hub = self.read_hub_txt()
            genomesFile = args_hub.get('genomesFile', None)
            hub_txt_name = re.sub('genomes.txt', 'hub.txt', genomesFile)
            hub_txt = os.path.join(self.remote_dir, hub_txt_name)
            hub_url = re.sub(self.http_root_dir, self.http_host, hub_txt)
        elif self.is_hub_local(self.hub):
            hub_url = re.sub(self.http_root_dir, self.http_host, self.hub)
        elif self.is_hub_url(self.hub):
            hub_url = self.hub
        else:
            raise ValueError('hub, illegal, {}'.format(self.hub))

        self.hub_url = hub_url


    def init_http(self):
        """
        Deposite the track files in the same server, running this program
        require http config for the server
        http_host:     the url of the server, open accessiable
        http_root_dir: the root dir of the http server
        """
        if not self.is_hub_url():
            # read http_config, update http_host, http_root_dir
            if isinstance(self.http_config, str):
                args_http = Json(self.http_config).reader()
                self.http_host = args_http.get('host', None)
                self.http_root_dir = args_http.get('root_dir', None)

            # check values
            if not isinstance(self.http_host, str):
                raise ValueError('http_host expect str, got {}'.format(
                    type(self.http_host).__name__))
                
            if not isinstance(self.http_root_dir, str):
                raise ValueError('http_root_dir expect str, got {}'.format(
                    type(self.http_root_dir).__name__))


    def is_hub_txt(self, x=None):
        """
        Check if hub is text of (hub.txt)
        template:

        hub hub
        shortLabel test
        longLabel test
        genomesFile test.genomes.txt
        email mail@abc.com
        descriptionUrl
        """
        if x is None:
            x = self.hub

        p = re.compile('^hub\s.*\ngenomesFile\s.*', re.DOTALL)
        return True if p.search(x) else False


    def is_hub_local(self, x=None):
        """
        Check if the 'url' is the local directory:
        """
        if x is None:
            x = self.hub

        p = re.compile('^/\w+', re.DOTALL)
        return True if p.search(x) else False


    def is_hub_url(self, x=None):
        """
        Check if the 'url' is url: http://ip, or ftp://ip
        """
        if x is None:
            x = self.hub

        # match
        p = re.compile('(^http|^https|^ftp)://')
        return True if p.search(x) else False

        # tag = False
        # if isinstance(x, str):
        #     if re.split(':', x)[0] in ['ftp', 'http', 'https']:
        #         tag = True

        # return tag


    def is_genome_file_txt(self, x=None):
        """
        Check if genome_file is text of (genomes.txt)
        template:

        hub hub
        shortLabel test
        longLabel test
        genomesFile test.genomes.txt
        email mail@abc.com
        descriptionUrl
        """
        if x is None:
            x = self.genome_file

        if x is None:
            return False

        p = re.compile('^genome\s\w+.*trackDb\s\w+', re.DOTALL)
        return True if p.search(x) else False


    def read_local(self, x):
        """
        Read text file: x
    
        format:
        key value
        """
        d = {}
        try:
            with open(x, 'rt') as r:
                s = re.split('\s', r.read())
            sx = iter(s)
            d = dict(zip(sx, sx))
        except:
            log.error('reading file failed, {}'.format(x))

        # update
        return d


    def read_url(self, url):
        """
        Read content from url
        """
        d = {}
        try:
            html = urlopen(url).read()
            text = BeautifulSoup(html, features='html.parser')
            x = re.split('\s', text.text)
            xit = iter(x)
            d = dict(zip(xit, xit))
            # # add attr 
            # for k, v in d.items():            
            #     if not hasattr(self, k):
            #         setattr(self, k, v)
        except:
            log.error('reading url failed: {}'.format(url))

        return d


    def read_hub_txt(self):
        """
        hub.txt
        hub, shortLabel, longLabel, genomesFile, descriptionUrl

        genomes.txt
        genome, trackDb
        """
        if self.is_hub_txt():
            x = re.split('\s', self.hub)
            xit = iter(x)
            args = dict(zip(xit, xit))
        elif self.is_hub_local():
            args = self.read_local(self.hub)
        elif self.is_hub_url():
            args = self.read_url(self.hub)
        else:
            args = {}

        return args


    def read_genome_txt(self):
        """
        Get the genome from hubUrl, genomes.txt
        """
        if self.is_hub_txt():
            # require genome_file
            if self.is_genome_file_txt():
                x = re.split('\s', self.genome_file)
                xit = iter(x)
                args = dict(zip(xit, xit))
            else:
                raise ValueError('genome_file, required, for example: \
                    [{}]'.format(r"'genome dm6\ntrackDb dm6/trackDb.txt\n\n'"))
        elif self.is_hub_local():
            args_hub = self.read_hub_txt()
            genomesFile = args_hub.get('genomesFile')
            genome_file = os.path.join(os.path.dirname(self.hub), genomesFile)
            args = self.read_local(genome_file)
        elif self.is_hub_url():
            args_hub = self.read_hub_txt()
            genomesFile = args_hub.get('genomesFile')
            genome_file = os.path.join(os.path.dirname(self.hub), genomesFile)
            args = self.read_url(genome_file)
        else:
            raise ValueError('genome_file, illegal, {}'.format(genome_file))

        return args


    def format_position(self, position=None):
        """
        Convert position to url format
        Input: chr1:1-1000
        Output: &position=chr1%3A1-1000

        see: UCSC Documentation
        https://genome.ucsc.edu/goldenPath/help/customTrack.html#optParams
        """
        if isinstance(position, str):
            position = re.sub('[^\w:-]', '', position) # sanitize chr
            p = re.compile('^(\w+):([\d,]+)-([\d,]+)$')
            
            if p.search(position):
                chr, start, end = re.split('[:-]', position)
                suffix = '&position={}%3A{}-{}'.format(chr, start, end)
            else:
                log.warning('illegal position, chr1:1-100 expected, \
                    got {}'.format(position))
                suffix = ''
        else:
            log.warning('illegal position, str expected, got {}'.format(
                type(position).__name__))
            suffix = ''

        return suffix


    def to_trackhub_url(self, genome_browser_mirror='usa', position=None):
        """
        Convert hub.txt url to trackhub Url for sharing
        output: 
        trackhub_url: http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://ip/hub.txt 

        region: usa, asia, euro

        usa:  http://genome.ucsc.edu
        asia: http://genome-asia.ucsc.edu
        euro: http://genome-euro.ucsc.edu
        custome: ''
        """
        mirrors = {
            'usa': 'http://genome.ucsc.edu',
            'asia': 'http://genome-asia.ucsc.edu',
            'euro': 'http://genome-euro.ucsc.edu'}

        # determine the mirror
        if self.is_hub_url(genome_browser_mirror):
            mirror = genome_browser_mirror
        elif genome_browser_mirror in mirrors:
            mirror = mirrors.get(genome_browser_mirror, None)
        else:
            raise ValueError('genome_browser_mirror, illegal value, \
                http://, or {usa|asia|euro} ecpect, got {}'.format(
                    type(genome_browser_mirror).__name__))

        # add genome
        trackhub_url = '{}/cgi-bin/hgTracks?db={}'.format(mirror, self.genome)
        
        # add position
        trackhub_url += self.format_position(position)

        # add hub_txt
        trackhub_url += '&hubUrl={}'.format(self.hub_url)

        return trackhub_url


class Subgroup(object):
    """
    Create subgroups, based on filename
    
    pre-defined groups: ATACseq, RNAseq, ChIPseq, HiSeq, ...
    config: json file, pre-defined subgroups, 

    Sub groups for the project:
    ATACseq: gene, stage, kind
    RNAseq: gene, stage, strand
    ChIPseq: rep, strand, kind

    check the dimension="dimX=strand dimY=rep dimA=kind"

    or load from config.json
    dimX: {keywords}  
    dimY: {keywords}
    dimA: {keywords}
    """
    def __init__(self, filelist=None, config=None, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.filelist = filelist
        self.config = config
        self.init_args()


    def init_args(self):
        """
        Read filenames from filelist
        Read config from config

        config:
        dimX: label: content: 
        dimY:
        dimA
        """
        if isinstance(self.filelist, str):
            stype = TrackFile(self.filelist)
        elif isinstance(self.filelist, list):
            # check files: ATACseq|ChIPseq|RNAseq ...
            stype = [TrackFile(i).seqtype for i in self.filelist]
            stype = list(set(stype)) # unique
            stype = stype.pop()
        else:
            stype = 'default'
        
        # load config
        if self.config is None:
            groups = self.init_group() # default
        else:
            groups = Json(self.config).reader()

        # choose group        
        if not stype in groups:
            stype = 'default'

        self.group = groups.get(stype, None)

        # dimensions
        self.dimensions = 'dimX={} dimY={} dimA={}'.format(
            self.group.get('dimX').get('name'),
            self.group.get('dimY').get('name'),
            self.group.get('dimA').get('name'))


    def init_group(self):
        """
        Predefined groups for samples

        required: name, label, mapping
        
        ChIPseq: 
        dimX: antibody
        dimY: gene/cell_line
        dimA: peak/signal

        """
        df = {
            "default": {
                "dimX": {
                    "name": "kind",
                    "label": "kind",
                    "mapping": {
                        "signal": "signal",
                        "peak": "peak"
                    }
                },
                "dimY": {
                    "name": "kind",
                    "label": "kind",
                    "mapping": {
                        "signal": "signal",
                        "peak": "peak"
                    }
                },
                "dimA": {
                    "name": "kind",
                    "label": "kind",
                    "mapping": {
                        "signal": "signal",
                        "peak": "peak"
                    }
                }
            }
        }

        return df


    def to_trackhub(self):
        """
        Convert subgroup to trackhub.SubGroupDefinition()
        """
        return [trackhub.SubGroupDefinition(
            name = v.get('name', None),
            label = v.get('label', None),
            mapping = v.get('mapping', None)) for k, v in self.group.items()]


class TrackColor(object):
    def __init__(self):
        pass

    def get_strand_color(self):
        """Pick color based on the filename
        f:      forward, '#2E3440'
        r:      reverse, '#6DBFA7'
        other:  others, '#lalala'

        example: 
        rnaseq.fwd.bigWig
        
        to-do:
        assign multiple colors for tracks/files    
        """
        # Due to how code is extracted from the docs and run during tests, we
        # need to import again inside a function. You don't normally need this.
        colors = {
            'f': '#2E3440',
            'r': '#6DBFA7',
            'other': '#1a1a1a',
        }
        return trackhub.helpers.hex2rgb(colors[self.strand])


class TrackFile(object):
    """
    Check files: bigWig, bigBed, ...
    filename
    strand
    replicate
    antibody
    cell-line

    filename format:
    {seqtype}_{source}_{gene}_{treat}_{rep}_{strand:f|r|other}
    """
    def __init__(self, file):
        self.file = file
        self.init_args()
        self.track = self.get_track()


    def init_args(self):
        """
        Default values
        """
        self.filename = os.path.basename(self.file)
        self.fname = trackhub.helpers.sanitize(self.filename) #
        self.strand = self.get_strand()
        self.filetype = self.get_filetype()
        self.seqtype = self.get_seqtype()
        self.groups = self.get_groups()
        self.color = '0,0,254' # red


    def get_strand(self):
        """
        Determine the strand of file based on filename
        watson, forward, fwd, plus
        crick, reverse, rev, minus
        """
        fname = os.path.splitext(self.file)[0].lower()

        x = re.split('[._]', fname)[-1] # in the tails

        if x in ['+', 'forward', 'fwd', 'plus', 'watson']:
            strand = 'f'
        elif x in ['-', 'reverse', 'rev', 'minus', 'crick']:
            strand = 'r'
        else:
            strand = 'other'

        return strand


    def get_filetype(self):
        """
        Determine the type of file
        bigWig: bigWig, bigwig, bw
        bigBed: bigBed, bigbed, bb
        bedGraph: bedGraph, bedgraph, bg
        """
        ftype = os.path.splitext(self.file)[1].lower().lstrip('.')

        if ftype in ['bigwig', 'bw']:
            ftype = 'bigWig'
        elif ftype in ['bigbed', 'bb']:
            ftype = 'bigBed'
        elif ftype in ['bedgraph', 'bg']:
            ftype = 'bedGraph'
        else:
            ftype = None
            raise ValueError('unknown filetype: {}'.format(self.filename))

        return ftype


    def get_seqtype(self):
        """
        Determine the sequencing type of the file
        ATACseq, RNAseq, ChIPseq
        default: ChIPseq
        """
        stype = re.split('[._]', os.path.basename(self.file))[0] # head
        stype = stype.lower()

        if stype in ['atac', 'atacseq', 'cutandtag']:
            stype = 'ATACseq'
        elif stype in ['rnaseq', 'mrnaseq']:
            stype = 'RNAseq'
        elif stype in ['chip', 'chipseq', 'cutandrun']:
            stype = 'ChIPseq'
        elif stype in ['smrna', 'smallrna', 'smrnaseq', 'smallrnaseq']:
            stype = 'smRNAseq',
        else:
            stype = 'HiSeq'

        return stype


    def signal_track(self):
        """
        Create signal track for the file: bigWig, bedGraph
        """
        track = trackhub.Track(
            name=self.fname,
            source=self.file,
            visibility='full',
            viewLimits='0:10',
            maxHeightPixels='8:40:128',
            subgroups=self.groups,
            color=self.color,
            tracktype=self.filetype)

        return track


    def region_track(self):
        """
        Create region track for the file: bigBed, bed, ...
        """
        track = trackhub.Track(
            name=self.fname,
            source=self.file,
            visibility='dense',
            subgroups=self.groups,
            color=self.color,
            tracktype=self.filetype)

        return track


    def get_track(self):
        """
        Create Track for file
        bigWig, bigBed
        """
        if self.filetype in ['bigWig', 'bedGraph']:
            track = self.signal_track()
        elif self.filetype in ['bigBed']:
            track = self.region_track()
        else:
            raise ValueError('unknown filetype: {}'.format(self.filename))

        return track


    def get_groups(self):
        """
        Extract group info from filename:
        groups: gene, stage, filetype
        
        !!! Warning !!!
        This function only for custome project:
        filename:
        ATACseq_DaGal4shGene_stage_rep1
        
        gene: shGene, shGene1_shGene2
        stage: 1h, 2h, 
        kind: peak/signal
        """
        fname = self.filename.lower()

        ## shGene_
        p1 = re.compile('gal4x.?sh([A-Za-z0-9]+(\_[^0-9][A-Za-z0-9]+)?)', 
            flags=re.IGNORECASE)
        p1x = p1.search(fname)
        if p1x is None:
            gene = p1x[1]
        else:
            gene = 'gene'

        ## stage (3h, ... L1, )
        p2 = re.compile('[-_]?(\d+h|L[123]|\d+h?_\d+h)_?', flags=re.IGNORECASE)
        p2x = p2.search(fname)
        if p2x:
            stage = p2x[1]
        else:
            stage = 'stage'

        # output
        groups = {
            'gene': gene, 
            'stage': stage, 
            'kind': self.filetype}

        ## for ChIPseq
        if self.filetype in ['ChIPseq', 'CUTNRUN']:
            ## antibody (ChIPseq, CUT&RUN)
            rbp_list = ['H3K4m1', 'H3K4me2', 'H3K4me3', 
                        'H3K9me1' 'H3K9me2', 'H3K9me3',
                        'H3K27me3', 'H3K27ac', 'CTCF',
                        'GFP', 'eGFP', 'IgG']
            # p3 = re.compile('(H3K(4|9|27|36)me\d|CTCF|)')
            p3 = [i for i in rbp_list if i in fname]
            if len(p3) == 1:
                rbp = p3.pop()
            elif len(p3) > 1:
                log.warning('multiple RBP detected in file: {}'.format(fname))
                rbp = p3.pop()
            else:
                rbp = 'Gene'
            
            # output
            groups = {
            'antibody': rbp,
            'stage': stage, 
            'kind': self.filetype}

        return groups


class TrackHubConfig(object):
    """
    Prepare files, args:

    data_dir: path-to-bw-files
    remote_dir: path-to-save-bw-files
    http_root_dir:

    working mode:
    local_files -> local_path -> remote_path -> hub.txt
    -> hub_url -> trackhub_url
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_http()


    def init_args(self):
        """
        required arguments for trackhub
        """
        args_init = {
            'data_dir': None,
            'hub_name': None,
            'remote_dir': None,
            'short_label': None,
            'long_label': None,
            'genome': None,
            'user': None,
            'email': None,
            'descriptionUrl': None,
            'recursive': False,
            'subgroup_config': None,
            'http_config': None,
            'http_host': None,
            'http_root_dir': None,
            'mirror': 'usa',
            'dry_run': False}
        self = update_obj(self, args_init, force=False)

        # check data_dir, remote_dir
        if not isinstance(self.data_dir, str):
            raise ValueError('data_dir, expect str, got {}'.format(
                type(self.data_dir).__name__))

        # hub name: required
        if not isinstance(self.hub_name, str):
            raise ValueError('hub_name, expect str, got {}'.format(
                type(self.hub_name).__name__))

        # remote dir
        if not isinstance(self.remote_dir, str):
            raise ValueError('remove_dir, expect str, got {}'.format(
                type(self.remote_dir).__name__))
        # update remote_dir #!!!!
        self.remote_dir = os.path.join(self.remote_dir, self.hub_name)

        # label
        if self.short_label is None:
            self.short_label = self.hub_name

        if self.long_label is None:
            self.long_label = self.short_label

        # genome
        if not isinstance(self.genome, str):
            raise ValueError('genome, expect str, got{}'.format(
                type(self.genome).__name__))

        # user,email
        if not isinstance(self.user, str):
            self.user = 'UCSC'

        if not isinstance(self.email, str):
            self.email = 'mail@abc.com'

        # descriptionUrl
        if not isinstance(self.descriptionUrl, str):
            self.descriptionUrl = '' # empty

        # available files in data_dir
        self.file_list = self.get_files()
        if len(self.file_list) < 1:
            raise ValueError('no files detected: {}'.format(self.data_dir))
        # groups
        self.bw_files = self.get_bw_files()
        self.bb_files = self.get_bb_files()
        self.peak_files = self.get_peak_files()
        self.bam_files = self.get_bam_files()

        # bw files(>0)
        if len(self.bw_files) < 1:
            raise ValueError('bigWig not detected: {}'.format(self.data_dir))


    def init_http(self):
        """
        Deposite the track files in the same server, running this program
        require http config for the server
        http_host:     the url of the server, open accessiable
        http_root_dir: the root dir of the http server
        """
        # read http_config, update http_host, http_root_dir
        if isinstance(self.http_config, str):
            args_http = Json(self.http_config).reader()
            self.http_host = args_http.get('host', None)
            self.http_root_dir = args_http.get('root_dir', None)

        # check values
        if not isinstance(self.http_host, str):
            raise ValueError('http_host expect str, got {}'.format(
                type(self.http_host).__name__))
            
        if not isinstance(self.http_root_dir, str):
            raise ValueError('http_root_dir expect str, got {}'.format(
                type(self.http_root_dir).__name__))


    def get_files(self):
        """
        Get all bigWig/bigBed files
        
        recursive:{True|False}

        dir1
          |- files
          |- dir2
          |    |- files
        """
        if self.recursive:
            flist = listdir(self.data_dir, recursive=True, include_dir=False)
        else:
            # level1 + level2 
            xlist = listdir(self.data_dir, recursive=False, include_dir=True)
            flist1 = [i for i in xlist if os.path.isfile(i)]
            flist2 = [listdir(i, recursive=False, include_dir=False) 
                for i in xlist if os.path.isdir(i)]

            # merge two list
            flist2.append(flist1)
            flist2 = [i for i in flist2 if len(i) > 0] # remove empty list

            # flat list
            flist = [i for s in flist2 for i in s]

        return flist


    def get_bw_files(self):
        """
        List the track files, bigWig, bigbed, bedgraph
        bigWig, bw, bigwig, bigWig
        bigBed, bb, bigbed, bigBed
        # bedGraph, bg, bedgraph, bedGraph
        """
        return [i for i in self.file_list if 
            os.path.splitext(i)[1] in ['.bigWig', '.bigwig', '.bw']]


    def get_bb_files(self):
        """
        List the track files, bigWig, bigbed, bedgraph
        bigWig, bw, bigwig, bigWig
        bigBed, bb, bigbed, bigBed
        # bedGraph, bg, bedgraph, bedGraph
        """
        return [i for i in self.file_list if 
            os.path.splitext(i)[1] in ['.bigBed', '.bigbed', '.bb']]


    def get_bam_files(self):
        """
        List the track files, bigWig, bigbed, bedgraph
        bigWig, bw, bigwig, bigWig
        bigBed, bb, bigbed, bigBed
        # bedGraph, bg, bedgraph, bedGraph
        """
        return [i for i in self.file_list if 
            os.path.splitext(i)[1] in ['.bam']]


    def get_peak_files(self):
        """
        List the track files, bigWig, bigbed, bedgraph
        bigWig, bw, bigwig, bigWig
        bigBed, bb, bigbed, bigBed
        # bedGraph, bg, bedgraph, bedGraph
        peak: .narrowPeak, .broadPeak
        """
        return [i for i in self.file_list if 
            os.path.splitext(i)[1] in ['.narrowPeak']]


class TrackHub(object):
    """
    Generate tracks for specific files
    
    Composite track (currently)
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = TrackHubConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        self.init_hub()
        self.init_composite()
        # save config

        hub = HubUrl(
            hub=str(self.hub),
            genome_file=str(self.genome_file),
            remote_dir=self.remote_dir,
            http_host=self.http_host,
            http_root_dir=self.http_root_dir)

        self.hub_url = hub.hub_url
        self.trackhub_url = hub.to_trackhub_url(self.mirror)


    def init_hub(self):
        self.hub, self.genome_file, _, self.trackdb = trackhub.default_hub(
            hub_name=self.hub_name,
            short_label=self.short_label,
            long_label=self.long_label,
            genome=self.genome,
            email=self.email,
            descriptionUrl=self.descriptionUrl)


    def init_composite(self):
        """
        Initiate the composite track
        default
        dimensions: dimX={} dimY={} dimA={}
        """
        # determine subgroup by seqtype
        subgroup = Subgroup(filelist=self.bw_files,
            config=self.subgroup_config)

        ## basic composite track 
        # Create the composite track
        composite = trackhub.CompositeTrack(
            name='composite',
            short_label=self.short_label,
            long_label=self.long_label,
            dimensions=subgroup.dimensions,
            filterComposite='dimA',
            tracktype='bigWig',
            visibility='full'
        )

        # Add subgroups to the composite track
        composite.add_subgroups(subgroup.to_trackhub())

        # Add the composite track to the trackDb
        self.trackdb.add_tracks(composite)

        self.composite = composite # CompositeTrack


    def signal_track_bw(self):
        """
        Default Track for bigWig

        # CompositeTracks compose different ViewTracks.
        # one for signal in bigWig, another one for bigBed regions.
        """
        return trackhub.ViewTrack(
            name='signalviewtrack',
            short_label='Signal',
            view='signal',
            visibility='full',
            tracktype='bigWig')


    def signal_track_bg(self):
        """
        Default Track for bedGraph

        # CompositeTracks compose different ViewTracks.
        # one for signal in bigWig, another one for bigBed regions.
        """
        return trackhub.ViewTrack(
            name='singalviewtrack2',
            short_label='Signal2',
            view='signal',
            visibility='full',
            tracktype='bedGraph')


    def region_track_bb(self):
        """
        Default Track for bigBed

        # CompositeTracks compose different ViewTracks.
        # one for signal in bigWig, another one for bigBed regions.
        """
        return trackhub.ViewTrack(
            name='regionsviewtrack',
            short_label='Regions',
            view='regions',
            visibility='dense',
            tracktype='bigBed')


    def add_tracks(self):
        """
        Composite Tracks compose different ViewTracks
        We will add ViewTrack for signal (bigWig) and for regions (bigBed)
        add bigWig, bigBed
        """
        signal_view = self.signal_track_bw()
        region_view = self.region_track_bb()
        self.composite.add_view(signal_view)
        self.composite.add_view(region_view)

        ## add files
        for f in sorted(self.bw_files + self.bb_files):
            f_track = TrackFile(f).track # create tracks
            if os.path.splitext(f)[1] in ['.bw', '.bigWig', '.bigwig']:
                signal_view.add_tracks(f_track)
            else:
                region_view.add_tracks(f_track)


    def run(self):
        # Create tracks
        self.add_tracks()

        # create directory rights: 711
        if not os.path.exists(self.remote_dir):
            try:
                os.makedirs(self.remote_dir, 0o711)
            except IOError:
                log.error('Create directory failed: {}'.format(self.remote_dir))

        # copy files
        if not self.dry_run:
            trackhub.upload.upload_hub(self.hub, host='localhost', 
               remote_dir=self.remote_dir)

        # save urls
        url_file = os.path.join(self.remote_dir, 'config.json')
        url_dict = {
            'hub_url': self.hub_url,
            'trackhub_url': self.trackhub_url
        }
        Json(url_dict).writer(url_file)

        # print to screen
        log.info('hub_url: {}'.format(self.hub_url))
        log.info('trackhub_url: {}'.format(self.trackhub_url))


def main():
    args = vars(get_args())

    # print_dict(args)
    a = TrackHub(**args)
    # print_dict(a.__dict__)
    a.run()


if __name__ == '__main__':
    main()



# class TrackFile(object):
#     """Extract properties of files for TrackHub, bigWig, bigBed, ...
#     filename, 
#     strand
#     replicates 
#     antibody
#     cell-line
#     """
#     def __init__(self, x):
#         self.file = x
#         self.strand = self.get_strand()
#         self.strand_color = self.get_strand_color()
#         self.ext = self.get_ext()
#         self.seqtype = self.get_seqtype()
#         self.subgroups = self.get_subgroups()


#     def get_strand(self):
#         """Determine the strand of file based on filename
#         watson, forward, fwd, plus
#         crick, reverse, rev, minus
#         """
#         fname = os.path.splitext(self.file)[0].lower()

#         x = re.split('[._]', fname)[-1] # in the tails

#         if x in ['+', 'forward', 'fwd', 'plus', 'watson']:
#             strand = 'f'
#         elif x in ['-', 'reverse', 'rev', 'minus', 'crick']:
#             strand = 'r'
#         else:
#             strand = 'other'

#         return strand


#     def get_strand_color(self):
#         """Pick color based on the filename
#         f:      forward, '#2E3440'
#         r:      reverse, '#6DBFA7'
#         other:  others, '#lalala'

#         example: 
#         rnaseq.fwd.bigWig
        
#         to-do:
#         assign multiple colors for tracks/files    
#         """
#         # Due to how code is extracted from the docs and run during tests, we
#         # need to import again inside a function. You don't normally need this.
#         colors = {
#             'f': '#2E3440',
#             'r': '#6DBFA7',
#             'other': '#1a1a1a',
#         }
#         return trackhub.helpers.hex2rgb(colors[self.strand])


#     def get_ext(self):
#         """Determine the type of file
#         bigWig: bigWig, bigwig, bw
#         bigBed: bigBed, bigbed, bb
#         bedGraph: bedGraph, bedgraph, bg
#         """
#         fext = os.path.splitext(self.file)[1].lower().lstrip('.')

#         if fext in ['bigwig', 'bw']:
#             ext = 'bigWig'
#         elif fext in ['bigbed', 'bb']:
#             ext = 'bigBed'
#         elif fext in ['bedgraph', 'bg']:
#             ext = 'bedGraph'
#         else:
#             ext = 'bigWig'

#         return ext


#     def get_seqtype(self):
#         """Determine the sequencing type of the file
#         ATACseq, RNAseq, ChIPseq
#         default: ChIPseq
#         """
#         f = self.file.lower().replace('-', '')

#         if 'atacseq' in f:
#             seqtype = 'ATACseq'
#         elif 'rnaseq' in f:
#             seqtype = 'RNAseq'
#         elif 'chipseq' in f:
#             seqtype = 'ChIPseq'
#         else:
#             seqtype = 'HiSeq'

#         return seqtype


#     def subgroups_ATACseq(self):
#         """Pick subgroups based on the filename
#         ATACseq: gene, stage, kind
#         RNAseq: gene, stage, strand
#         ChIPseq: rep, strand, kind

#         # to-do
#         ChIPseq: gene, cell-line, kind

#         This functions figures out subgroups based on the number in the
#         filename.  Subgroups provided to the Track() constructor is
#         a dictionary where keys are `rep` attributes from the subgroups added
#         to the composite above, and values are keys of the `mapping` attribute
#         of that same subgroup.

#         Might be easier to cross-reference with the subgroups above, but an
#         example return value from this function would be:
#         """
#         ## ATACseq
#         ## gene, stage, kind
#         fname = os.path.splitext(self.file)[0].lower()
#         ## shGene_
#         p1 = re.compile('gal4x.?sh([A-Za-z0-9]+(\_[^0-9][A-Za-z0-9]+)?)', flags=re.IGNORECASE)
#         p1x = p1.search(fname)
#         if p1x is None:
#             gene = 'gene'
#         else:
#             gene = p1x[1]

#         ## stage (3h, ... L1, )
#         p2 = re.compile('[-_]?(\d+h|L[123]|\d+h?_\d+h)_?', flags=re.IGNORECASE)
#         p2x = p2.search(fname)
#         if p2x is None:
#             stage = 'stage'
#         else:
#             stage = p2x[1]

#         ## kind
#         kind = 'peak' if self.ext == 'bigBed' else 'signal'

#         return {
#             'gene': gene,
#             'stage': stage,
#             'kind': kind
#         }


#     def subgroups_RNAseq(self):
#         """Pick subgroups based on the filename
#         ATACseq: gene, stage, kind
#         RNAseq: gene, stage, strand
#         ChIPseq: rep, strand, kind
#         """
#         ## RNAseq
#         ## gene, stage, strand
#         fname = os.path.splitext(self.file)[0].lower()
#         ## shGene_
#         p1 = re.compile('gal4x.?sh([A-Za-z0-9]+(\_[^0-9][A-Za-z0-9]+)?)', flags=re.IGNORECASE)
#         p1x = p1.search(fname)
#         if p1x is None:
#             gene = 'gene'
#         else:
#             gene = p1x[1]

#         ## stage (3h, ... L1, )
#         p2 = re.compile('[-_]?(\d+h|L[123]|\d+h?_\d+h)_?', flags=re.IGNORECASE)
#         p2x = p2.search(fname)
#         if p2x is None:
#             stage = 'stage'
#         else:
#             stage = p2x[1]

#         ## strand
#         return {
#             'gene': gene,
#             'stage': stage,
#             'strand': self.strand
#         }


#     def subgroups_HiSeq(self):
#         """Pick subgroups based on the filename
#         ATACseq: gene, stage, kind
#         RNAseq: gene, stage, strand
#         ChIPseq: rep, strand, kind
#         HiSeq: rep, strand, kind
#         """
#         ## replicate
#         p1 = re.compile('[-_.](rep|r)(\d+)', flags=re.IGNORECASE)
#         p1x = p1.search(fname)
#         if p1x is None:
#             rep = 'merged'
#         else:
#             rep = p1x[2] # 1, 2

#         ## kind
#         kind = 'peak' if self.ext == 'bigBed' else 'signal'

#         return {
#             'rep': rep,
#             'strand': self.strand,
#             'kind': kind
#         }


#     def get_subgroups(self):
#         if self.seqtype in ['ATACseq']:
#             subgroups = self.subgroups_ATACseq()
#         elif self.seqtype in ['RNAseq']:
#             subgroups = self.subgroups_RNAseq()
#         elif self.seqtype in ['HiSeq', 'ChIPseq']:
#             subgroups = self.subgroups_HiSeq()
#         else:
#             subgroups = self.subgroups_HiSeq()

#         return subgroups


#     def to_track(self):
#         """
#         Create Track for file
#         bigWig, bigBed
#         """
#         if self.ext == 'bigWig':
#             track = trackhub.Track(
#                 name=trackhub.helpers.sanitize(os.path.basename(self.file)),
#                 source=self.file,
#                 visibility='full',
#                 viewLimits='0:10',
#                 maxHeightPixels='8:40:128',
#                 subgroups=self.subgroups,
#                 color=self.strand_color,
#                 tracktype=self.ext)
#         elif self.ext == 'bigBed':
#             track = trackhub.Track(
#                 name=trackhub.helpers.sanitize(os.path.basename(self.file)),
#                 source=self.file,
#                 visibility='dense',
#                 subgroups=self.subgroups,
#                 color=self.strand_color,
#                 tracktype='bigBed')
#         else:
#             raise ValueError('unknown filetype, [bigWig, bigBed] expected, got {}'.format(self.ext))

#         return track


# class Subgroup(object):
#     """
#     group: ATACseq, RNAseq, ChIPseq, HiSeq, ...
#     config: json file, pre-defined subgroups, 

#     ## Sub groups for the project:
#     ## ATACseq: gene, stage, kind
#     ## RNAseq: gene, stage, strand
#     ## ChIPseq: rep, strand, kind

#     ## check the dimension="dimX=strand dimY=rep dimA=kind"
#     """
#     def __init__(self, group="default", config=None, **kwargs):
#         self.update(kwargs, force=True)
#         self.group = group
#         self.config = config
#         self.groups = self.load_group() # groups
#         self.group_list = list(self.groups.keys())
#         self.subgroup = self.get_subgroup() # after load_group()
#         self.subgroup_keys = list(self.subgroup.keys())
#         print('!A-1', self.subgroup, self.subgroup_keys)
#         self.dimensions = 'dimX={} dimY={} dimA={}'.format(
#             self.subgroup_keys[0],
#             self.subgroup_keys[1],
#             self.subgroup_keys[2])


#     def update(self, d, force=True, remove=False):
#         """
#         d: dict
#         force: bool, update exists attributes
#         remove: bool, remove exists attributes
#         Update attributes from dict
#         force exists attr
#         """
#         # fresh start
#         if remove is True:
#             for k in self.__dict__:
#                 # self.__delattr__(k)
#                 delattr(self, k)
#         # add attributes
#         if isinstance(d, dict):
#             for k, v in d.items():
#                 if not hasattr(self, k) or force:
#                     setattr(self, k, v)    


#     def default(self):
#         #  rep strand kind
#         df = {
#             "rep": {
#                 "name": "rep",
#                 "label": "Replicate",
#                 "mapping": {
#                     "0": "merged",
#                     "1": "rep1",
#                     "2": "rep2",
#                     "3": "rep3",
#                     "4": "rep4",
#                     "5": "rep5",
#                     "6": "rep6"
#                 }
#             },
#             "strand": {
#                 "name": "strand",
#                 "label": "Strand",
#                 "mapping": {
#                     "f": "fwd",
#                     "r": "rev",
#                     "other": "non"
#                 }
#             },
#             "kind": {
#                 "name": "kind",
#                 "label": "kind",
#                 "mapping": {
#                     "signal": "signal",
#                     "peak": "peak"
#                 }
#             }
#         }

#         return df
        

#     def load_group(self):
#         if self.config is None:
#             groups = self.default()
#         else:
#             groups = Json(self.config).reader()

#         return groups


#     def example(self):
#         """
#         Create subgroups, from input para
#         name, label, mapping
#         """
#         df = {
#             "example": {
#                 "rep": {
#                     "name": "rep",
#                     "label": "Replicate",
#                     "mapping": {
#                         "1": "rep1",
#                         "2": "rep2"
#                     }
#                 },
#                 "strand": {
#                     "name": "strand",
#                     "label": "Strand",
#                     "mapping": {
#                         "f": "fwd",
#                         "r": "rev"
#                     }
#                 },
#                 "kind": {
#                     "name": "kind",
#                     "label": "kind",
#                     "mapping": {
#                         "signal": "signal",
#                         "peak": "peak"
#                     }
#                 }
#             }
#         }

#         tmp = Json(df).writer()

#         # read text
#         with open(tmp) as r:
#             return r.read()


#     def get_subgroup(self, group=None):
#         if group is None:
#             group = self.group
#         # self.load_group() # self.groups
#         return self.groups.get(group, self.default())


#     def to_trackhub(self):
#         """
#         Convert subgroup to trackhub.SubGroupDefinition()
#         """
#         return [trackhub.SubGroupDefinition(
#                         name = v.get('name', None),
#                         label = v.get('label', None),
#                         mapping = v.get('mapping', None)) 
#                 for k, v in self.subgroup.items()]



# class HubUrl(object):
#     """
#     Construct Hub_url: http://path/hub.txt

#     input: local_path: hub.txt

#     output: http: hub.txt
#     """
#     def __init__(self, url, http_config=None, host=None, root_dir=None):
#         self.url = url
#         self.http_config = http_config
#         self.host = host
#         self.root_dir = root_dir
#         self.config() # update url


#     def config(self):
#         """
#         Determine the url, from host or local directory

#         1. update url from local path
#         """
#         if self.url is None:
#             raise ValueError('failed, url illegal: str expected, got None type')

#         if os.path.exists(self.url):
#             # url is the local path to hub.txt
#             log.info('Convert path to hub_url: {}'.format(self.url))
#             self.url = os.path.abspath(self.url) # update to absolute path
#             self.local_path = self.url
#             self.from_local() # update self.url

#         # check url
#         if self.url is None:
#             raise ValueError('failed, url illegal: str expected, got None type')
#         elif isinstance(self.url, str):
#             if not self.url.split(':')[0] in ['http', 'https', 'ftp']:
#                 raise ValueError('failed, url illegal: [http|ftp] expected, got {}:'.format(self.url))
#         else:
#             raise ValueError('failed, url illegal: str expected, got {}'.format(type(self.url).__name__))


#     def read_url(self, url=None):
#         """
#         hub.txt             
#         hub, shortLabel, longLabel, genomesFile, descriptionUrl

#         genomes.txt
#         genome, trackDb
#         """
#         if url is None:
#             url = self.url

#         ## read url content
#         try:
#             html = urlopen(url).read()
#             text = BeautifulSoup(html, features='html.parser')
#             x = re.split('\s', text.text)
#             xit = iter(x)
#             d = dict(zip(xit, xit))
#             # add attr 
#             for k, v in d.items():            
#                 if not hasattr(self, k):
#                     setattr(self, k, v)
#         except:
#             log.error('reading url failed: {}'.format(url))


#     def hub_checker(self, hubCheck=None):
#         """Check the hubUrl by hubCheck
#         specify the path of hubCheck
#         """
#         if hubCheck is None:
#             hubCheck = shutil.which('hubCheck') # from env
#         cmd = hubCheck + ' ' + self.url

#         log.info('Checking Track Hub: {}'.format(self.url))
#         try:
#             run_shell_cmd(cmd)
#             return True
#         except:
#             log.error('TrackHub failed: {}'.format(self.url))


#     def get_genome(self):
#         """
#         Get the genome from hubUrl, genomes.txt
#         """
#         self.read_url() # read url, update genomesFile
#         genome_txt = getattr(self, 'genomesFile', None)

#         if not genome_txt is None:
#             genome_url = os.path.join(os.path.dirname(self.url), genome_txt)
#             self.read_url(genome_url) # update genome

#         # update genome
#         return getattr(self, 'genome', None) # update genome


#     def add_position(self, position=None):
#         """
#         Convert position to url format
#         Input: chr1:1-1000
#         Output: &position=chr1%3A1-1000

#         see: UCSC Documentation
#         https://genome.ucsc.edu/goldenPath/help/customTrack.html#optParams
#         """
#         if isinstance(position, str):
#             position = re.sub('[^\w:-]', '', position) # remove characters
#             p = re.compile('^(\w+):([\d,]+)-([\d,]+)$')
            
#             if p.search(position):
#                 chr, start, end = re.split('[:-]', position)
#                 suffix = '&position={}%3A{}-{}'.format(chr, start, end)
#             else:
#                 log.warning('illegal position, chr1:1-100 expected, got {}'.format(position))
#                 suffix = ''
#         else:
#             log.warning('illegal position, str expected, got {}'.format(type(position).__name__))
#             suffix = ''

#         return self.url + suffix


#     def from_local(self):
#         """
#         Convert the local path to hub_txt
#         x, the directory of the hub.txt file
#         config is a json file, contains: 'host' and 'root_dir' at least
        
#         host: the ip or domain of the host, eg: 123.123.123.123 (http)
#         root_dir; the root dir of the website

#         return the url hub.txt
#         """
#         # read {host|root_dir} from config
#         if not self.http_config is None:
#             args = Json(self.http_config).reader()
#             self.host = args.get('host', None)
#             self.root_dir = args.get('root_dir', None)

#         # check {host|root_dir}
#         if self.host is None or self.root_dir is None:
#             log.warning('http_config failed: host={}, root_dir={}'.format(self.host, self.root_dir))
#             return None # failed

#         # trim the root_dir from x
#         if os.path.exists(self.url) and self.url.startswith(self.root_dir):
#             # Trim the common prefix from x
#             url_new = self.url.replace(self.root_dir, '')
#             url_new = url_new.lstrip('/')
#             # update the url
#             # hub_url = host + '/' + x2
#             self.url = self.host + '/' + url_new
#         else:
#             log.error('file path not consistent with root_dir: {}, {}'.format(self.url, self.root_dir))


#     def to_trackhub(self, region='usa', mirror=None, position=None):
#         """Convert hub.txt url to trackhub Url for sharing
#         hub.txt: http://ip/hub.txt
#         trackhub_url: http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://ip/hub.txt 

#         region: usa, asia, euro

#         usa:  http://genome.ucsc.edu
#         asia: http://genome-asia.ucsc.edu
#         euro: http://genome-euro.ucsc.edu

#         """
#         regions = {
#             'usa': 'http://genome.ucsc.edu',
#             'asia': 'http://genome-asia.ucsc.edu',
#             'euro': 'http://genome-euro.ucsc.edu'}

#         # determine the mirror
#         base_url = None
#         if not mirror is None:
#             if mirror.split(':')[0] in ['http', 'https']:
#                 base_url = mirror
#         elif not region is None:
#             region = region.lower()
#             base_url = regions.get(region, None)
#         else:
#             pass

#         ## add cgi
#         if base_url is None:
#             raise ValueError('mirror failed: ucsc_url: {}\
#                 mirror={}, \
#                 region={}'.format(base_url, mirror, region))

#         ## add genome
#         trackhub_url = '{}/cgi-bin/hgTracks?db={}'.format(base_url, self.get_genome())
        
#         ## add trackhub url
#         if not position is None:
#             hub_url = self.add_position(position)
#         else:
#             hub_url = self.url

#         trackhub_url += '&hubUrl={}'.format(hub_url)

#         return trackhub_url


# class TrackHubPort(object):
#     """
#     The function/class for all trackhub generator

#     data,   
#     genome, hg19, hg38, mm9, mm10, dm3, dm6
#     hub_name, 
#     user,
#     email
#     remote dir
#     config
#     dry-run
#     """
#     def __init__(self, **kwargs):
#         self.args = kwargs
#         self.update(kwargs, force=True) # fresh new
#         self.config()
#         self.init_hub() # 
#         self.composite = self.init_composite(group=self.seqtype, config=self.subgroups_config)


#     def update(self, d, force=True, remove=False):
#         """
#         d: dict
#         force: bool, update exists attributes
#         remove: bool, remove exists attributes
#         Update attributes from dict
#         force exists attr
#         """
#         # fresh start
#         if remove is True:
#             for k in self.__dict__:
#                 # self.__delattr__(k)
#                 delattr(self, k)
#         # add attributes
#         if isinstance(d, dict):
#             for k, v in d.items():
#                 if not hasattr(self, k) or force:
#                     setattr(self, k, v)    


#     def config(self):
#         args_init = {
#             'data_dir': str(pathlib.Path.cwd()),
#             'user': 'UCSC',
#             'email': 'abc@abc.com',
#             'short_label': 'hiseq',
#             'long_label': 'hiseq',
#             'descriptionUrl': None
#         }
#         self.update(args_init, force=False)

#         # 1. labels
#         if self.short_label is None:
#             self.short_label = self.hub_name
#         if self.long_label is None:
#             self.long_label = self.short_label

#         # 2. list files
#         data_dir_list = [i for i in listdir(self.data_dir, include_dir=True) if os.path.isdir(i)]

#         self.bw_files = []
#         self.bb_files = []

#         # list all bw files, 
#         [self.bw_files.extend(self.get_bw_files(i)) for i in data_dir_list]

#         # list all bb files
#         [self.bb_files.extend(self.get_bb_files(i)) for i in data_dir_list]

#         # 3. seqtype
#         if len(self.bw_files) > 0:
#             self.seqtype = TrackFile(self.bw_files[0]).seqtype # by the first one
#         else:
#             self.seqtype = 'HiSeq'

#         # 4. subgroups config
#         self.subgroup_config = '/data/biodata/mydb/subgroups_config.json'
#         self.http_config = '/data/biodata/mydb/http_config.json'

#         # 5. read config
#         http_df = Json(self.http_config).reader()
#         self.http_host = http_df.get('host', None)
#         self.http_root_dir = http_df.get('root_dir', None)


#     def init_hub(self):
#         self.hub, self.genome_file, _, self.trackdb = trackhub.default_hub(
#             hub_name=self.hub_name,
#             short_label=self.short_label,
#             long_label=self.long_label,
#             genome=self.genome,
#             email=self.email,
#             descriptionUrl=self.descriptionUrl)


#     def init_composite(self, group, config):
#         """Initiate the composite track
#         default
#         frame
#         """
#         # determine subgroup by seqtype
#         subgroup = Subgroup(group=self.seqtype, config=self.subgroup_config)

#         ## basic composite track 
#         # Create the composite track
#         composite = trackhub.CompositeTrack(
#             name='composite',
#             short_label=self.short_label,
#             long_label=self.short_label,
#             dimensions=subgroup.dimensions,
#             filterComposite='dimA',
#             # sortOrder='rep=+ kind=-',
#             tracktype='bigWig',
#             visibility='full'
#         )

#         # Add subgroups to the composite track
#         composite.add_subgroups(subgroup.to_trackhub())

#         # Add the composite track to the trackDb
#         self.trackdb.add_tracks(composite)

#         return composite # CompositeTrack


#     def get_bw_files(self, path):
#         """List the track files, bigWig, bigbed, bedgraph
#         bigWig, bw, bigwig, bigWig
#         bigBed, bb, bigbed, bigBed
#         # bedGraph, bg, bedgraph, bedGraph
#         """
#         bw1 = listfile(path, "*.bigWig")
#         bw2 = listfile(path, "*.bigwig")
#         bw3 = listfile(path, "*.bw")

#         return bw1 + bw2 + bw3


#     def get_bb_files(self, path):
#         """List the track files, bigWig, bigbed, bedgraph
#         bigWig, bw, bigwig, bigWig
#         bigBed, bb, bigbed, bigBed
#         # bedGraph, bg, bedgraph, bedGraph

#         """
#         bb1 = listfile(path, "*.bigBed")
#         bb2 = listfile(path, "*.bigbed")
#         bb3 = listfile(path, "*.bb")

#         return bb1 + bb2 + bb3


#     def bw2bb(self, bw):
#         """
#         get the bb file by bw name, prefix
#         """
#         bw_name = os.path.basename(os.path.splitext(bw)[0])
#         # [i for i in self.bb_files if os.path.basename(i).startswith(bw_name)]
#         bb_list = [os.path.splitext(i)[0] + '.bigBed' for i in bw_name]
#         return [i for i in self.bb_files if os.path.basename(i) == bw_name]


#     def bb2bw(self, bb):
#         """
#         get the bb file by bw name, prefix
#         """
#         bb_name = os.path.basename(os.path.splitext(bb)[0])
#         # return [i for i in self.bw_files if os.path.basename(i).startswith(bb_name)]
#         bw_list = [os.path.splitext(i)[0] + '.bigWig' for i in bb_name]
#         return [i for i in self.bw_files if os.path.basename(i) == bb_name]


#     def signal_track_bw(self):
#         """
#         Default Track for bigWig

#         # CompositeTracks compose different ViewTracks.
#         # one for signal in bigWig, another one for bigBed regions.
#         """
#         return trackhub.ViewTrack(
#             name='signalviewtrack',
#             short_label='Signal',
#             view='signal',
#             visibility='full',
#             tracktype='bigWig')


#     def signal_track_bg(self):
#         """
#         Default Track for bedGraph

#         # CompositeTracks compose different ViewTracks.
#         # one for signal in bigWig, another one for bigBed regions.
#         """
#         return trackhub.ViewTrack(
#             name='singalviewtrack2',
#             short_label='Signal2',
#             view='signal',
#             visibility='full',
#             tracktype='bedGraph')


#     def region_track_bb(self):
#         """
#         Default Track for bigBed

#         # CompositeTracks compose different ViewTracks.
#         # one for signal in bigWig, another one for bigBed regions.
#         """
#         return trackhub.ViewTrack(
#             name='regionsviewtrack',
#             short_label='Regions',
#             view='regions',
#             visibility='dense',
#             tracktype='bigBed')


#     def add_tracks(self):
#         # Define the tracks
#         if len(self.bw_files) < 1:
#             raise ValueError('bigWig files not found: {}'.format(self.data_dir))
        
#         # Composite Tracks compose different ViewTracks
#         # We will add ViewTrack for signal (bigWig)
#         # and for regions (bigBed)
#         ## add bigWig, bigBed
#         signal_view = self.signal_track_bw()
#         region_view = self.region_track_bb()
#         self.composite.add_view(signal_view)
#         self.composite.add_view(region_view)

#         bb_out = self.bb_files
#         for bw in self.bw_files:
#             signal_view.add_tracks(TrackFile(bw).to_track())
#             bb_list = self.bw2bb(bw)
#             if isinstance(bb_list, list):
#                 if len(bb_list) > 0:
#                     for bb in bb_list:
#                         region_view.add_tracks(TrackFile(bb).to_track())
#                         bb_out.remove(bb)

#         ## add bigBed (extra)
#         # region_view = self.region_track_bb()
#         # self.composite.add_view(region_view)
#         # if len(self.bb_files) > 0:
#         #     for bb in self.bb_files:
#         #         x = TrackFile(bb).ext
#         #         region_view.add_tracks(TrackFile(bb).to_track())

#         if len(bb_out) > 0:
#             for bb in bb_out:
#                 region_view.add_tracks(TrackFile(bb).to_track())


#     def run(self):
#         self.add_tracks()

#         # create directory rights: 711
#         if not os.path.exists(remote_dir):
#             try:
#                 os.makedirs(remote_dir, 0o711)
#             except IOError:
#                 log.error('Create directory failed: {}'.format(remote_dir))
        
#         # copy files
#         trackhub.upload.upload_hub(self.hub, host='localhost', remote_dir=remote_dir)

#         # hub.txt to url
#         hub_url = HubUrl(hub_txt,  http_config=self.http_config).url

#         # convert to trackhub url
#         trackhub_url = HubUrl(hub_url).to_trackhub()

#         print('hub_url:', hub_url)
#         print('trackhub_url:', trackhub_url)


# def local_to_hub(x, config=None, host=None, root_dir=None, hub_check=False):
#     """
#     Convert the local path to hub_txt

#     x absolute path of hub.txt, eg: /data/hub.txt 
#     host the ip or domain of the host, eg: 123.123.123.123 (http)
#     root_dir the root dir of the website

#     or read the args from config

#     return the list of hub.txt
#     """
#     # read from config
#     if not config is None:
#         args = Json(config).reader()
#         host = args.get('host', None)
#         root_dir = args.get('root_dir', None)

#     # check input
#     if host is None or root_dir is None:
#         return None # failed

#     # turn the x to absolute path
#     x = os.path.abspath(x)

#     # trim the root_dir from x
#     if x.startswith(root_dir):
#         x2 = x.replace(root_dir, '')
#         x2 = x2.lstrip('/')
#         hub_url = host + '/' + x2

#         # check
#         if hub_check:
#             hub_checker(hub_url)

#         return hub_url
#     else:
#         log.error('file path not consistent with root_dir: {}, {}'.format(x, root_dir))


# def pos_formatter(pos):
#     """Convert position to url format
#     Input: chr1:1-1000
#     Output: &position=chr1%3A1-1000
#     """
#     assert isinstance(pos, str)
#     pos = re.sub('[^\w:-]', '', pos) # remove characters
#     p = '^chr\w+\:\d+\-\d+$' # pattern
    
#     if p.fullmatch(pos):
#         chr, start, end = re.spilt('[:-]', pos)
#         return '{}:{}-{}'.format(chr, start, end)
#     else:
#         return None


# def hub_checker(url, hubCheck=None):
#     """Check the hubUrl by hubCheck"""
#     if hubCheck is None:
#         hubCheck = shutil.which('hubCheck') # from env
#     cmd = hubCheck + ' ' + url

#     log.info('Checking Track Hub: {}'.format(url))
#     try:
#         run_shell_cmd(cmd)
#         return True
#     except:
#         log.error('TrackHub failed: {}'.format(url))



# def list_local_hub(x, config, recursive=False, check=False):
#     """
#     List the hub.txt files from the path
#     config, json file, save "host" and "root_dir" of the server

#     return the url of the hub
#     """
#     # all hub.txt files
#     x_list = listfile(x, "*hub.txt", recursive=recursive)

#     # convert to url
#     hub_list = [local_to_hub(i, config) for i in x_list]
    
#     # check hub
#     hub_status = [hub_checker(i) for i in hub_list]
#     hub_ok = ['ok' if i else 'failed' for i in hub_status]

#     for url, status in zip(hub_list, hub_ok):
#         log.info('Trackhub : {:6s} : {}'.format(status, url))

#     if check:
#         hub_list = [i for i,j in zip(hub_list, hub_status) if j]

#     return hub_list


# def hub_to_trackhub(hub_url, genome, position=None):
#     """Create url for hub.txt, for sharing
#     hub_txt, the accessiable url of hub.txt file
#     genome, UCSC supported name of genome, eg: hg19, dm3
#     position, chr1:1-1000

#     return http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://url-to-hub.txt 
#     """
#     assert isinstance(hub_url, str)
#     assert isinstance(genome, str)
#     base_url = 'http://genome.ucsc.edu'
#     trackhub_url = '{}/cgi-bin/hgTracks?db={}'.format(base_url, genome)
#     trackhub_url += '&hubUrl={}'.format(hub_url)

#     if not position is None:
#         position = pos_formatter(position)
#         if position:
#             trackhub_url += '&position={}'.format(position)

#     return trackhub_url

