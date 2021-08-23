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
import re
import signal
import time
import logging
import trackhub
import tempfile
import argparse
import subprocess
from shutil import which
from urllib.request import urlopen
from bs4 import BeautifulSoup
from hiseq.utils.file import file_exists, file_abspath, list_dir, \
    check_path, remove_path
from hiseq.utils.utils import log, update_obj, Config, get_date, run_shell_cmd
from hiseq.get_trackhub.hub_url import HubUrl
from hiseq.get_trackhub.http_server import HttpServer
from hiseq.get_trackhub.track_file import TrackFile
from hiseq.get_trackhub.utils import trackhub_config_template


## run
def get_ticks():
    """Returns ticks.
       Mark the time stamp
        - Python3: Use time.perf_counter().
    """
    return getattr(time, 'perf_counter', getattr(time, 'time'))()


# Stage4. Auto recognize composite tracks, subgroups
#
# 1. if subgroups.json exists, the files in the folder should organize into one
#    compositeTrack()
#
# 2. According to TrackFile, add signal_view or region_view, or both
#
# 3. Create, Init compositeTrack, based on subgroups


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
        # required arguments:
        # 1. subgroups, YAML file, saving subgroups structure
        # 2. dimensions, dimX, dimY, ...
        # 3. sortOrder, gene=+ stage=+ (default: dimX, dimY, ...)
        # 4. filterComposite, dimA
        # load subgroups.yaml to dict
        # dimX, dimY, dimA, ...
        subgroups_dict = Config().load(composite_subgroup)

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
        """Sanitize a string, for shortLabel, letters, numbers 
        allowed only [A-Za-z0-9]

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
                    print('!AAAA-2', self.colorPal)
                    tk = TrackFile(bw,
                        subgroups=subgroups_dict,
                        colorDim=self.colorDim, # dimX, dimY, dimA, ...
                        colorPal=self.colorPal,
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
            check_path(self.remote_hub_dir, create_dirs=True)
            # construct tracks()
            tmp_dir = self.render()
            remove_path(tmp_dir, ask=False)
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
        Config().dump(hub_dict, hub_yaml)
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
        Config().dump(hub_dict, hub_yaml)
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
            args_config = Config().load(self.config)
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
            d_lv1 = [i for i in list_dir(self.data_dir, include_dir=True)
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
            A yaml file, specify the subgroups for the track files in the dir
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
                    f_lv1 = [i for i in list_dir(data_dir)
                        if TrackFile(i).is_track_file()]
                    # files: level=2
                    d_lv1 = [i for i in list_dir(data_dir, include_dir=True)
                        if os.path.isdir(i)]
                    f_lv2 = [] # empty
                    if len(d_lv1) > 0:
                        [f_lv2.extend([i for i in list_dir(d)
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


def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog='get_trackhub',
        description='Generate trackhub for bigWig and bigBed files',
        epilog='Example: \n\
               python get_trackhub.py --config config.yaml')
    parser.add_argument('-d', '--demo', action='store_true',
        help='Show the tutorial')
    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
        help='Generate trackhub files, Do not copy track files to remote_dir')
    parser.add_argument('-c', '--config', default=None,
        help='Config file, run -d, generate the template')
    return parser
    
    
def main():
    args = vars(get_args().parse_args())
    TrackHub(**args).run()


if __name__ == '__main__':
    main()

#