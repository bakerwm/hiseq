#!/usr/bin/env python
"""
This script is designed to create track_hub for bigWig + bigBed
composite tracks
make subgroups: strand (fwd/rev), replicate (if spcified), kind (signal, peak)


"""

import glob, os, re, argparse
import subprocess, shlex, logging
import trackhub

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog='get_trackhub',
        description='Generate trackhub for bigWig and bigBed files',
        epilog='Example: \
               python ~/work/wmlib/wmlib/make_trackhub.py --source clip_bigWig/ --hub PTB_CLIP \
               --genome hg19 --server yylab_tracks ')
    parser.add_argument('--source', required=True, dest='source',
        help='The directory of bigWig and bigBed files')
    parser.add_argument('--hub', required=True, metavar='hub_name',
        default=None, help='hub name (remove blanks)')
    parser.add_argument('--genome', required=True, metavar='GENOME', 
        default='hg19', help='genome for the trackhub')
    parser.add_argument('--short-label', default=None, dest='short_label',
        help='short label for the hub, default: [--hub]')
    parser.add_argument('--long-label', default=None, dest='long_label',
        help='long label for the hub, default: [--hub]')
    parser.add_argument('--user', default='Ming Wang',
        help='Who maintain the trackhub')
    parser.add_argument('--email', default='wangm08@hotmail.com',
        help='email for the hub')
    parser.add_argument('--server', default='yulab_tracks',
        help='Name of directory to save track files, default: [yulab_tracks]')
    parser.add_argument('--server-dir', dest='server_dir',
        default='/data/public/upload/ucsc_trackhub/',
        help='Directory to save the track files, default: /data/public/upload/ucsc_trackhub/')
    parser.add_argument('--server-url', default='http://159.226.118.232/upload/ucsc_trackhub',
        help='Server to deposite the data, default: http://159.226.118.232/upload/ucsc_trackhub')
    args=parser.parse_args()
    return args


##----------------------------------------------------------------------------##
## helper functions
# snippet is placed into public domain by
# anatoly techtonik <techtonik@gmail.com>
# http://stackoverflow.com/questions/8151300/ignore-case-in-glob-on-linux

import fnmatch
import os
import re

def findfiles(which, where='.'):
    """Returns list of filenames from `where` path matched by 'which'
    shell pattern. Matching is case-insensitive.
    # findfiles('*.ogg')
    """    
    # TODO: recursive param with walk() filtering
    rule = re.compile(fnmatch.translate(which), re.IGNORECASE)
    fn_names = [name for name in os.listdir(where) if rule.match(name)]
    return [os.path.join(where, f) for f in fn_names]


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def hub_url_checker(url):
    """Check the hubUrl by hubCheck"""
    hubCheck = which('hubCheck')
    if hubCheck is None:
        logging.info('hubCheck not found, please check the url manually.')
        return False
    else:
        p1 = subprocess.run(shlex.split('%s %s' % (hubCheck, url)),
            stderr=subprocess.PIPE)
        if p1.stderr.decode() == '':
            return True
        else:
            logging.error('hub-url error: %s' % url)
            return False


def source_checker(path):
    """Check source directory, whether .bigWig or .bigBed files exists"""
    check_bw1 = glob.glob(os.path.join(path, '*.bigWig'))
    check_bw2 = glob.glob(os.path.join(path, '*.bw'))
    check_bg1 = glob.glob(os.path.join(path, "*.bigBed"))
    check_bg2 = glob.glob(os.path.join(path, "*.bb"))
    check_bw1.extend(check_bw2)
    if len(check_bw1) == 0:
        return False
    else:
        return True


def file_strand(fn):
    """Determine the strand of file based on filename"""
    assert isinstance(fn, str)
    fn_name = os.path.splitext(fn)[0].lower()
    if fn_name.endswith('fwd') or fn_name.endswith('watson') or fn_name.endswith('plus'):
        strand = 'f'
    elif fn_name.endswith('rev') or fn_name.endswith('crick') or fn_name.endswith('minus'):
        strand = 'r'
    else:
        strand = 'other'
    return strand


def pos_to_url(pos):
    """Convert position to url format
    Input: chr1:1-1000
    Output: &position=chr1%3A1-1000
    """
    assert isinstance(pos, str)
    pos = re.sub(',', '', p) # remove comma from string
    flag = re.fullmatch('^chr\w+\:\d+\-\d+$', pos)
    if flag is None:
        return None
    else:
        chr, coords = pos.split(':')
        start, end = coords.split('-')
        pos_url = '&position=%s%3A%s-%s' % (chr, start, end)
        return pos_url


def hub_to_url(hub_txt, genome, position=None):
    """Create url for hub.txt
    hub_txt, the accessiable url of hub.txt file
    genome, UCSC supported name of genome, eg: hg19, dm3
    position, chr1:1-1000
    """
    assert isinstance(hub_txt, str)
    assert isinstance(genome, str)
    base_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks'
    track_url = '%s?db=%s&hubUrl=%s' % (base_url, genome, hub_txt)
    if position is None:
        return track_url
    else:
        pos_url = pos_to_url(position)
        return base_url + pos_url


def get_subgroups():
    # Subgroups provide a way of tagging tracks.
    subgroups = [
        # A subgroup for replicates
        trackhub.SubGroupDefinition(
            name='rep',
            label='rep',
            mapping={
                '0': 'merged',
                '1': 'rep1',
                '2': 'rep2',
                '3': 'rep3',
                '4': 'rep4',
                '5': 'rep5',
                '6': 'rep6',
            }
        ),

        # Forward/Reverse strand
        trackhub.SubGroupDefinition(
            name='strand',
            label='strand',
            mapping={
                'f': 'Forward',
                'r': 'Reverse',
            }
        ),

        # Turning on/off the signal or regions in bulk.
        trackhub.SubGroupDefinition(
            name='kind',
            label='kind',
            mapping={
                'signal': 'signal',
                'peak': 'peak',
            }
        ),
    ]
    return subgroups


# def subgroups_from_filename(fn):
def subgroups_picker(fn):
    """Pick subgroups based on the filename
    This functions figures out subgroups based on the number in the
    filename.  Subgroups provided to the Track() constructor is
    a dictionary where keys are `rep` attributes from the subgroups added
    to the composite above, and values are keys of the `mapping` attribute
    of that same subgroup.

    Might be easier to cross-reference with the subgroups above, but an
    example return value from this function would be:
    """
    assert isinstance(fn, str)
    # 1. strand
    strand = file_strand(fn)

    # 2. kind, peak/signal
    fn_ext = os.path.splitext(fn)[1].lower()
    if fn_ext == '.bw' or fn_ext == '.bigwig':
        kind = 'signal'
    elif fn_ext == '.bb' or fn_ext == '.bigbed':
        kind = 'peak'
    else:
        kind = 'peak'
    
    # 3. replicate
    fn_name = os.path.splitext(fn)[0].lower()
    fn_flag = re.compile('_rep[0-9]').search(fn_name)
    if fn_flag is None:
        rep = 0
    else:
        fn_flag_text = fn_name[fn_flag.start():fn_flag.end()]
        rep = re.sub('_rep', '', fn_flag_text) # the number
    track_subgroup = {
        'rep': rep,
        'strand': strand,
        'kind': kind,
    }
    return track_subgroup


# def color_from_filename(fn):
def color_picker(fn):
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
    import trackhub
    assert isinstance(fn, str)
    strand = file_strand(fn)

    colors = {
        'f': '#2E3440',
        'r': '#6DBFA7',
        'other': '#1a1a1a',
    }
    return trackhub.helpers.hex2rgb(colors[strand])


def create_composite(trackdb, short_label, bigwig_files, 
    bigbed_files=None, bedgraph_files=None):
    """ Create composite tracks, bigWig, bigBed,
    Add CompositeTrack to trackdb
    """
    assert isinstance(trackdb, trackhub.trackdb.TrackDb)
    assert isinstance(short_label, str)
    assert isinstance(bigwig_files, list)

    # Create the composite track
    composite = trackhub.CompositeTrack(
        name='composite',
        short_label=short_label,
        long_label=short_label,
        dimensions='dimX=strand dimY=rep dimA=kind',
        filterComposite='dimA',
        sortOrder='rep=+ kind=-',
        tracktype='bigWig',
        visibility='full'
    )

    # Add subgroups to the composite track
    composite.add_subgroups(get_subgroups())

    # Add the composite track to the trackDb
    trackdb.add_tracks(composite)

    # CompositeTracks compose different ViewTracks.
    # one for signal in bigWig, another one for bigBed regions.
    signal_view = trackhub.ViewTrack(
        name='signalviewtrack',
        short_label='Signal',
        view='signal',
        visibility='full',
        tracktype='bigWig')

    signal_view2 = trackhub.ViewTrack(
        name='singalviewtrack2',
        short_label='Signal2',
        view='signal',
        visibility='full',
        tracktype='bedGraph')

    regions_view = trackhub.ViewTrack(
        name='regionsviewtrack',
        short_label='Regions',
        view='regions',
        visibility='dense',
        tracktype='bigBed')

    # for bigWig files
    # signal view
    composite.add_view(signal_view)
    for bigwig in bigwig_files:
        track=trackhub.Track(
            name=trackhub.helpers.sanitize(os.path.basename(bigwig)),
            source=bigwig,
            visibility='full',
            tracktype='bigWig',
            viewLimits='-2:2',
            maxHeightPixels='8:40:128',
            subgroups=subgroups_picker(bigwig),
            color=color_picker(bigwig))

        # Note that we add the track to the *view* rather than the trackDb as
        # we did in the README example.
        signal_view.add_tracks(track)


    # # for bedgraph
    # if bedgraph_files is None:
    #     pass
    # elif len(bedgraph_files) > 0:
    #     composite.add_view(signal_view2)
    #     for bg in bedgraph_files:
    #         track = trackhub.Track(
    #             name=trackhub.helpers.sanitize(os.path.basename(bg)),
    #             source=bg,
    #             visibility='full',
    #             subgroups=subgroups_picker(bg),
    #             color=color_picker(bg),
    #             tracktype='bedGraph',
    #             viewLimits='-2:2',
    #             maxHeightPixels='8:40:128')
    #         signal_view2.add_tracks(track)
    # else:
    #     pass


    # for bigBed files
    # region view
    # bigbed_files = glob.glob(os.path.join(data_path, '*.bigBed'))
    if bigbed_files is None:
        pass
    elif len(bigbed_files) > 0 :
        composite.add_view(regions_view)
        for bigbed in bigbed_files:
            track = trackhub.Track(
                name=trackhub.helpers.sanitize(os.path.basename(bigbed)),
                source=bigbed,
                visibility='full',
                subgroups=subgroups_picker(bigbed),
                color=color_picker(bigbed),
                tracktype='bigBed')
            regions_view.add_tracks(track)
    else:
        pass

    return composite # CompositeTrack


def main():
    args = get_args()

    # remove spaces for hub_name
    hub_name = re.sub('\s+', '', args.hub)

    # check labels for hub
    if args.short_label is None:
        args.short_label = hub_name
    if args.long_label is None:
        args.long_label = args.long_label

    # check server name
    if args.server is None:
        args.server = hub_name
    
    # server dir, save files
    server_path = os.path.join(args.server_dir, args.server, hub_name)
    if not os.path.exists(server_path):
        try:
            os.makedirs(server_path, 0o711)
        except IOError:
            print('cannot create directory: ' + server_path)

    # source directory
    # option-1
    # source/*.bigWig
    #
    # option-2
    # source/*/*.bigWig
    if source_checker(args.source) is True:
        bw_dirs = [args.source, ]
    else:
        bw_dirs = [d for d in glob.glob(args.source + '/*') if os.path.isdir(d)]
        bw_dirs = [d for d in bw_dirs if source_checker(d) is True]

    # list files
    bigwig_files = []
    bigbed_files = []
    for bw_dir in bw_dirs:
        bw = findfiles('*.bigwig', bw_dir)
        bw.extend(findfiles('*.bw', bw_dir))
        bb = findfiles('*.bigbed', bw_dir)
        bb.extend(findfiles('*.bb', bw_dir))

        bigwig_files.extend(bw)
        bigbed_files.extend(bb)

    if len(bigwig_files) == 0:
        raise ValueError('no *bigWig files found: %s' % args.source)


    #################
    # Create tracks #
    #################
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=args.short_label,
        long_label=args.long_label,
        genome=args.genome,
        email=args.email)

    #########################
    # Create CompositeTrack #
    #########################
    composite = create_composite(
        trackdb=trackdb, 
        short_label=args.short_label, 
        bigwig_files=bigwig_files, 
        bigbed_files=bigbed_files)

    trackhub.upload.upload_hub(hub=hub, host='localhost', 
        remote_dir=server_path)

    # hub_txt
    # server:     hub_name
    # server_dir: /data/ucsc_gb/trackhub_open/
    # server_url: http://159.226.118.232/open/
    # 
    # server_path: /data/ucsc_gb/trackhub_open/<server>/
    # server_url:  http://159.226.118.232/open/<server>/
    # hub_txt_url: http://159.226.118.232/open/<server>/<hub_name>.hub.txt
    hub_txt_url = os.path.join(args.server_url, args.server, hub_name, hub_name + '.hub.txt')
    
    if hub_url_checker(hub_txt_url) is True:
        status = 'ok'
    else:
        status = 'failed'
    print('%10s : %s' % (status, hub_txt_url))


    # Example uploading to a web server (not run):
    if 0:
        trackhub.upload.upload_hub(
            hub=hub, host='example.com', user='username',
            remote_dir='/var/www/example_hub')

if __name__ == '__main__':
    main()


## EOF
