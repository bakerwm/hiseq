#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Generate config.toml template
"""

import os
import tempfile
from hiseq.utils.utils import log, Config, get_date


def trackhub_config_template():
    """Generate config.toml template
    set default values
    save args to file: config_template.toml
    save subgroups to file: subgroups_ytemplate.toml
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
    config_f = 'config_template.toml'
    subgroups_f = 'subgroups_template.toml'
    Config().dump(default_args, config_f)
    Config().dump(default_subgroups, subgroups_f)
    # show message
    cf = '{:>16s} : {}'.format('config.tomll', config_f),
    gf = '{:>16s} : {}'.format('subgroups.toml', subgroups_f),
    msg = '\n'.join([
        '{}'.format('#'*80),
        '# {:<77s}#'.format('Mini-tutorial [run_trackhub]'),
        '# {:<77s}#'.format('1. Generating template files:'),
        '# {:<77s}#'.format('$ python run_trackhub --demo'),
        '# {:<77s}#'.format(str(cf)),
        '# {:<77s}#'.format(str(gf)),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('2. Move the toml files to {data_dir}'),
        '# {:<77s}#'.format('$ mv subgroups.toml {data_dir}/subgroups.toml'),
        '# {:<77s}#'.format('$ mv config.toml {data_dir}/config.toml'),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('Update the toml files, according to your data'),
        '# {:<77s}#'.format('Required fields - config.toml'),
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
        '# {:<77s}#'.format('Required fields - subgroups.toml'),
        '# {:<77s}#'.format('  - mapping'),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('3. Generating trackhub files'),
        '# {:<77s}#'.format('$ hiseq run_trackhub -c {data_dir}/config.toml'),
        '# {:<77s}#'.format(''),
        '# {:<77s}#'.format('4. Find the hub.txt file'),
        '# {:<77s}#'.format('{remote_dir}/hub_name/{hub_name}_hub.txt'),
        '{}'.format('#'*80)
    ])
    print(msg)
    # save to tempfile: get_trackhub.readme.txt
    f = 'get_trackhub.readme.txt'
    if os.path.exists(f):
        t = get_date() # '2021-05-31 01:26:40'
        t = t.replace('-', '').replace(':', '').replace(' ', '_')
        f = 'get_trackhub.readme.{}.txt'.format(t)
    with open(f, 'wt') as w:
        w.write(msg+'\n')


def main():
    trackhub_config_template()


if __name__ == '__main__':
    main()
    
#