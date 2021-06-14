#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Check bacteria content using Kraken2, http://ccb.jhu.edu/software/kraken2/

Notes:
1. Use clean data for assessment. (why?)

1. bowtie method (not recommended)

2. kraken2/centrifuge method

## NCBI Taxonomy classification

## standard output of kraken2
Like Kraken 1, Kraken 2 offers two formats of sample-wide results. 
Kraken 2's standard sample report format is tab-delimited with one line 
per taxon. The fields of the output, from left-to-right, are as follows:

1. Percentage of fragments covered by the clade rooted at this taxon
2. Number of fragments covered by the clade rooted at this taxon
3. Number of fragments assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
   (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
   Taxa that are not at any of these 10 ranks have a rank code that is
   formed by using the rank code of the closest ancestor rank with
   a number indicating the distance from that rank.  E.g., "G2" is a
   rank code indicating a taxon is between genus and species and the
   grandparent taxon is at the genus rank.
5. NCBI taxonomic ID number
6. Indented scientific name

see here for output format: https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats
"""

import os
import sys
import argparse
import pathlib
import subprocess
from multiprocessing import Pool
from shutil import which
import pandas as pd
import hiseq
from hiseq.trim.trim_r1 import TrimR1
from hiseq.utils.fastx import Fastx
from hiseq.utils.utils import update_obj, log, Config, run_shell_cmd
from hiseq.utils.file import file_exists, check_path, file_abspath, fx_name, \
    file_prefix, check_fx, symlink_file


class Kraken2(object):
    """Run taxonomic classification using program kraken2. (version 2.0.7-beta)
    ## extra: remove adapters first
    
    $ kraken2 --db <path-to-db> 
    mission:

    1. Construct database
    download NCBI data (human, bacteria, archaea, and virus)
    require >130 GB disk space to process the database,
    final size is ~100 GB.

    2. Run program

    3. Return report
    # fq = '/home/wangming/work/yu_2019/projects/20190312_zp_goldclip/results/goldclip_v2/06.unmap_reads/demo.fq'
    # out = '/home/wangming/work/yu_2019/projects/20190312_zp_goldclip/results/goldclip_v2/06.unmap_reads/demo'
    # db = '/home/wangming/data/custom_db/kraken2'
    # Kraken2(fq, out, db).install_kraken2()
    # Kraken2(fq, out, db).install_std_db()
    # Kraken2(fq, out, db).db_checker()
    # Kraken2(fq, out, db).run()
    # Kraken2(fq, out, db).report()
    # Kraken2(fq, out, db).stat()
    """
    def __init__(self, **kwargs):
        """retrieve the input fastq file, outdir director,
        check kraken2 from your PATH
        specify db_path
        """
        self = update_obj(self, kwargs, force=True)
        self.init_args()

        
    def init_args(self):
        args_init = {
            'fq': None,
            'outdir': None,
            'kraken2_db': None,
            'kraken2': None,
            'threads': 1,
            'parallel_jobs': 1,
            'save_out': False,
            'unmap_file': None,
            'trimmed': False,
            'sub_sample': 0,
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        if isinstance(self.fq, str):
            self.fq = [self.fq]
        elif isinstance(self.fq, list):
            pass
        else:
            raise ValueError('fx not valide, expect list or str, got {}'.format(
                type(self.fq).__name__))
        if not all(check_fx(self.fq)):
            raise ValueError('fx not exists: {}'.format(self.fq))
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        # clean dir
        self.clean_dir = os.path.join(self.outdir, 'clean_data')
        self.subset_dir = os.path.join(self.clean_dir, '01_subset')
        self.trim_dir = os.path.join(self.clean_dir, '02_trim_ad')
        check_path([self.subset_dir, self.trim_dir])
        # args
        self.init_kraken2()
        self.init_kraken2_db()
        
        
    def init_kraken2(self):
        # check command + db
        if not isinstance(self.kraken2, str):
            self.kraken2 = which('kraken2')
        if self.kraken2 is None:
            self.install_kraken2()
            raise ValueError('install kraken2 with above message')
        

    def install_kraken2(self):
        msg = """
        Please go throught the following instructions:
        ----------------------------------------------------------------------
        Here is brief installation instruction for Kraken2,
        See official website for more details:
        https://github.com/DerrickWood/kraken2/wiki/Manual#installation
        ----------------------------------------------------------------------

        1. Requirements

        100 GB disk space, 29 GB RAM

        2. Using conda
        $ conda install -c bioconda kraken2

        3. Using source code
        $ cd ~/biosoft/
        $ git clone https://github.com/DerrickWood/kraken2.git
        $ cd kraken2
        $ ./install_kraken2.sh <path_to_kraken2>
        
        # <path_to_kraken2> if the directory where you want to install
        # I choose "./bin"
        # copy main Kraken2 scripts to a directory in your PATH
        
        $ cp <path_to_kraken2>/kraken2{,-build,-inspect} <PATH_directory>

        $ which kraken2

        your installation is correct, if above command returns the path of kraken2.
        """
        log.info(msg)


    def init_kraken2_db(self, check_db=False):
        """
        Validate the kraken2 database
        $ kraken2-inspect --db <db_path>/standard --skip-counts
        # Database options: nucleotide db, k = 35, l = 31
        # Spaced mask = 11111111111111111111111111111111110011001100110011001100110011
        # Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
        # Total taxonomy nodes: 23303
        # Table size: 6689453365
        # Table capacity: 9554994102
        # Min clear hash value = 0

        $ kraken2-inspect --db <db_path>/minikraken2_v1_8GB --skip-counts
        # Database options: nucleotide db, k = 35, l = 31
        # Spaced mask = 11111111111111111111111111111111111111001100110011001100110011
        # Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
        # Total taxonomy nodes: 21101
        # Table size: 1399914077
        # Table capacity: 2000000000
        # Min clear hash value = 13727554816041021440

        error:
        kraken2-inspect: database ("minikraken2_v1_8GB/") does not contain necessary file taxo.k2d
        """
        # check required files: hash.k2d, opts.k2d, taxo.k2d
        f_list = [os.path.join(self.kraken2_db, i) for i in \
            ['hash.k2d', 'opts.k2d', 'taxo.k2d']]
        if not all(file_exists(f_list)):
            self.install_kraken2_db()
            raise ValueError('not a kraken2_db: {}'.format(self.kraken2_db))
        # check the details of kraken2_db
        if check_db:
            stdout = os.path.join(self.outdir, 'kraken2_inspect.stdout')
            stderr = os.path.join(self.outdir, 'kraken2_inspect.stderr')
            cmd_file = os.path.join(self.outdir, 'kraken2_inspect.cmd.sh')
            cmd = ' '.join([
                '{}'.format(which('kraken2-inspect')),
                '--threads {}'.format(self.threads),
                '--skip-counts',
                '--db {}'.format(self.kraken2_db),
                '1> {}'.format(stdout),
                '2> {}'.format(stderr),
            ])
            # save command
            with open(cmd_file, 'wt') as w:
                w.write(cmd+'\n')
            # run 
            try:
                run_shell_cmd(cmd)
            except:
                log.error('kraken2-inspect() failed')
            # check
            with open(stderr) as r:
                header = next(r)
            if header.startswith('kraken2-inspect'):
                self.install_kraken2_db()
                raise ValueError('kraken2_db not valid. see {}'.format(stderr))
        

    def install_kraken2_db(self):
        """
        Install Kraken2 database

        Option-1: Download pre-built index from Kraken2 website
        see: https://benlangmead.github.io/aws-indexes/k2
        Standard (46.8GB): https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz
        Standard-8 (7.5GB): https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20201202.tar.gz 
        Standard-16 (14.9GB): https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20201202.tar.gz 
        
        Option-2: Build index on local server
        1. Standard Database
        kraken2-build --standard --db <db_name>
        kraken2-build --standard --threads 16 --db <db_name>

        2. Custom Database
            2.1 install a taxonomy
            $ kraken2-build --download-taxonomy --db $DBNAME
        """
        msg = """
        Install standard Kraken2 database:
        -------------------------------------------------------------------------------
        Standard database is fine for most research
        (Kraken2 also support custom database, see official website for details)

        $ kraken2-build --standard --threads 16 --db <db_name>

        The program requires ~140 GB of disk space during creation, and hours to download
        NCBI datasets depends on your network.

        !!!                                                                !!!
        !!! The program stopped while processing the third library         !!!
        !!! you need to type "y"                                           !!!
        !!! to confirm the file replacement, x to "assembly_summary.txt"   !!!
        !!!                                                                !!!
        mv: replace 'assembly_summary.txt', overriding mode 0444 (r--r--r--)? y

        The full log looks like this:

        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 297 projects (449 sequences, 755.01 Mbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Masking low-complexity regions of downloaded library... done.
        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 15151 projects (32846 sequences, 60.75 Gbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Masking low-complexity regions of downloaded library... done.
        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 8609 projects (11060 sequences, 278.49 Mbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Masking low-complexity regions of downloaded library... done.
        mv: replace 'assembly_summary.txt', overriding mode 0444 (r--r--r--)? y
        Step 1/2: Performing rsync file transfer of requested files
        Rsync file transfer complete.
        Step 2/2: Assigning taxonomic IDs to sequences
        Processed 1 project (594 sequences, 3.26 Gbp)... done.
        All files processed, cleaning up extra sequence files... done, library complete.
        Downloading UniVec_Core data from server... done.
        Adding taxonomy ID of 28384 to all sequences... done.
        Masking low-complexity regions of downloaded library... done.
        Creating sequence ID to taxonomy ID map (step 1)...
        Sequence ID to taxonomy ID map complete. [0.302s]
        Estimating required capacity (step 2)...
        Estimated hash table requirement: 38219976408 bytes
        Capacity estimation complete. [9m22.571s]
        Building database files (step 3)...
        Taxonomy parsed and converted.
        CHT created with 15 bits reserved for taxid.
        Completed processing of 48085 sequences, 65038298666 bp
        Writing data to disk...  complete.
        Database files completed. [50m55.800s]
        Database construction complete. [Total: 1h0m19.486s]
        """
        log.info(msg)


    def parse_log(self, x):
        """
        parse Kraken2 log (stderr)
        : total 
        : classified
        : unclassified
        Example: 
        Loading database information... done.
        4000000 sequences (424.30 Mbp) processed in 9.158s (26207.0 Kseq/m, 2779.91 Mbp/m).
          434506 sequences classified (10.86%)
          3565494 sequences unclassified (89.14%)
          
        output: total, hit, un-hit, hit_pct
        """
        n_total = 1
        n_hit = 0
        n_unhit = 1
        try:
            with open(x) as r:
                for line in r:
                    n = line.strip().split()[0]
                    if 'processed' in line:
                        n_total = int(n)
                    elif 'sequences classified' in line:
                        n_hit = int(n)
                    elif 'sequences unclassified' in line:
                        n_unhit = int(n)
                    else:
                        pass
        except IOError as e:
            log.error(e)
        n_hit_pct = n_hit / n_total * 100
        msg = '{} ({:.1f}%) of {} sequences were classified: {}'.format(
            n_hit, n_hit_pct, n_total, file_prefix(x))
        log.info(msg)
        return (n_total, n_hit, n_unhit, n_hit_pct)


    def parse_report(self, in_report, out_stat, topN=10, tax_level='G'):
        """
        Parameters
        ----------
        in_report: str
            The --report output of kraken2
            
        out_stat: str
            The file saving the simplified report, by specific tax level
            
        topN: int
            The top N taxon, default: 100
            
        tax_level: str
            The specific tax level to show, default: ['G'],
            options: [D, K, P, C, O, F, G, S]
        
        parse Kraken2 report (--report)
        Inspect the output of kraken2 report file
        pandas table
        
        col-1: Percentage of fragments covered by the clade rooted at this taxon
        col-2: Number of fragments covered by the clade rooted at this taxon
        col-3: Number of fragments assigned directly to this taxon
        col-4: A rank code, indicating (U)nclassified, (R)oot, (D)omain, 
               (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus,
               or (S)pecies. Taxa that are not at any of these 10 ranks have 
               a rank code that is formed by using the rank code of the closest
               ancestor rank with a number indicating the distance from that rank.
               E.g., "G2" is a rank code indicating a taxon is between genus and 
               species and the grandparent taxon is at the genus rank.
        col-5: NCBI taxonomic ID number
        col-6: Indented scientific name
        
        # choose the top species, by col-3, read directly in taxon
        """
        s = ['pct', 'reads_in_clade', 'reads_in_tax', 'code', 'taxid', 'name']
        try:
            df1 = pd.read_csv(in_report, '\t', names=s)
        except IOError as e:
            log.error(e)
            return None
        # choose Geneus (G); root, unclassified
        df_un = df1.loc[df1['code'] == 'U'] # unclassified
        df_root = df1.loc[df1['code'] == 'R'] # root
        df_tax = df1.loc[df1['code'].isin([tax_level]), ] # eg: G
        # choose top ranked taxon
        df_tax = df_tax.sort_values(['reads_in_tax'], ascending=False)
        # sub-sample
        df2 = df_tax.iloc[0:topN,]
        # combine
        df = pd.concat([df_un, df2])
        ## remove white spaces
        df['name'] = df['name'].str.strip()
        # add sample name
        df['sample'] = file_prefix(in_report) # prefix
        # get total
        n_hit = df_root.loc[:, 'reads_in_clade']
        n_unhit = df_un.loc[:, 'reads_in_clade']
        # adt pct
        df['reads_hit'] = int(n_hit) # root
        df['reads_total'] = int(n_hit) + int(n_unhit) # root + unclassified
        df['hit_pct'] = df['reads_in_tax'] / (int(n_hit) + int(n_unhit)) * 100
         # sub-sample
        df3 = df.loc[:,['sample', 'name', 'reads_in_tax', 'hit_pct', 'reads_hit', 'reads_total']]        
        # custome options
#         with pd.option_context('display.expand_frame_repr', False, 'display.max_colwidth', 10):
#             print(df3)
        # save to file
        df3.to_csv(out_stat, sep='\t', index=False)
        return df3


    def run_single_fx(self, fx):
        prefix = fx_name(fx)
        stdout = os.path.join(self.outdir, prefix+'.kraken2.stdout')
        stderr = os.path.join(self.outdir, prefix+'.kraken2.stderr')
        report = os.path.join(self.outdir, prefix+'.kraken2.report')
        output = os.path.join(self.outdir, prefix+'.kraken2.out')
        stat = os.path.join(self.outdir, prefix+'.kraken2.stat')
        cmd_file = os.path.join(self.outdir, prefix+'.kraken2.cmd.sh')
        # output files
        sub_fx = os.path.join(self.subset_dir, os.path.basename(fx))
        trimmed_fx = os.path.join(self.trim_dir, os.path.basename(fx))
        # subsample the reads
        if self.sub_sample > 0:
            Fastx(fx).sample(out=sub_fx, n=self.sub_sample)
        else:
            symlink_file(fx, sub_fx)
        # update fx, trim adapters
        if self.trimmed:
            symlink_file(sub_fx, trimmed_fx)
        else:
            args_trim = {
                'fq1': sub_fx,
                'outdir': self.trim_dir,
                'threads': self.threads,
            }
            t = TrimR1(**args_trim)
            t.run()
            symlink_file(t.clean_fq1, trimmed_fx)
        # args
        args_zip = '--gzip-compressed' if fx.endswith('.gz') else \
            '--bzip2-compressed' if fx.endswith('.bz2') else ''
        cmd = ' '.join([
            '{}'.format(which('kraken2')),
            '--db {}'.format(self.kraken2_db),
            '--threads {}'.format(self.threads),
            '--use-names --quick',
            args_zip, # gzip or bzip2
            '--report {}'.format(report),
            '--output {}'.format(output),
            '{}'.format(fx),
            '1> {}'.format(stdout),
            '2> {}'.format(stderr),
        ])
        # save command
        with open(cmd_file, 'wt') as w:
            w.write(cmd+'\n')
        # run
        if os.path.exists(output) and self.overwrite is False:
            log.info('run_kraken2() skipped, file exists')
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('run_kraken2() failed, see {}'.format(stderr))
        # wrap
        self.parse_log(stderr)
        self.parse_report(report, stat)
    
    
    def run_multiple_fx(self):        
        if self.parallel_jobs > 1 and len(self.fq) > 1:
            ## Pool() run in parallel
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fx, self.fq)
        else:
            [self.run_single_fx(fx) for fx in self.fq]
            

    def report(self):
        # the template: R package
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_kraken2_report.R')
        report_html = os.path.join(self.outdir, 'kraken2_report.html')
        log_stdout = os.path.join(self.outdir, 'report.stdout')
        log_stderr = os.path.join(self.outdir, 'report.stderr')
        cmd_file = os.path.join(self.outdir, 'report.sh')
        cmd = ' '.join([
            which('Rscript'),
            qc_reportR,
            self.outdir,
            self.outdir,
            '1> {}'.format(log_stdout),
            '2> {}'.format(log_stderr),
        ])
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(report_html) and self.overwrite is False:
            log.info('report() skipped, file exists')
        else:
            run_shell_cmd(cmd)
        if not os.path.exists(report_html):
            log.error('report() failed, check {}'.format(log_stderr))
        # finish
        log.info('output - {}'.format(report_html))

        
    def run(self):
        self.run_multiple_fx()
        self.report()


def get_args():
    parser = argparse.ArgumentParser(prog='hiseq bacteria', 
                                     description='check bacteria content of fastq files')
    parser.add_argument('-i', '--fq', nargs='+', required=True, dest='fq',
        type=str, help='fastq files')
    parser.add_argument('-o', '--outdir', default=None, dest='outdir',
        metavar='OUTPUT',  help='The directory to save results. (default:\
        current directory)')
    parser.add_argument('--db', dest='kraken2_db', required=True, 
        help='Path to Kraken2 database')
    parser.add_argument('--trimmed', action='store_true',
        help='The input file was trimmed')
    parser.add_argument('--sub-sample', dest='sub_sample', type=int, 
        default=0,
        help='Sub sample the input fastq, default: [0], total reads')
    parser.add_argument('--top-n', default=10, type=int, dest='topN',
        help='Show topN species. (default: 10)')
    parser.add_argument('--save-out', dest='save_out', action='store_true',
	help='Save alignment output')
    parser.add_argument('--unmap-file', default=None, dest='unmap_file',
        help='Save unclassified reads to file')
    parser.add_argument('--kraken2', default=None, 
        help='Path to the Kraken2 script. default, search in $PATH')
    parser.add_argument('-p', '--threads', default=16, type=int, dest='threads',
        help='Number of processors to use. (default: 16)')    
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Kraken2(**args).run()


if __name__ == '__main__':
    main()

#