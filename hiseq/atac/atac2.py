"""
ATAC-seq pipeline

Input config file:

genome    outdir    name    fq1    fq2

Output:

- raw_data 
- clean_data
- align  
- bam_files 
- bw_files
- peak 
- motif 
- report
  - reads num
  - peaks num 
  - peaks annotation 
  - qc 
 
####

compare replicates:

- bam_files
- bw_files
- peak
- motif
- report
"""

import os
import re
from hiseq.utils.helper import *
from hiseq.qc.trimmer import Trimmer
from hiseq.align.alignment import Alignment
from hiseq.peak.call_peak import Macs2
from hiseq.atac.atac_utils import *
import hiseq


class AtacConfig(object):
    """
    UPDATED: 
    directory type: single, merge, multiple
    config: init_atac, single, merge, multiple
    list all files (objects)

    required args:
    fq1, fq2, genome, outdir, ...

    options:
    """
    def __init__(self, **kwargs):
        self.args = kwargs
        self.atac_type = self.mission()

        # for Atac checker
        create_dirs = kwargs.get('create_dirs', True)

        if self.atac_type == 'single':
            # self.atac_type = self.mission()
            if os.path.exists(self.args.get('design', None)):
                # for --fq1, --fq2, ...
                args_in = Json(self.args['design']).dict
                self.args.update(args_in)
            self.init_atac_single(create_dirs)
        elif self.atac_type == 'merge':
            self.init_atac_merge(create_dirs)
        elif self.atac_type == 'multiple':
            self.init_atac_multiple(create_dirs)
        else:
            pass


    def mission(self):
        """
        Determine the purpose of the ATAC, 
        1. single + merge + multiple 
        2. single
        3. merge
        4. mulltiple
        """
        args = self.args.copy()

        config = args.get('config', None)
        design = args.get('design', None)
        fq1 = args.get('fq1', None)
        fq2 = args.get('fq2', None)
        genome = args.get('genome', None)
        outdir = args.get('outdir', None)
        rep_list = args.get('rep_list', None)
        input_dir = args.get('input_dir', None)

        # atac single
        if not design is None:
            print('!BBB1')
            args_in = Json(args['design']).dict
            self.args.update(args_in) # update global self. args
            flag = 'single'
        elif all([not i is None for i in [fq1, fq2, genome, outdir]]):
            print('!BBB2')
            flag = 'multiple'
        # atac merge
        elif isinstance(rep_list, list):
            print('!BBB3')
            flag = 'merge'
        # atac multiple
        elif isinstance(input_dir, str):
            flag = 'multiple'
        # atac, parse all
        else:
            raise Exception("""unknown ATAC() arguments;
                ATAC single: design=, or fq1, fq2, genome, outdir;
                ATAC merge: rep_list=;
                ATAC multiple: input_dir=;
                """)
        return flag


    def init_atac_single(self, create_dirs=True):
        """
        initiate the config, directories, files for atac single
        update self.
        fq1, fq2, genome, outdir: args, 
        or:
        design
        """
        args = self.args.copy() # global

        self.fq1 = args['fq1']
        self.fq2 = args['fq2']
        self.genome = args['genome']
        self.outdir = args['outdir']

        ## absolute path
        self.fq1 = os.path.abspath(self.fq1)
        self.fq2 = os.path.abspath(self.fq2)
        self.outdir = os.path.abspath(self.outdir)

        ## sample name
        args['smp_name'] = args.get('smp_name', None)
        fqname = file_prefix(args['fq1'])[0]
        fqname = re.sub('[._][rR]?1$', '', fqname)
        fqname = re.sub('_1$', '', fqname)
        if args['smp_name']:
            fqname = args['smp_name']
        self.fqname = fqname

        ## outdir
        self.configdir = os.path.join(self.outdir, 'config')
        self.rawdir = os.path.join(self.outdir, 'raw_data')
        self.cleandir = os.path.join(self.outdir, 'clean_data')
        self.aligndir = os.path.join(self.outdir, 'align')
        self.bamdir = os.path.join(self.outdir, 'bam_files')
        self.bwdir = os.path.join(self.outdir, 'bw_files')
        self.peakdir = os.path.join(self.outdir, 'peak')
        self.motifdir = os.path.join(self.outdir, 'motif')
        self.qcdir = os.path.join(self.outdir, 'qc')
        self.reportdir = os.path.join(self.outdir, 'report')
        self.out_prefix = os.path.join(self.outdir, fqname)

        # self.raw_fq_list = [os.path.join(self.rawdir, i) for i in fqnames]
        # self.clean_fq_list = [os.path.join(self.cleandir, i) for i in fqnames]
        ## fastq files
        ## consider input:
        ## input: fastq, fq, fastq.gz, fq.gz
        ## raw = input
        ## clean = *.fq.gz # gzip if required.
        # fqnames = list(map(os.path.basename, [self.fq1, self.fq2]))
        self.raw_fq_list = [
            os.path.join(self.rawdir, os.path.basename(self.fq1)),
            os.path.join(self.rawdir, os.path.basename(self.fq2))]
        ## clean = fq.gz
        self.clean_fq_list = [
            os.path.join(self.cleandir, file_prefix(self.fq1)[0] + '.fq.gz'),
            os.path.join(self.cleandir, file_prefix(self.fq2)[0] + '.fq.gz')]
        ## 
        self.trim_stat = os.path.join(self.cleandir, fqname + '.qc.stat')
        self.bam_raw = os.path.join(self.aligndir, fqname, '2.*', fqname + '.bam')
        self.align_stat = os.path.join(self.aligndir, fqname + '.align.txt')
        self.bam_rmdup = os.path.join(self.bamdir, fqname + '.rmdup.bam')
        self.bam_proper_pair = os.path.join(self.bamdir, fqname + '.proper_pair.bam')
        self.peak = os.path.join(self.peakdir, fqname + '_peaks.narrowPeak')
        self.bw = os.path.join(self.bwdir, fqname + '.bigWig')
        ## qc files
        self.lendist_txt = os.path.join(self.qcdir, 'length_distribution.txt')
        self.lendist_pdf = os.path.join(self.qcdir, 'length_distribution.pdf')
        self.frip_txt = os.path.join(self.qcdir, 'FRiP.txt')

        ## update args
        args['overwrite'] = args.get('overwrite', False)
        args['threads'] = args.get('threads', 8)
        args['fqname'] = self.fqname

        ## trimming, cutadapt
        adapter3 = hiseq.utils.args.Adapter('nextera').adapters[0] #
        args['len_min']  = args.get('len_min', 20)
        args['adapter3'] = args.get('adapter3', adapter3) # nextera
        args['AD3'] = args.get('AD3', adapter3) # nextera

        ## alignment
        args['aligner'] = args.get('aligner', 'bowtie2') # bowtie alignment
        args['n_map'] = args.get('n_map', 2)
        args['align_to_chrM'] = True
        args['extra_para'] = '-X 2000'

        ## update args
        args['rawdir'] = self.rawdir
        args['cleandir'] = self.cleandir
        args['aligndir'] = self.aligndir
        args['bamdir'] = self.bamdir
        args['bwdir'] = self.bwdir
        args['peakdir'] = self.peakdir
        args['motifdir'] = self.motifdir
        args['reportdir'] = self.reportdir
        args['out_prefix'] = self.out_prefix 
        args['raw_fq_list'] = self.raw_fq_list
        args['clean_fq_list'] = self.clean_fq_list
        args['trim_stat'] = self.trim_stat
        args['bam_raw'] = self.bam_raw
        args['align_stat'] = self.align_stat
        args['bam_rmdup'] = self.bam_rmdup
        args['bam_proper_pair'] = self.bam_proper_pair
        args['peak'] = self.peak
        args['bw'] = self.bw

        # update, global
        self.args = args

        ## create directories
        if create_dirs is True:
            check_path([
                self.configdir,
                self.rawdir, 
                self.cleandir, 
                self.aligndir, 
                self.bamdir, 
                self.bwdir, 
                self.peakdir, 
                self.motifdir, 
                self.qcdir, 
                self.reportdir])


    def init_atac_merge(self, create_dirs=True):
        """
        initiate the config, directories, files, for atac merge
        update self.
        """
        args = self.args.copy() # global
        assert 'rep_list' in args # merge

        if len(args['rep_list']) < 2:
            raise Exception('require at least 2 replicats, get {} files'.format(
                len(args['rep_list'])))

        for rep in args['rep_list']:
            assert self.is_atac_single(rep) # should be ATAC single dir
        
        ## sample name
        self.smpname = args.get('smpname', None)
        self.outdir = args.get('outdir', None)
        
        if self.smpname is None:
            self.smpname = merge_names(args['rep_list'])

        if self.outdir is None:
            self.outdir = os.path.join(
                os.path.dirname(args['rep_list'][0]), 
                self.smpname)
        
        ## absolute path
        self.outdir = os.path.abspath(self.outdir)
        self.rep_list = [os.path.abspath(i) for i in args['rep_list']]

        ## outdir
        self.configdir = os.path.join(self.outdir, 'config')
        self.aligndir = os.path.join(self.outdir, 'align')
        self.bamdir = os.path.join(self.outdir, 'bam_files')
        self.bwdir = os.path.join(self.outdir, 'bw_files')
        self.peakdir = os.path.join(self.outdir, 'peak')
        self.motifdir = os.path.join(self.outdir, 'motif')
        self.qcdir = os.path.join(self.outdir, 'qc')
        self.reportdir = os.path.join(self.outdir, 'report')

        ## names
        self.bam = os.path.join(self.bamdir, self.smpname + '.bam')
        self.peak = os.path.join(self.peakdir, self.smpname + '_peaks.narrowPeak')
        self.bw = os.path.join(self.bwdir, self.smpname + '.bigWig')
        self.lendist_txt = os.path.join(self.qcdir, 'length_distribution.txt')
        self.lendist_pdf = os.path.join(self.qcdir, 'length_distribution.pdf')
        self.frip_txt = os.path.join(self.qcdir, 'FRiP.txt')
        self.cor_npz = os.path.join(self.qcdir, 'cor.bam.npz')
        self.cor_counts = os.path.join(self.qcdir, 'cor.bam.counts.tab')
        self.idr_txt = os.path.join(self.qcdir, 'idr.txt')
        self.idr_log = os.path.join(self.qcdir, 'idr.log')
        self.idr_png = self.idr_txt + '.png'
        self.peak_overlap = os.path.join(self.qcdir, 'peak_overlap.pdf')
     
        if create_dirs:
            check_path([
                self.configdir,
                self.aligndir,
                self.bamdir,
                self.bwdir,
                self.peakdir,
                self.motifdir,
                self.qcdir,
                self.reportdir])


    def init_atac_multiple(self, create_dirs=True):
        """
        initiate the config, directories, files, for atac multiple samples
        **always for summary report for a project**
        update self.
        """
        args = self.args.copy() # global
        assert in_dict(args, 'input_dir') # merge

        # outdir, input_dir + summary
        self.outdir = args.get('outdir', None)
        if self.outdir is None:
            self.outdir = os.path.join(args['input_dir'], 'summary')

        ## outdir
        # print('!xxxx', self.outdir)
        # sys.exit('!mmmm')
        self.configdir = os.path.join(self.outdir, 'config')
        self.reportdir = os.path.join(self.outdir, 'report')

        if create_dirs:
            check_path([
                self.configdir,
                self.reportdir])

        # list all atac_single, atac_merge directories
        dirs = listfiles(args['input_dir'], include_dir=True)
        dirs = [os.path.abspath(i) for i in dirs if os.path.isdir(i)]
        dirs = sorted(dirs)

        # rep, merge
        self.single_list = []
        self.merge_list = []
        for x in dirs:
            if self.is_atac_single(x):
                print('single', x)
                self.single_list.append(x)
            elif self.is_atac_merge(x):
                print('merge', x)
                self.merge_list.append(x)
            elif self.is_atac_multiple(x):
                print('multiple', x)
            else:
                print('unknown', x)
                continue


    def is_atac_single(self, x):
        """
        The directory is a atac single
        check required attribtes:
        ...
        """
        x_design = x + '.config.json'
        x_pickle = x + '.config.pickle'
        if check_file([x_design, x_pickle]):
            chk = AtacConfig(design=x_design, create_dirs=False)
            args = pickle_to_dict(x_pickle)

            tags = in_attr(chk, [
                'rawdir',
                'cleandir',
                'aligndir',
                'bamdir',
                'peakdir',
                'qcdir',
                'reportdir',
                'frip_txt',
                'lendist_txt',
                'align_stat',
                'peak',
                'bam_proper_pair'], True)

            log.info('Check files for ATACseq single: {}'.format(x))
            check_file(tags, show_log=True)
            return all(map(os.path.exists, tags))
        else:
            logging.error('config files not found: {}, {}'.format(
                x_design, x_pickle))
        

    def is_atac_merge(self, x):
        """
        The directory is a atac merge replicates
        check required attributes:
        ...
        """
        sample_list = os.path.join(x, 'config', 'sample_list.txt')

        if os.path.exists(sample_list):
            n = 0
            rep_list = []
            with open(sample_list) as r:
                for line in r:
                    n += 1
                    rep_list.append(line.rstrip())
                    if n > 100:
                        break
            # further config
            chk = AtacConfig(rep_list=rep_list, create_dirs=False)
            # print(dir(chk))

            tags = in_attr(chk, [
                'configdir',
                'aligndir',
                'bamdir',
                'peakdir',
                'qcdir',
                'reportdir',            
                'bam',
                'peak',
                'frip_txt',
                'cor_counts',
                'idr_txt',
                'peak_overlap'], True)
            log.info('Check files for ATACseq merge: {}'.format(x))
            check_file(tags, show_log=True)
            return all(map(os.path.exists, tags))


    def is_atac_multiple(self, x):
        """
        The directory is a atac multiple samples report directories
        check required attributes:
        ...
        """        
        single_list_file = os.path.join(x, 'config', 'sample_list.txt')
        merge_list_file = os.path.join(x, 'config', 'merge_list.txt')

        ## file exists
        if not check_file([single_list_file, merge_list_file]):
            return False

        # for single
        n = 0
        single_list = []
        with open(single_list_file) as r:
            for line in r:
                n += 1
                single_list.append(line.rstrip())
                if n > 100:
                    break
        chk1 = all([self.is_atac_single(i) for i in single_list])

        # for merge
        n = 0
        merge_list = []
        with open(merge_list_file) as r:
            for line in r:
                n += 1
                merge_list.append(line.rstrip())
                if n > 100:
                    break
        chk2 = all([self.is_atac_merge(i) for i in merge_list])
        log.info('Check files for ATACseq multiple: {}'.format(x))

        return chk1 and chk2


class AtacSingle(object):
    """
    ATAC-seq analysis pipeline
    replicates: 1,2,3,...
    paired-end:

    require: 
    design
    or
    fq1, fq2, genome, outdir
    """
    def __init__(self, **kwargs):
        """
        Parsing design file: Json
        fq1, fq2, genome, outdir
        """
        self.args = kwargs
        self.status = self.init_atac() # update all variables: *.config, *.args
        

    def init_atac(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        # args = self.args.copy()
        self.config = AtacConfig(**self.args) # update, global
        self.args.update(self.config.args) # update, global
        assert in_attr(self.config, ['fq1', 'fq2', 'genome', 'outdir'])
        assert self.config.atac_type == 'single'

        # required
        self.fq1 = os.path.abspath(self.args.get('fq1', None))
        self.fq2 = os.path.abspath(self.args.get('fq2', None))
        self.genome = self.args.get('genome', None)
        self.outdir = os.path.abspath(self.args.get('outdir', None))
        self.args['fq1'] = self.fq1
        self.args['fq2'] = self.fq2
        self.args['outdir'] = self.outdir

        # save config to files
        args_pickle = os.path.join(self.config.configdir, 'arguments.pickle')
        args_json = os.path.join(self.config.configdir, 'arguments.json')
        args_txt = os.path.join(self.config.configdir, 'arguments.txt')
        # check arguments
        chk1 = args_checker(self.args, args_pickle)
        Json(self.args).writer(args_json)
        args_logger(self.args, args_txt)
        self.args['overwrite'] = self.args.get('overwrite', False)
        chk2 = self.args['overwrite'] is False

        # status
        return all([chk1, chk2])


    #######################################
    ## main pipeline
    def prep_raw(self, copy=False):
        """
        Copy raw data to dest dir
        if not: create a symlink
        
        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.config.raw_fq_list

        # copy
        if copy is True:
            shutil.copy(self.fq1, raw_fq1)
            shutil.copy(self.fq2, raw_fq2)
        else:
            symlink(self.fq1, raw_fq1, absolute_path=True)
            symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=False):
        """
        Trim reads:
        hiseq.qc.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        # args = self.args.copy()
        fq1, fq2 = self.config.raw_fq_list
        clean_fq1, clean_fq2 = self.config.clean_fq_list

        # update args
        args = self.config.args # all
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.cleandir
        
        if trimmed is True:
            print('!xxxx')
            # create symlink from rawdir
            # if raw is not gzipped, do it
            if is_gz(fq1) and is_gz(fq2):
                symlink(fq1, clean_fq1)
                symlink(fq2, clean_fq2)
            else:
                # gzip files
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
                gzip_cmd(fq2, clean_fq2, decompress=False, rm=False)
            # for s, d in zip(self.config.raw_fq_list, self.config.clean_fq_list):
            #     self.symlink(s, d)
        else:
            print('!yyyy')
            if check_file(self.config.clean_fq_list):
                log.info('trim() skipped, file exists: {}'.format(
                    self.config.clean_fq_list))
            else:
                Trimmer(**args).run()


    def align(self):
        """
        Alignment PE reads to reference genome, using bowtie2

        -X 2000, ...
        """
        fq1, fq2 = self.config.clean_fq_list

        # update arguments
        args = self.config.args
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.aligndir
        args['extra_para'] = '-X 2000' # repeat, 
        args['aligner'] = 'bowtie2' # repeat,

        if check_file(self.config.bam_rmdup):
            log.info('align() skipped, file exists: {}'.format(
                self.config.bam_rmdup))
        else:
            Alignment(**args).run()


    def get_raw_bam(self):
        """
        Get the align bam file
        from bamdir/1., 2., ...
        !!! specific: 2.genome/*.bam
        """
        # bamdir = os.path.join(self.config.aligndir, '2.*')
        bamdir = self.config.align_stat.rstrip('.align.txt')
        bamlist = listfiles2('*.bam', bamdir, recursive=True)
        # [chrM, genome]
        return(bamlist[1]) # genome


    def bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        bam_raw = self.get_raw_bam()

        # rmdup
        if rmdup is True:
            Bam(bam_raw).rmdup(self.config.bam_rmdup)
        else:
            shutil.copy(bam_raw, self.config.bam_rmdup)

        if check_file(self.config.bam_proper_pair):
            log.info('bam_rmdup() skipped, file exists: {}'.format(
                self.config.bam_proper_pair))
        else:
            # save proper pair reads
            Bam(self.config.bam_rmdup).proper_pair(self.config.bam_proper_pair)
            # index
            Bam(self.config.bam_proper_pair).index()


    def callpeak(self):
        """
        Call peaks using MACS2
        """
        args = self.config.args.copy()
        bam = self.config.bam_proper_pair
        genome = args.pop('genome', None)
        output = args.pop('peakdir', None)
        prefix = args.pop('fqname', None)

        if check_file(self.config.peak):
            log.info('callpeak() skipped, file exists: {}'.format(
                self.config.peak))
        else:
            Macs2(bam, genome, output, prefix, atac=True, **args).callpeak()


    def qc_lendist(self):
        """
        Create length distribution plot
        """
        # create plot
        pkg_dir = os.path.dirname(hiseq.__file__)
        lendistR = os.path.join(pkg_dir, 'bin', 'atac_qc_lendist.R')
        cmd = 'Rscript {} {} {}'.format(
            lendistR,
            self.config.lendist_txt,
            self.config.lendist_pdf)

        if check_file(self.config.lendist_txt):
            log.info('qc_lendist() skipped, file exists: {}'.format(
                self.config.lendist_txt))
        else:
            # compute length distribution
            frag_length(self.config.bam_proper_pair, self.config.lendist_txt)
            # create plot
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('qc_lendist() failed.')


    def qc_frip(self):
        """
        Compute frip
        """
        if check_file(self.config.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.config.frip_txt))
        else:
            frip, n, total = cal_FRiP(self.config.peak, 
                self.config.bam_proper_pair)

            hd = ['FRiP', "peak_reads", "total_reads", "id"]
            n = list(map(str, [frip, n, total]))
            n.append('self.config.fqname')
            with open(self.config.frip_txt, 'wt') as w:
                w.write('\t'.join(hd) + '\n')
                w.write('\t'.join(n) + '\n')



    def qc_mito(self):
        """
        Mito reads, percentage
        """
        pass


    def qc_tss(self):
        """
        TSS enrichment
        require:
        BAM, peak, TSS, ...
        """
        pass


    def qc(self):
        """
        QC for ATACseq sample
        FRiP
        Fragment size (length distribution)
        TSS
        """
        self.qc_frip()
        self.qc_lendist()
        self.qc_mito()
        self.qc_tss()


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir    = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atac_report_single.R')
        atac_report_html = os.path.join(
            self.config.reportdir, 
            'atac_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.config.outdir,
            self.config.reportdir)
        
        if check_file(atac_report_html):
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('report() failed.')


    def run(self):
        """
        Run all steps for ATACseq pipeline
        """
        # init dir
        args = self.config.args.copy()

        trimmed = args.get('trimmed', False)
        copy_raw_fq = args.get('copy_raw_fq', False)
        rmdup = args.get('rmdup', True) # remove PCR dup

        # 1. copy raw data
        self.prep_raw(copy_raw_fq)

        # 2. trim
        self.trim(trimmed)

        # 3. align
        self.align()

        # 4. rmdup
        self.bam_rmdup(rmdup)

        # 5. peak 
        self.callpeak()
        
        # 6. motif
        # self.motif()

        # 7. qc
        self.qc()

        # 8.report
        self.report()


class AtacMerge(object):
    """
    Merge replicates for ATAC samples:
    require: rep_list
    """
    def __init__(self, **kwargs):
        self.args = kwargs
        assert in_dict(self.args, ['rep_list'])
        self.status = self.init_atac() # update all variables: *.config, *.args


    def init_atac(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        self.config = AtacConfig(**self.args) # update, global
        self.args.update(self.config.args) # update, global
        # assert in_attr(self., ['rep_list'])        
        assert self.config.atac_type == 'merge'
        self.args['rep_list'] = [os.path.abspath(i) for i in self.args['rep_list']]

        # config for all replicates
        self.rep_config = []
        # print(self.args['rep_list'])
        for i in self.args['rep_list']:
            i_design = os.path.join(i, 'config', 'arguments.json')
            self.rep_config.append(AtacConfig(design=i_design))

        # global variables, replicates
        self.align_txt = [c.align_stat for c in self.rep_config]
        self.bam_list = [c.bam_proper_pair for c in self.rep_config]
        self.peak_list = [c.peak for c in self.rep_config]
        self.frip_list = [c.frip_txt for c in self.rep_config]

        # save config to files
        args_pickle = os.path.join(self.config.configdir, 'arguments.pickle')
        args_json = os.path.join(self.config.configdir, 'arguments.json')
        args_txt = os.path.join(self.config.configdir, 'arguments.txt')
        replist_txt = os.path.join(self.config.configdir, 'sample_list.txt')
        # check arguments
        chk1 = args_checker(self.args, args_pickle)
        Json(self.args).writer(args_json)
        args_logger(self.args, args_txt)
        self.args['overwrite'] = self.args.get('overwrite', False)
        chk2 = self.args['overwrite'] is False
        # save replist
        with open(replist_txt, 'wt') as w:
            w.write('\n'.join(self.args['rep_list']) + '\n')

        # status
        return all([chk1, chk2])


    def mergebam(self):
        """
        Merge replicates, BAM
        """
        cmd = 'samtools merge {} {} && \
            samtools sort -o {} {} && \
            samtools index {}'.format(
            self.config.bam + '.tmp',
            ' '.join(self.bam_list),
            self.config.bam,
            self.config.bam + '.tmp',
            self.config.bam)

        if os.path.exists(self.config.bam):
            log.info('mergebam() skipped, file exists: {}'.format(
                self.config.bam))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('mergebam() failed.')


    def callpeak(self):
        """
        Call peaks using MACS2
        """        
        bam = self.config.bam
        genome = self.rep_config[0].genome
        output = self.config.peakdir
        prefix = self.config.smpname
        if os.path.exists(self.config.peak):
            log.info('callpeak() skipped, file exists: {}'.format(
                self.config.peak))
        else:
            m = Macs2(bam, genome, output, prefix, atac=True).callpeak()


    def bam_cor(self, window=500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp
        
        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
            --outRawCounts *counts.tab -b bam
        """
        multiBamSummary = shutil.which('multiBamSummary')

        cmd = '{} bins --binSize {} --smartLabels -o {} --outRawCounts {} \
            -b {}'.format(
                multiBamSummary, 
                window,
                self.config.cor_npz,
                self.config.cor_counts,
                ' '.join(self.bam_list))

        if os.path.exists(self.config.cor_counts):
            log.info('bam_cor() skipped, file.exsits: {}'.format(
                self.config.cor_counts))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('bam_cor() failed.')


    def peak_overlap(self):
        """
        Compute the overlaps between overlaps
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        peak_overlapR = os.path.join(pkg_dir, 'bin', 'atac_peak_overlap.R')

        cmd = 'Rscript {} {} {}'.format(
            peak_overlapR,
            self.config.qcdir,
            ' '.join(self.peak_list))

        if os.path.exists(self.config.peak_overlap):
            log.info('peak_overlap() skipped, file exists: {}'.format(
                self.config.peak_overlap))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('peak_overlap() failed.')


    def rep_idr(self):
        """
        Calculate IDR for replicates
        1 vs 1
        peak files
        """
        idr = shutil.which('idr') # command
        # sort peakfile by p.value
        cmd = 'sort -k8,8nr -o {} {} && \
            sort -k8,8nr -o {} {} && \
            {} --input-file-type narrowPeak --rank p.value --plot \
            --output-file {} --log-output-file {} --samples {}'.format(
            self.peak_list[0],
            self.peak_list[0],
            self.peak_list[1],
            self.peak_list[1],
            idr,
            self.config.idr_txt,
            self.config.idr_log,
            ' '.join(self.peak_list))

        if os.path.exists(self.config.idr_txt):
            logging.info('rep_idr() skipped, file exists: {}'.format(
                self.config.idr_txt))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('rep_idr() failed.')


    def get_frip(self):
        """
        Save all FRiP.txt file to one
        """
        if not os.path.exists(self.config.frip_txt):
            with open(self.config.frip_txt, 'wt') as w:
                for f in self.frip_list:
                    with open(f) as r:
                        for line in r:
                            w.write(line)


    def get_align_txt(self):
        """
        Copy align.txt files to align/
        """
        pass


    def tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        pass


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir    = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atac_report.R')
        atac_report_html = os.path.join(
            self.config.reportdir, 
            'atac_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.config.outdir,
            self.config.reportdir)
        
        if check_file(atac_report_html):
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('report() failed.')


    def run(self):
        ## question:
        ## merge bam or merge peaks?!
        # merge
        self.mergebam()
        self.callpeak()
        # qc
        self.bam_cor()
        self.peak_overlap()
        self.rep_idr()
        self.get_frip()
        # self.get_align_txt()
        # self.tss_enrich()
        self.report()


class AtacMultiple(object):
    """
    Organize multiple ATAC-seq directories
    create a summary report
    """
    def __init__(self, **kwargs):
        self.args = kwargs
        assert in_dict(self.args, ['input_dir'])
        self.status = self.init_atac() # update all variables: *.config, *.args


    def init_atac(self):
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        self.config = AtacConfig(**self.args) # update, global
        self.args.update(self.config.args) # update, global
        assert self.config.atac_type == 'multiple'
        self.args['input_dir'] = os.path.abspath(self.args['input_dir'])

        # save config to files
        args_pickle = os.path.join(self.config.configdir, 'arguments.pickle')
        args_json = os.path.join(self.config.configdir, 'arguments.json')
        args_txt = os.path.join(self.config.configdir, 'arguments.txt')
        single_list = os.path.join(self.config.configdir, 'single_list.txt')
        merge_list = os.path.join(self.config.configdir, 'merge_list.txt')

        # check arguments
        chk1 = args_checker(self.args, args_pickle)
        Json(self.args).writer(args_json)
        args_logger(self.args, args_txt)
        self.args['overwrite'] = self.args.get('overwrite', False)
        chk2 = self.args['overwrite'] is False
        
        # save replist
        with open(single_list, 'wt') as w:
            w.write('\n'.join(self.config.single_list) + '\n')

        # save mergelist
        with open(merge_list, 'wt') as w:
            w.write('\n'.join(self.config.merge_list) + '\n')

        # status
        return all([chk1, chk2])


    def report(self):
        """
        Create report for multiple samples within the directories
        html
        """
        pkg_dir    = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atac_report.R')
        atac_report_html = os.path.join(
            self.config.reportdir, 
            'atac_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            os.path.dirname(self.config.outdir),
            self.config.reportdir)
        
        if check_file(atac_report_html):
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('report() failed.')


    def run(self):
        """
        Run AtacMultiple
        """
        self.report()


class Atac(object):
    """
    Main port for ATAC pipeline
    :
    design : one file
    config : multiple fq files, one fq per line 
    --fq1, --fq2 : multiple fq files, 
    ...
    """
    def __init__(self, **kwargs):
        """
        options: 
        1. config: *.txt
        2. design: *.json
        3. --fq1, --fq2, --genome, --outdir: 
        ...
        """
        self.args = kwargs
        self.init_atac()


    def init_atac(self):
        args = self.args.copy()

        config = args.get('config', None)
        design = args.get('design', None)
        fq1 = args.get('fq1', None)
        fq2 = args.get('fq2', None)
        genome = args.get('genome', None)
        outdir = args.get('outdir', None)
        
        assert is_path(outdir)

        if check_file(config, show_log=False):
            # update config to list
            if isinstance(config, str):
                self.args['config'] = [config]
            print('!aaaa')
            self.json_list = self.config_to_json()
        elif check_file(design, show_log=False):
            # udpate design to list
            if isinstance(design, str):
                self.args['design'] = [design]
            print('!bbbb')
            self.json_list = self.design_to_json()
        elif all([not i is None for i in [fq1, fq2, genome, outdir]]):
            print('!cccc')
            self.json_list = self.args_to_json()
        else:
            raise Exception('Check arguments: config=, design=, or \
                fq1=, fq2=, genome=, outdir=')


    def config_to_json(self):
        """
        arguments: --config
        convert config to multiple *.json files
        """
        print('!AAAA1')
        args = self.args.copy()
        assert in_dict(args, 'config')

        # update arguments
        json_list = []
        nline = 0
        for config in args['config']:
            with open(config) as r:
                for line in r:
                    nline += 1
                    if line.startswith('#') or line.strip() == '':
                        continue # comment, blank rows
                    fields = line.strip().split('\t')
                    if len(fields) < 4:
                        log.warning('line-{}, 4 fields required, get {}'.format(
                            nline, len(fields)))
                        continue
                    # required args
                    fq1, fq2, genome, outdir = fields[:4]
                    assert check_path(outdir)
                    # prefix
                    fqname = file_prefix(fq1)[0] # name
                    fqname = re.sub('_1$', '', fqname) #

                    # sample dir
                    suboutdir = os.path.join(outdir, fqname)
                    config_json = suboutdir + '.config.json' # file
                    config_pickle = suboutdir + '.config.pickle' # file
                    d = {
                        'fq1': fq1, 
                        'fq2': fq2, 
                        'genome': genome, 
                        'outdir': suboutdir}
                    # d.update(args) # global config
                    Json(d).writer(config_json)
                    dict_to_pickle(d, config_pickle)
                    json_list.append(config_json)
                    print('!AAAA2', config_json)

        self.json_list = json_list
        return json_list


    def design_to_json(self):
        """
        Save design to files in outdir
        for: single fastq file
        """
        print('!BBBB1')
        args = self.args.copy()
        assert isinstance(args['design'], list)

        ## save json to outdir
        json_list = []
        for i in args['design']:
            d = Json(i).dict
            check_path(d['outdir']) # output dir
            fq1 = d.get('fq1', None)
            fq2 = d.get('fq2', None)
            genome = d.get('genome', None)
            outdir = d.get('outdir', None)
            assert isinstance(fq1, str)
            assert isinstance(fq2, str)
            assert isinstance(genome, str)
            assert isinstance(outdir, str)

            fqname = file_prefix(d['fq1'])[0]
            fqname = re.sub('_1$', '', fqname)            

            # sample dir
            # suboutdir = os.path.join(outdir, fqname)
            config_json = outdir + '.config.json'
            config_pickle = outdir + '.config.pickle'
            d.update(args) # global config

            Json(d).writer(config_json)
            dict_to_pickle(d, config_pickle)
            json_list.append(config_json)
            print('!BBBB2', config_json)

        self.json_list = json_list
        return json_list


    def args_to_json(self):
        """
        fq1, fq2, list
        genome, str
        outdir, str
        ...
        """
        print('!CCCC1')
        args = self.args.copy()
        assert in_dict(args, ['fq1', 'fq2', 'genome', 'outdir'])
        # assert isinstance(args['fq1'], list)
        # assert isinstance(args['fq2'], list)
        if isinstance(args['fq1'], str):
            args['fq1'] = [args['fq1']]
        if isinstance(args['fq2'], str):
            args['fq2'] = [args['fq2']]
        assert isinstance(args['genome'], str)
        assert isinstance(args['outdir'], str)
        if not len(self.args['fq1']) == len(self.args['fq2']):
            raise Exception('fq1 and fq2 are not in same length')

        ## convert to design (Json file)
        json_list = []
        for fq1, fq2 in zip(args['fq1'], args['fq2']):
            fqname = file_prefix(fq1)[0]
            fqname = re.sub('_1$', '', fqname)

            # sample dir
            suboutdir = os.path.join(args['outdir'], fqname)
            config_json = suboutdir + '.config.json'
            config_pickle = suboutdir + '.config.pickle'
            d = {
                'fq1': fq1, 
                'fq2': fq2, 
                'genome': args['genome'], 
                'outdir': suboutdir}
            # d.update(args) # global config
            args.update(d)
            Json(args).writer(config_json)
            dict_to_pickle(args, config_pickle)
            json_list.append(config_json)

        self.json_list = json_list
        return json_list

        
    def unique_names(self, name_list):
        """
        Get the name of replicates
        common in left-most
        """
        name_list = [os.path.basename(i) for i in name_list]
        name_list = [re.sub('.rep[0-9].*', '', i) for i in name_list]
        return list(set(name_list))


    def get_samples(self, l, x):
        """
        get samples from list: x from l
        """
        return [i for i in l if x in os.path.basename(i)]


    def run(self):
        # for each sample:
        smp_dirs = []
        # sys.exit(self.json_list)
        self.args.pop('design', None) # 
        self.args.pop('rep_list', None) # 

        for i in self.json_list:
            s = AtacSingle(design=i, **self.args).run()
            smp_dirs.append(re.sub('\.config.json', '', i))

        # for merge sample:
        smp_names = self.unique_names(smp_dirs)
        print('!aaaa', smp_names)
        for x in smp_names:
            print('!MMMM1')
            x_list = self.get_samples(smp_dirs, x)
            if len(x_list) < 2:
                continue
            else:
                print('!MMMM2')
                AtacMerge(rep_list=x_list, **self.args).run()

        # for multiple sample:
        # summary report for a project
        # dirname of smp_dirs
        # prj_dir = os.path.dirname(smp_dirs[0])
        self.args['input_dir'] = os.path.dirname(smp_dirs[0])
        if len(self.json_list) > 1:
            AtacMultiple(**self.args).run()



