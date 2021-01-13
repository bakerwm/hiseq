
"""
Call peaks using macs2

Install macs2 for python3

"""

import os
from hiseq.utils.helper import * # all help functions

class Macs2(object):
    """
    Run macs2 for BAM files
    1. macs2 callpeak -f BAM -t {IP.bam} -c {input.bam} -g {gsize} --outdir {out_dir} 
        -n {prefix} -B --SPMR {--broad} {--keep-dup auto}
    2. macs2 bdgcmp -t {prefix}_treat_pileup.bdg -c {prefix}_control_lambda.bdg -o {prefix}.ppois.bdg -m ppois
    3. macs2 bdgcmp -t {prefix}_treat_pileup.bdg -c {prefix}_control_lambda.bdg -o {prefix}.FE.bdg -m FE
    4. macs2 bdgcmp -t {prefix}_treat_pileup.bdg -c {prefix}_control_lambda.bdg -o {prefix}.logLR.bdg -m logLR -p 0.00001
    """
    def __init__(self, ip, genome, output, prefix=None, control=None, 
        atac=False, overwrite=False, genome_size=0, gsize_file=None, **kwargs):
        """Parse the parameters
        venv, the virtualenv created for macs2, running in Python2
        """
        self.ip = ip
        self.genome = genome
        self.output = output
        self.control = control
        self.overwrite = overwrite
        self.genome_size = genome_size
        self.gsize_file = gsize_file
        self.prefix = prefix
        self.atac = atac

        if prefix is None:
            prefix = file_prefix(ip)[0]
            # prefix = os.path.splitext(os.path.basename(ip))[0]

        self.gsize = self.get_gsize()
        if self.gsize is None:
            raise ValueError('unknown genome: {}'.format(genome))

        # get gsize_file, only for genome=None
        if not isinstance(genome, str):
            if genome_size == 0 or gsize_file is None:
                raise ValueError('genome_size and gsize_file required, if genome=None')

        is_path(self.output)


    def get_gsize(self):
        """
        Genome size for macs2 -g option
        """
        d = {
          'dm6': 'dm',
          'dm3': 'dm',
          'mm9': 'mm',
          'mm10': 'mm',
          'hg19': 'hs',
          'hg38': 'hs',
          'GRCh38': 'hs'}

        # prefer: if input genome_size is int
        if isinstance(self.genome_size, int) and self.genome_size > 0:
            gsize = self.genome_size # from args
        else:
            gsize = d.get(self.genome, None) # from macs2

        return gsize


    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.out',
                                            delete=False)
        return tmpfn.name


    def callpeak(self):
        """Call peaks using MACS"""
        callpeak_log = os.path.join(self.output, self.prefix + '.macs2.callpeak.out')
        if self.atac is True:
            # ATAC-seq
            macs2_cmd = 'macs2 callpeak --nomodel --shift -100 --extsize 200 \
                -t {} -g {} --outdir {} -n {} 2> {}'.format(
                    self.ip,
                    self.gsize,
                    self.output,
                    self.prefix,
                    callpeak_log)
            peak = os.path.join(self.output, self.prefix + '_peaks.narrowPeak')
        else:
            # ChIP-seq
            macs2_cmd = 'macs2 callpeak --nomodel --extsize 150 -t {} -g {} --outdir {} -n {} \
                --keep-dup auto -B --SPMR'.format(
                    self.ip, 
                    self.gsize,
                    self.output,
                    self.prefix)
            if not self.control is None:
                macs2_cmd += ' -c {}'.format(self.control)

            macs2_cmd += ' 2> {}'.format(callpeak_log)
            # peak file
            peak = os.path.join(self.output, self.prefix + '_peaks.narrowPeak')
            # save cmd
            cmd_txt = os.path.join(self.output, 'cmd.sh')
            with open(cmd_txt, 'wt') as w:
                w.write(macs2_cmd + '\n')

       # output file
        # macs2_out = os.path.join(self.output, self.prefix + '_peaks.xls')
        try:
            if os.path.exists(peak) and self.overwrite is False:
                logging.info('file exists, skip macs2 callpeak')
            else:
                logging.info('run macs2 callpeak')
                run_shell_cmd(macs2_cmd)
        except:
            logging.error('callpeak() failed')

        return peak


    def bdgcmp(self, opt='ppois'):
        """Options for -m:
        {ppois,qpois,subtract,logFE,FE,logLR,slogLR,max}
        """
        # use output of callpeak
        ip_bdg = os.path.join(self.output, self.prefix + '_treat_pileup.bdg')
        input_bdg = os.path.join(self.output, self.prefix + '_control_lambda.bdg')

        if not os.path.exists(ip_bdg) or not os.path.exists(input_bdg):
            raise ValueError('*.bdg file not found, need to run .callpeak() first')
        if not opt in ['ppois', 'qpois', 'subtract', 'logFE', 'FE', 'logLR', 'slogLR', 'max']:
            raise ValueError('unknown option: opt={}'.format(opt))

        if opt == 'logLR':
            opt_ext = '-p 0.00001'
        else:
            opt_ext = ''

        out_bdg = os.path.join(self.output, self.prefix + '.' + opt + '.bdg')
        log = os.path.join(self.output, self.prefix + '.macs2.bdgcmp.out')
        cmd1 = 'macs2 bdgcmp -t {} -c {} -o {} -m {} {}'.format(
            ip_bdg, input_bdg, out_bdg, opt, opt_ext)
        
        if os.path.exists(out_bdg) and self.overwrite is False:
            logging.info('file exists, skip macs2 bdgcmp')
        else:
            run_shell_cmd(cmd1)

        # sort output *.bdg
        cmd2 = 'sort -k1,1 -k2,2n -o {} {}'.format(out_bdg, out_bdg)

        # cnvert *.bdg to *.bigWig
        gsize_file = self.gsize_file
        out_bw  = os.path.join(self.output, self.prefix + '.' + opt + '.bigWig')
        cmd3 = 'bedGraphToBigWig {} {} {}'.format(out_bdg, gsize_file, out_bw)

        try:
            if os.path.exists(out_bw) and self.overwrite is False:
                logging.info('file exists, skip bg2bw')
            else:
                run_shell_cmd(cmd2)
                run_shell_cmd(cmd3)
        except:
            logging.error('bdgcmp() failed')


    def bdgpeakcall(self):
        pass


    def bdgopt(self):
        pass


    def cmbreps(self):
        pass


    def bdgdiff(self):
        pass


    def filterdup(self):
        pass


    def predicted(self):
        pass


    def pileup(self):
        pass


    def randsample(self):
        pass


    def refinepeak(self):
        pass


    def annotation(self):
        """Annotate narrowpeak using HOMER annotatePeaks.pl script
        BED foramt
        """
        anno_exe = shutil.which('annotatePeaks.pl')
        if not os.path.exists(anno_exe):
            logging.error('command not exists, skip annotation - {}'.format(anno_exe))
            return None

        # # narrow peak file
        # peak_listA = listfiles2('*_peaks.narrowPeak', self.output)
        # peak_listB = listfiles2('*_peaks.broadPeak', self.output)
        peak_listA = listfille(self.output, '*_peaks.narrowPeak')
        peak_listB = listfille(self.output, '*_peaks.broadPeak')
        peak_list = peak_listA + peak_listB

        anno_list = []
        for peak_file in peak_list:
            peak_anno = peak_file + '.annotation'
            peak_log = peak_anno + '.log'
            cmd = 'perl {} {} {} 1> {} 2> {}'.format(
                anno_exe, peak_file, self.genome, peak_anno, peak_log)
            
            try:
                # check existence
                if os.path.exists(peak_anno) and self.overwrite is False:
                    log.info('file exists, skip annotation: {}'.format(peak_anno))
                else:
                    run_shell_cmd(cmd)

                anno_list.append(peak_anno)
            except:
                logging.error('annotation() failed')

        return anno_list


    def broadpeak_annotation(self):
        """Annotate broadpeak using HOMER annotatePeaks.pl script
        BED foramt
        """
        anno_exe = shutil.which('annotatePeaks.pl')
        if not os.path.exists(anno_exe):
            logging.error('command not exists, skip annotation - {}'.format(anno_exe))
            return None

        # broad peak file
        broadpeak = os.path.join(self.output, self.prefix + '_peaks.broadPeak')
        if not os.path.exists(broadpeak):
            raise ValueError('require .callpeak() first - {}'.format(broadpeak))

        # run
        anno_peak = os.path.join(self.output, self.prefix + '_peaks.broadPeak.annotation')
        anno_log = os.path.join(self.output, self.prefix + '_peaks.broadPeak.annotation.log')
        cmd = 'perl {} {} {} 1> {} 2> {}'.format(
            anno_exe, broadpeak, self.genome, anno_peak, anno_log)

        try:
            run_shell_cmd(cmd)
        except:
            logging.error('broadpeak_annotation() failed')
        
        return anno_peak


    def get_effect_size(self):
        """
        Extract the effective depth of macs2 files
        parse the file: output/*_peaks.xls
        tags after filtering in treatment
        tags in treatment
        tags after filtering in control
        tags in control
            tag size is determined as 100 bps
            total tags in treatment: 7978071
            tags after filtering in treatment: 2384854
            maximum duplicate tags at the same position in treatment = 1
            Redundant rate in treatment: 0.70
            total tags in control: 10555283
            tags after filtering in control: 6639591
            maximum duplicate tags at the same position in control = 1
            Redundant rate in control: 0.37
            d = 122
            alternative fragment length(s) may be 122 bps
        """
        # search the xls file
        f = listfile(self.output, '*peaks.xls')
        if not os.path.exists(f[0]):
            raise ValueError('file missing in macs2 callpeak output: {}'.format(f))

        # top
        topN = 100
        counter = 0
        dep = {}
        # ip_depth = ip_scale = input_depth = input_scale = 0
        with open(f[0], 'rt') as fi:
            for line in fi:
                if not line.startswith('#'): 
                    continue
                if counter > 100: # nrows
                    break # stop
                num = line.strip().split()[-1]
                if 'tags after filtering in treatment' in line:
                    dep['ip_depth'] = num
                if 'tags in treatment' in line:
                    s = 1e6 / int(num)
                    dep['ip_scale'] = '{:.6f}'.format(s)
                if 'tags after filtering in control' in line:
                    dep['input_depth'] = num
                if 'tags in control' in line:
                    s = 1e6 / int(num)
                    dep['input_scale'] = '{:.6f}'.format(s)
                counter += 1

        return dep



