
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
        atac=False, overwrite=False):
        """Parse the parameters
        venv, the virtualenv created for macs2, running in Python2
        """
        self.ip = ip
        self.genome = genome
        self.output = output
        self.control = control
        self.overwrite = overwrite
        self.prefix = prefix
        self.atac = atac

        if prefix is None:
            prefix = file_prefix(ip)[0]
            # prefix = os.path.splitext(os.path.basename(ip))[0]

        self.gsize = self.get_gsize()
        if self.gsize is None:
            raise ValueError('unknown genome: {}'.format(genome))

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

        return d.get(self.genome, None)


    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.out',
                                            delete=False)
        return tmpfn.name


    def callpeak(self):
        """Call peaks using MACS"""
        # genome_size = self.get_gsize()
        log = os.path.join(self.output, self.prefix + '.macs2.callpeak.out')
        if self.atac is True:
            # ATAC-seq
            macs2_cmd = 'macs2 callpeak --nomodel --shift -100 --extsize 200 \
                -f BAM -t {} -g {} --outdir {} -n {} 2> {}'.format(
                    self.ip,
                    self.gsize,
                    self.output,
                    self.prefix,
                    log)
        else:
            # ChIP-seq
            macs2_cmd = 'macs2 callpeak -f BAM -t {} -g {} --outdir {} -n {} \
                --keep-dump auto -B --SPMR'.format(
                    self.ip, 
                    self.gsize,
                    self.output,
                    self.prefix)
            if not self.control is None:
                macs2_cmd += ' -c {}'.format(self.control)

            macs2_cmd += ' 2> {}'.format(log)

       # output file
        macs2_out = os.path.join(self.output, self.prefix + '_peaks.xls')
        if os.path.exists(macs2_out) and self.overwrite is False:
            log.info('file exists, skip macs2 callpeak')
        else:
            logging.info('run macs2 callpeak')
            run_shell_cmd(macs2_cmd)


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
            raise ValueError('unknown option: opt=%s' % opt)

        if opt == 'logLR':
            opt_ext = '-p 0.00001'
        else:
            opt_ext = ''

        out_bdg = os.path.join(self.output, self.prefix + '.' + opt + '.bdg')
        log = os.path.join(self.output, self.prefix + '.macs2.bdgcmp.out')
        c = "macs2 bdgcmp -t %s -c %s -o %s -m %s %s" % (ip_bdg, input_bdg, out_bdg, opt, opt_ext)
        
        if os.path.exists(out_bdg) and self.overwrite is False:
            logging.info('file exists, skip macs2 bdgcmp')
        else:
            self.python2_run(c)

        # sort output *.bdg
        c1 = 'sort -k1,1 -k2,2n -o %s %s' % (out_bdg, out_bdg)

        # cnvert *.bdg to *.bigWig
        gsize_file = Genome(self.genome).get_fasize()
        out_bw  = os.path.join(self.output, self.prefix + '.' + opt + '.bigWig')
        c2 = 'bedGraphToBigWig %s %s %s' % (out_bdg, gsize_file, out_bw)
        if os.path.exists(out_bw) and self.overwrite is False:
            logging.info('file exists, skip bg2bw')
        else:
            subprocess.run(shlex.split(c1)) # sort bdg
            subprocess.run(shlex.split(c2)) # convert bdg to bw


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

        # narrow peak file
        peak_listA = listfiles2('*_peaks.narrowPeak', self.output)
        peak_listB = listfiles2('*_peaks.broadPeak', self.output)
        peak_list = peak_listA + peak_listB

        anno_list = []
        for peak_file in peak_list:
            peak_anno = peak_file + '.annotation'
            peak_log = peak_anno + '.log'
            cmd = 'perl {} {} {} 1> {} 2> {}'.format(
                anno_exe, peak_file, self.genome, peak_anno, peak_log)
            # check existence
            if os.path.exists(peak_anno) and self.overwrite is False:
                log.info('file exists, skip annotation: {}'.format(peak_anno))
            else:
                run_shell_cmd(cmd)

            anno_list.append(peak_anno)

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
            raise ValueError('file not found, need to run .callpeak() first - %s' % broadpeak)

        # run
        anno_peak = os.path.join(self.output, self.prefix + '_peaks.broadPeak.annotation')
        anno_log = os.path.join(self.output, self.prefix + '_peaks.broadPeak.annotation.log')
        cmd = 'perl {} {} {} 1> {} 2> {}'.format(
            anno_exe, broadpeak, self.genome, anno_peak, anno_log)

        run_shell_cmd(cmd)
        
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
        f = listfiles2("*peaks.xls", self.output)
        if not os.path.exists(f[0]):
            raise ValueError('file missing in macs2 callpeak output: %s' % f)

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
                    dep['ip_scale'] = '%.6f' % s
                if 'tags after filtering in control' in line:
                    dep['input_depth'] = num
                if 'tags in control' in line:
                    s = 1e6 / int(num)
                    dep['input_scale'] = '%.6f' % s
                counter += 1

        return dep



