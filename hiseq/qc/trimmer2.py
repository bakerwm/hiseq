
# -*- coding: utf-8 -*-

"""Quality control for fastq files

1. Trim N-bases (optional)
2. Trim adapter (guess adapter, see: trim_galore)
3. Trim N-bases (optional)

optional
1. --rm-untrim, --save-too-short, ...

"""


import os
import sys
import re
import shutil
import logging

from hiseq.utils.args import args_init
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions


@Logger('INFO')
class Trimmer(object):
    """Processing fastq files

    # for only 1 input fastq: (for smp_name)

    General
      a. low quality bases at 3'
      b. trim-n

    1. CLIP reads

      a. trim adapter (sliding)
      b. cut inline-barcode
      c. cut N-bases (5', 3' ends)
      e. collapse (remove duplicates)
      e. cut random barcode

    2. RNAseq reads

      a. trim adapter
      b. cut N-bases

    3. ChIPseq reads
      a. trim adapter
      b. cut N-bases (5', 3' ends)

    4. smRNAseq reads
      a. trim adapter
      b. discard untrimmed reads
      c. cut N-bases

    5. ATACseq reads
      a. trim adapter (sliding)

    """

    def __init__(self, **kwargs):
        """
        Main function
        """
        self.args = args_init(kwargs, trim=True)
        self.init_trimmer()


    def init_trimmer(self):
        # global vars
        self.fq1 = self.args.get('fq1', None) # checked in ArgumentsInit
        self.fq2 = self.args.get('fq2', None) # as above
        self.outdir = self.args.get('outdir', None) # as above
        assert isinstance(self.fq1, str)
        assert isinstance(self.outdir, str)
        assert is_path(self.outdir)

        # optional
        self.fqname = file_prefix(self.fq1)[0]
        if not self.fq2 is None:
            self.fqname = re.sub('[._][rR]?1$', '', self.fqname)
            self.fqname = re.sub('_1$', '', self.fqname)
        if self.args.get('smp_name', None):
            self.fqname = self.args['smp_name']

        self.fq_out_prefix = os.path.join(self.outdir, self.fqname) #self.args['fq_out_prefix']
        self.fq_out = self.fq_out_prefix + '.fq'    # SE
        self.fq_out1 = self.fq_out_prefix + '_1.fq' # PE1
        self.fq_out2 = self.fq_out_prefix + '_2.fq' # PE2

        # gzip
        if self.args.get('gzip', False):
            self.fq_out += '.gz'
            self.fq_out1 += '.gz'
            self.fq_out2 += '.gz'

        # files
        se_out = [self.fq_out, None]
        pe_out = [self.fq_out1, self.fq_out2]
        self.out_files = se_out if self.fq2 is None else pe_out


    def report(self):
        """
        Report the number of reads in each step
        """
        args = self.args.copy()
        rpt_file = self.fq_out_prefix + '.qc.stat'
        n_input = Fastx(self.fq1).count
        n_output = Fastx(self.out_files[0]).count # first one of output

        rpt_line = '#sample\tinput\toutput\tpercent\n'
        rpt_line += '{}\t{:d}\t{:d}\t{:.2f}%'.format(
            self.fqname,
            int(n_input),
            int(n_output),
            n_output / n_input * 100)

        with open(rpt_file, 'wt') as w:
            w.write(rpt_line + '\n')

        return [n_input, n_output]


    def wrap(self):
        """
        1. save log files: log/
        2. remove temp files: outdir/temp
        3. gzip output files:
        """
        args = self.args.copy()

        ## output
        log_dir = os.path.join(self.outdir, 'log')
        is_path(log_dir)

        log_cutadapt, log_rmdup, log_cut2 = [
            os.path.join(log_dir, self.fqname + i) for i in ['.cutadapt.log', '.rmdup.log', '.rmdup_cut.log']]

        ## 1. log file, cutadapt log        
        log_cutadapt_from = os.path.join(self.outdir, 'temp', '01.cutadapt_output', self.fqname + '.cutadapt.log')
        shutil.copy(log_cutadapt_from, log_cutadapt)

        ## 2. rmdup, cut
        with open(log_rmdup, 'wt') as w:
            w.write('remove_duplicate: {}'.format(args['rmdup']))

        ## 3. cut2
        with open(log_cut2, 'wt') as w:
            w.write('cut_after_trim: {}'.format(args['cut_after_trim']))

        ## 4. temp files
        temp_dir = os.path.join(self.outdir, 'temp')
        if not args['save_temp'] is True:
            shutil.rmtree(temp_dir, ignore_errors=True)


    def check(self):
        """
        Check arguments, target file exists, 
        """
        args = self.args.copy()

        ## save parameters
        args_pickle = self.fq_out_prefix + '.arguments.pickle'

        ## chk
        chk1 = args_checker(args, args_pickle)
        chk2 = args['overwrite'] is False

        if self.fq2 is None:
            chk3 = os.path.exists(self.out_files[0]) # SE
        else:
            chk3 = os.path.exists(self.out_files[0]) and os.path.exists(self.out_files[1]) # PE

        return all([chk1, chk2, chk3])


    def run(self):
        args = self.args.copy()

        if self.check():
            log.warning('{:>20} : file exists, skipped ...'.format(self.fqname))
            return self.out_files
        else:
            args_file = self.fq_out_prefix + '.arguments.txt'
            args_logger(args, args_file, True) # update arguments.txt

        ## 1.cutadapt directory
        args['outdir'] = os.path.join(self.outdir, 'temp', '01.cutadapt_output')
        f1, f2 = Cutadapt(**args).run()

        ## 2.rmdup
        rmdup_dir = os.path.join(self.outdir, 'temp', '02.rmdup_output')
        rmdup_f1 = os.path.join(rmdup_dir, os.path.basename(f1))
        rmdup_f2 = None if f2 is None else os.path.join(rmdup_dir, os.path.basename(f2))

        ## 3.cut-after-trim
        cut2_dir = os.path.join(self.outdir, 'temp', '03.cut_after_trim')
        cut2_f1 = os.path.join(rmdup_dir, 'cut.' + os.path.basename(f1))
        cut2_f2 = None if f2 is None else os.path.join(rmdup_dir, 'cut.' + os.path.basename(f2))

        ## create dir
        assert is_path(rmdup_dir)
        assert is_path(cut2_dir)

        # SE
        if self.fq2 is None:
            ## rmdup
            if args['rmdup']:
                Fastx(f1).collapse(rmdup_f1, fq_out=True)
            else:
                rmdup_f1 = f1

            ## cut2
            if eval(args['cut_after_trim']):
                args['cut'] = args.get('cut_after_trim', 0) # eg: 5,-3
                args['cut_to_length'] = args.get('cut_to_length', 0) # eg: 5,-3
                args['discard_tooshort'] = args.get('discard_tooshort', True)
                Fastx(rmdup_f1).cut(cut2_f1, **args)
            else:
                cut2_f1 = rmdup_f1

            ## clean
            if args['gzip']:
                gzip_cmd(cut2_f1, self.fq_out, decompress=False)
            else:
                shutil.move(cut2_f1, self.fq_out)
        # PE
        else:
            ## rmdup
            if args['rmdup']:
                Fastx(f1).collapse(rmdup_f1, fq_out=True)
                Fastx(f2).collapse(rmdup_f2, fq_out=True)
            else:
                rmdup_f1 = f1
                rmdup_f2 = f2

            # cut2
            if eval(args['cut_after_trim']): # default: '0'
                args['cut'] = args.get('cut_after_trim', 0) # eg: 5,-3
                args['cut_to_length'] = args.get('cut_to_length', 0) # eg: 5,-3
                args['discard_tooshort'] = args.get('discard_tooshort', True)
                if args['rmdup']:
                    ## if rmdup + cut2
                    Fastx(rmdup_f1).cut(cut2_f1, **args)
                    Fastx(rmdup_f2).cut(cut2_f2, **args)
                else:
                    ## if only cut2, require: PE pairs
                    Fastx(rmdup_f1).cut_pe(input2=rmdup_f2, out1=cut2_f1, out2=cut2_f2, **args)
            else:
                cut2_f1 = rmdup_f1
                cut2_f2 = rmdup_f2

            ## clean
            if args['gzip']:
                gzip_cmd(cut2_f1, self.fq_out1, decompress=False)
                gzip_cmd(cut2_f2, self.fq_out2, decompress=False)
            else:
                shutil.move(cut2_f1, self.fq_out1)
                shutil.move(cut2_f2, self.fq_out2)

        # wrap
        self.wrap()

        # report
        self.n_input, self.n_output = self.report()

        return self.out_files


@Logger('INFO')
class Cutadapt(object):

    # def __init__(self, fq, outdir, fq2=None, **kwargs):
    def __init__(self, **kwargs):
        """
        Trim adapter and low quality bases for fastq file(s) using ``cutadapt``
        program version 2.7. the program will guess the adapters: TruSeq,
        Nextera, small RNAseq, ...


        *fq* is typically the name of fastq file

        *outdir* is the path to the directory saving the results

        *fq2* is the second read file of PE sequencing, optinoal

        Usage is to specify the fastq file, adapter and other arguments:

            Cutadapt('read1.fq', 'output', 'read2.fq')

        """
        # prepare args
        self.args = ArgumentsInit(kwargs, trim=True).dict.__dict__
        self.args.pop('args_input', None)
        self.args.pop('cmd_input', None)
        self.args.pop('dict', None)
        self.init_cut()


    def init_cut(self):
        """
        Check arguments: fq, fq2, outdir
        """
        args = self.args.copy()

        # global vars
        self.fq1 = self.fq = args.get('fq1', None) # checked in ArgumentsInit
        self.fq2 = args.get('fq2', None) # as above
        self.outdir = args.get('outdir', None) # as above

        assert isinstance(self.fq1, str)
        assert isinstance(self.outdir, str)
        assert is_path(self.outdir)

        # optional
        self.fqname = file_prefix(self.fq1)[0]
        if not self.fq2 is None:
            self.fqname = re.sub('[._][rR]?1$', '', self.fqname)
            self.fqname = re.sub('_1$', '', self.fqname)
        if args.get('smp_name', None):
            self.fqname = args['smp_name']

        self.fq_out_prefix = os.path.join(self.outdir, self.fqname) #self.args['fq_out_prefix']
        self.fq_out = self.fq_out_prefix + '.fq'    # SE
        self.fq_out1 = self.fq_out_prefix + '_1.fq' # PE1
        self.fq_out2 = self.fq_out_prefix + '_2.fq' # PE2
        self.log  = self.fq_out_prefix + '.cutadapt.log'
        self.cutadapt = shutil.which('cutadapt') # in PATH

        # files
        se_out = [self.fq_out, None]
        pe_out = [self.fq_out1, self.fq_out2]
        self.out_files = se_out if self.fq2 is None else pe_out


    def cut_cut(self, cutadapt_arg=True, read2=False):
        """
        Number of bases to cut from left or right of sequence
        5: from left
        -3: from right
        for cutadapt command:
        --cut={cut}
        """
        args = self.args.copy()

        cut_before_trim = args['cut_before_trim']
        cut_opts = cut_before_trim.split(',') if ',' in cut_before_trim else [cut_before_trim]

        if cutadapt_arg:
            if read2:
                p = '-U' # read2
            else:
                p = '-u' # read1
            out = ' '.join(['{} {}'.format(p, i) for i in cut_opts])
        else:
            out = cut_opts

        return out


    def adapter_sliding(self, step=2, window=15):
        """
        For some situation, the adapter (3') in library was differ, especially
        for ligation strategy libraries, (eg: iCLIP, eCLIP).
        We need to make sliding windows of the adatpers for adapter removing
        """
        args = self.args.copy()

        ## always, the 3' adapter
        adapter = args['adapter3']

        ## sliding
        adapter_sliders = [adapter[i:i+window] for i in range(0, len(adapter)-window, step) if i >= 0]
        if not adapter_sliders:
            adapter_sliders = [adapter] # full length

        return adapter_sliders


    def get_cmd(self):
        """
        Create basic command line for cutadapt program: SE read

        -a <ad> -m 15 -q 20 --trim-n --max-n=0.1 --error-rate=0.1

        """
        args = self.args.copy()

        # 3' adapter
        ad3_list = self.adapter_sliding() if args['adapter_sliding'] else [args['adapter3']]
        arg_ad3 = ' '.join(['-a {}'.format(i) for i in ad3_list])

        # SE mode
        if self.fq2 is None:
            arg_AD3 = '' # SE
            arg_out = '-o {}'.format(self.fq_out)
            # cut
            if args['cut_before_trim']:
                arg_cut = self.cut_cut()

            # untrimmed
            if args['save_untrim']:
                arg_untrim = '--untrimmed-output={}'.format(self.fq_out_prefix + '.untrimmed.fq')
                args['threads'] = 1
            elif args['rm_untrim']:
                arg_untrim = '--discard-untrimmed'
            else:
                arg_untrim = '' # save untrimmed

            # too-short
            if args['save_too_short']:
                arg_short = '--too-short-output={}'.format(self.fq_out_prefix + '.too_short.fq')
            else:
                arg_short = ''

            # too-long
            if args['save_too_long']:
                arg_long = '--too-long-output={}'.format(self.fq_out_prefix + '.too_long.fq')
            else:
                arg_long = ''

            # cut to length
            if args['cut_to_length'] >= args['len_min']:
                arg_cut_to_length = '--length={}'.format(args['cut_to_length'])
            else:
                arg_cut_to_length = ''

        else:
            AD3_list = self.adapter_sliding(args['AD3']) if args['adapter_sliding'] else [args['AD3']]
            arg_AD3 = ' '.join(['-A {}'.format(i) for i in AD3_list])
            arg_out = '-o {} -p {}'.format(self.fq_out1, self.fq_out2)

            # cut
            if args['cut_before_trim']:
                arg_cut1 = self.cut_cut()
                arg_cut2 = self.cut_cut(read2=True)
                arg_cut = arg_cut1 + ' ' + arg_cut2

            # untrimmed
            if args['save_untrim']:
                arg_untrim = '--untrimmed-output={} \
                    --untrimmed-paired-output={}'.format(
                        self.fq_out_prefix + '.untrimmed_1.fq',
                        self.fq_out_prefix + '.untrimmed_2.fq')
                args['threads'] = 1
            elif args['rm_untrim']:
                arg_untrim = '--discard-untrimmed'
            else:
                arg_untrim = '' # save untrimmed

            # too-short
            if args['save_too_short']:
                arg_short = '--too-short-output={} \
                    --too-short-paired-output={}'.format(
                        self.fq_out_prefix + '.too_short_1.fq',
                        self.fq_out_prefix + '.too_short_2.fq')
            else:
                arg_short = ''

            # too-long
            if args['save_too_long']:
                arg_long = '--too-long-output={} \
                    --too-long-paired-output={}'.format(
                        self.fq_out_prefix + '.too_long_1.fq',
                        self.fq_out_prefix + '.too_long_2.fq')
            else:
                arg_long = ''

            # cut to length
            if args['cut_to_length'] >= args['len_min']:
                arg_cut_to_length = '--length={}'.format(args['cut_to_length'])
            else:
                arg_cut_to_length = ''

        ## save log
        arg_log = ''
        ## command line

        # command line
        arg_main = '{} -m {} -q {} --trim-n --max-n=0.1 --error-rate={} \
            --times={} --cores {}'.format(
                self.cutadapt,
                args['len_min'],
                args['qual_min'],
                args['error_rate'],
                args['trim_times'],
                args['threads'])

        cmd_line = ' '.join([
            arg_main,
            arg_ad3,
            arg_AD3,
            arg_cut,
            arg_untrim,
            arg_short,
            arg_long,
            arg_cut_to_length])

        return cmd_line


    def run_se(self):
        args = self.args.copy()

        cmd_line = self.get_cmd()

        if os.path.exists(self.fq_out) and args['overwrite'] is False:
            log.warning('{:>20} : file exists, skipped ...'.format(self.fqname))
        else:
            cmd_line += ' -o {} {} 1>{}'.format(
                self.fq_out,
                self.fq,
                self.log)
            run_shell_cmd(cmd_line)

        return [self.fq_out, None]


    def run_pe(self):
        args = self.args.copy()

        cmd_line = self.get_cmd()

        if os.path.exists(self.fq_out1) and os.path.exists(self.fq_out2) and args['overwrite'] is False:
            log.warning('{:>20} : file exists, skipped ...'.format(self.fqname))
        else:
            cmd_line += ' -o {} -p {} {} {} 1>{}'.format(
                self.fq_out1,
                self.fq_out2,
                self.fq1,
                self.fq2,
                self.log)
            run_shell_cmd(cmd_line)

        return [self.fq_out1, self.fq_out2]


    def run(self):
        args = self.args.copy()

        # SE
        if self.fq2 is None:
            fq_out = self.run_se()
        # PE
        else:
            fq_out = self.run_pe()

        return fq_out

