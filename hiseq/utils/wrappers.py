
# -*- coding: utf-8 -*-

"""
Wrappers for qc, alignment, ..., create html report
"""

from hiseq.utils.helper import *

@Logger('INFO')
class ReportFastqc(object):
    """
    Organize the fastqc output files into one single html report
    """
    def __init__(self, path):
        """
        Input the directory, containing fastq files
        """
        assert os.path.isdir(path)

        ## search all fastq files, in current dir
        pattern = re.compile('.f(ast)?[aq](.gz)?$')        
        fq_files = listfiles(self.path, recursive=False)
        fq_files = [f for f in fq_files if re.search(pattern, f)]

        self.path = path
        self.fq_files = fq_files


    def fastqc(self, run_fastqc=False, path=None, overwrite=False):
        """
        *_fastqc.zip and *_fastqc.html files
        if not, run the fastqc command, with default parameters
        """
        ## output of fastqc
        if isinstance(path, str) and os.path.isdir(path):
            fqpath = path
        else:
            fqpath = self.path

        if not os.path.exists(fqpath):
            os.makedirs(fqpath)

        qc_files = []
        for fq in self.fq_files:
            fqname = file_prefix(fq)[0]
            z = os.path.join(fqpath, fqname + '_fastqc.zip')
            h = os.path.join(fqpath, fqname + '_fastqc.html')

            if os.path.exists(z) and os.path.exists(h) and overwrite is False:
                log.info('{:>20} : file exists, fastqc skipped ...'.format(fqname))
            else:
                fastqc_exe = shutil.which('fastqc')
                cmd = '{} -o {} {}'.format(fastqc_exe, fqpath, fq)
                if run_fastqc is True:
                    run_shell_cmd(cmd)
            qc_files.append(z)

        return fqpath


    def fastqcr(self, qcpath=None, template=''):
        """
        Generate html report for all the fastqc files
        using fastqcr R package
        
        package_dir/
        """
        ## output of fastqc
        if isinstance(path, str) and os.path.isdir(path):
            qcpath = path
        else:
            qcpath = self.path

        ## template
        ## convert to absolute path
        if isinstance(template, str) and os.path.exists(template):
            template = os.path.abspath(template)
        else:
            template = ''

        ## out html
        out_html = os.path.join(qcpath, 'qc_report.html')

        ## current script, (hiseq.utils)
        this_file = os.path.realpath(__file__)
        qc_r = os.path.join(os.path.dirname(this_file), 'qc_report.R')

        ## command line
        rscript_exe = shutil.which('Rscript')
        cmd = '{} {} {}'.format([
            rscript_exe,
            qc_r,
            qcpath,
            template])

        if not os.path.exists(out_html):
            run_shell_cmd(cmd)

        return out_html


    def run(self, run_fastqc=False, outdir=None, overwrite=False):
        """
        Run fastqc and fastqcr for all fastq files
        """
        ## run fastqc for all fastq files
        qcpath = os.path.join(self.path, 'qc')
        self.fastqc(True, qcpath, overwrite)

        ## run fastqcr 
        out_html = self.fastqcr(qcpath)

        log.info('{:>20} : save fastqc results to file')

        return out_html


@Logger('INFO')
class ReportAlignment(object):
    """
    Generate html report for alignment directories
    search for *.stat files
    """

    def __init__(self, path):
        """
        search *mapping_stat.csv files, or list of files
        """
        pass

        
        