
"""
Calculate the fragment size of PE reads, 

Available tools: CollectInsertSizeMetrics (Picard), bamPEFragmentSize (deeptools)

Here are the samtools method:


Statistics:
1. length count
2. mean/medium/qrt. 

"""

import os
import sys
import pathlib
import numpy as np
import pandas as pd
import pysam
import tempfile
import logging
import hiseq
from multiprocessing import Pool


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)

log.setLevel('INFO')


class BamPEFragSize(object):
    """Calculate insert size of PE fragment
    single BAM input

    sample size = 10000
    """
    def __init__(self, bam, labels=None, maxRecords=0):
        self.bam = bam
        # self.labels = labels
        self.maxRecords = maxRecords
        # name
        bam_name = os.path.splitext(os.path.basename(bam))[0]
        self.labels = labels if labels else bam_name
        # calculation
        self.freqTable = self.calFragSize() # dataframe


    def isBam(self):
        """Check input is BAM file
        update: self.bam

        string
        *.bam
        """
        flag = False
        if isinstance(self.bam, str):
            flag = self.bam.endswith('.bam') and os.path.exists(self.bam)

        if not flag:
            raise ValueError('bam file expected, {} got'.format(type(self.bam).__name__))

        return flag


    def isBamPE(self, topn=10000):
        """Check input bam is Paired end alignment"""
        flag = False
        if self.isBam():
            samfile = pysam.AlignmentFile(self.bam)
            flag = [read.is_paired for read in samfile.head(topn)]
            samfile.close()

        return all(flag)


    def calFreq(self, x, return_dataframe=True):
        """Calculate the frequency of list
        return dataframe

        index count
        """
        if isinstance(x, list):
            var, freq = np.unique(x, return_counts=True)
            df = pd.DataFrame(data=freq, index=var, columns=['count'])
        else:
            df = pd.DataFrame(columns=['count'])

        return df


    def calFragSize(self, bam=None, chunk=1000000):
        """Extract the fragment length of paired end alignment
        return

        length count id
        """
        if bam is None:
            bam = self.bam

        if not self.isBamPE():
            raise ValueError('PE bam expected, failed')

        sam = pysam.AlignmentFile(bam)

        # temp
        # dataframe
        flag = 0
        fragSizes = []
        df = pd.DataFrame(columns=['count'])

        for read in sam:
            if read.is_proper_pair \
            and not read.is_unmapped \
            and not read.mate_is_unmapped \
            and not read.is_read1 \
            and not read.is_duplicate \
            and read.template_length > 0:
                flag += 1
                fragSizes.append(read.template_length)

            # get sample size
            if self.maxRecords > 0 and flag > self.maxRecords:
                log.info('stop iteration at: {}'.format(flag))
                break # stop

            # chunk
            if flag > 0 and flag % chunk == 0:
                dfn = self.calFreq(fragSizes)
                df = pd.concat([df, dfn], axis=1).sum(axis=1)
                fragSizes = [] # empty
                log.info('{} : {} {}'.format('Processed', flag, self.labels))

        # last chunk
        if len(fragSizes) > 0:
            dfn = self.calFreq(fragSizes)
            df = pd.concat([df, dfn], axis=1).sum(axis=1)
            fragSizes = [] # empty
            log.info('{} : {} {}'.format('Processed', flag, self.labels))

        # convert to data.frame
        df = df.reset_index()
        df.columns = ['length', 'count']
        # bam_name = os.path.splitext(os.path.basename(bam))[0]
        # df['id'] = self.labels if self.labels else bam_name
        df['id'] = self.labels

        return df


    def distribution(self):
        """Basic statistics values
        value + freq

        mean, medium, mode, std, min, max, Q1, Q2, Q3
        """
        val = self.freqTable['length']
        freq = self.freqTable['count']
        inserts = np.repeat(val, freq)
        # inserts = np.repeat(self.freqTable['length'], self.freqTable['count'])

        # statistics
        q_mean = np.mean(inserts)
        q_median = np.median(inserts)
        q_median_dev = np.median(np.absolute(inserts - q_median))
        q_mode = val[np.argmax(freq)]
        q_std = np.std(inserts)
        q_min = np.min(inserts)
        q_max = np.max(inserts)
        q_qual = np.quantile(inserts, [0.25, 0.5, 0.75], axis=0)

        # core distribution
        s = np.array([q_mean, q_median, q_mode, q_std, q_min, q_max])
        s = np.append(s, q_qual)

        return s


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.csv',
            delete=False)
        return tmp.name


    def saveas(self, csv_file=None):
        """Save to file"""
        if csv_file is None:
            csv_file = self._tmp()

        log.info('saving to file: {}'.format(csv_file))

        try:
            # save table
            self.freqTable.to_csv(csv_file, index=False)

            # save statistics
            stat_file = csv_file + '.stat'
            da = self.distribution()
            pd.DataFrame(da).T.to_csv(stat_file, sep="\t", index=False, 
                header=['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1', 'Q2', 'Q3'])
        except:
            log.warning('failed saving file: {}'.format(csv_file))


    def plot(self, plot_file):
        """Generate freq table plot
        line plot, hist
        """
        # save to csv file
        csv_file = os.path.splitext(plot_file)[0] + '.csv'
        self.saveas(csv_file) # save csv, stat

        # save to plot
        # to-do: matplotlib function
        hiseq_dir = os.path.dirname(hiseq.__file__)
        fragPlotR = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
        cmd = ' '.join(['Rscript', fragPlotR, plot_file, csv_file])

        try:
            os.system(cmd)
        except:
            log.warning('fragsize.py failed')


    def run(self):
        """Generate dataframe, pdf
        bam, str or list
        """
        if isinstance(self.bam, str):
            bam_list = [self.bam]
        elif isinstance(self.bam, list):
            bam_list = self.bam
        else:
            raise Exception('str and list expected, {} got'.format(type(self.bam).__name__))

        # create df
        frames = [self.calFragSize(b) for b in bam_list]
        df = pd.concat(frames, axis=0)

        # save to txt
        df.to_csv(frag_stat, index=False)

        # save to plot
        hiseq_dir = os.path.dirname(hiseq.__file__)
        fragPlotR = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
        cmd = ' '.join(['Rscript', fragPlotR, frag_stat, frag_pdf])

        try:
            os.system(cmd)
        except:
            log.warning('fragsize.py failed')


class BamPEFragSize2(object):
    """for command-line 
    multiple BAM files
    rename (labels)
    """
    def __init__(self, **kwargs):
        """Options:
        bam
        labels
        pdf
        csv
        """
        self.update(kwargs, force=True) # fresh new
        self.config() # update, default args


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)    


    def config(self):
        """Check input bam files, labels, ...
        bam, 
        labels,
        pdf
        csv
        maxRecords
        maxLength
        ...
        """
        args_init = {
            'bam': None,
            'outdir': str(pathlib.Path.cwd()),
            'labels': None,
            'maxRecords': 0, # all
            'maxLength': 1000, #
            'pdf_out': None,
            'csv_out': None,
            'threads': 4,
            'overwrite': False
        }
        self.update(args_init, force=False) # default

        # bam_list
        if isinstance(self.bam, str):
            self.bam = [self.bam]
        elif isinstance(self.bam, list):
            pass
        else:
            raise ValueError('bam, str or list expected, {} got'.format(type(self.bam)))

        # outdir
        if not os.path.exists(self.outdir):
            try:
                os.makedirs(self.outdir)
            except IOError:
                log.error('failed to create directory: {}'.format(self.outdir))

        # labels
        bam_names = [os.path.splitext(os.path.basename(b))[0] for b in self.bam]
        if self.labels is None:
            self.labels = bam_names
        else:
            self.labels = self.labels[:len(bam_names)] + bam_names[len(self.labels):]


    def bamSingle(self, bam):
        """Run for single bam
        stat
        pdf
        require outdir
        """
        bam_i = self.bam.index(bam) # index for bam file
        labels = self.labels[bam_i] # labels for bam

        # output file
        # csv_file = os.path.join(self.outdir, labels + '.frag_size.csv')
        pdf_file = os.path.join(self.outdir, labels + '.frag_size.pdf')

        if os.path.exists(pdf_file) and not self.overwrite:
            log.info('file exists, frag_size() skipped, {}'.format(pdf_file))
        else:
            BamPEFragSize(bam, labels=labels).plot(pdf_file) # save to csv


    def run(self):
        # run for single
        with Pool(processes=self.threads) as pool:
            pool.map(self.bamSingle, self.bam)

        # plot merge
        csv_file_list = [os.path.join(self.outdir, i + '.frag_size.csv') for i in self.labels]
        csv_file_list = [i for i in csv_file_list if os.path.exists(i)]
        pdf_merge = os.path.join(self.outdir, 'frag_size.pdf')

        # save to plot
        # to-do: matplotlib function
        hiseq_dir = os.path.dirname(hiseq.__file__)
        fragPlotR = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
        cmd = ' '.join(['Rscript', fragPlotR, pdf_merge] + csv_file_list)

        try:
            os.system(cmd)
        except:
            log.warning('fragsize.py failed')


def main():
    if len(sys.argv) < 3:
        sys.exit('fragsize.py <outdir> <bam1> [bam2, ...]')

    outdir = sys.argv[1]
    bam_list = sys.argv[2:]

    BamPEFragSize2(bam=bam_list, outdir=outdir).run()


if __name__ == '__main__':
    main()


def frag_size_picard(bam, out):
    """Calculate PE fragment size, insertion size using Picard
    Example:

    java -jar picard.jar CollectInsertSizeMetrics \
      I=input.bam \
      O=insert_size_metrics.txt \
      H=insert_size_histogram.pdf \
      M=0.5
    """
    pass


def frag_size_deeptools(bam, out):
    """Calculate PE fragment size, using deeptools
    Example:

    deepTools2.0/bin/bamPEFragmentSize \
    -hist fragmentSize.png \
    -T "Fragment size of PE RNA-seq data" \
    --maxFragmentLength 1000 \
    -b testFiles/RNAseq_sample1.bam testFiles/RNAseq_sample2.bam \
    testFiles/RNAseq_sample3.bam testFiles/RNAseq_sample4.bam \
    -samplesLabel sample1 sample2 sample3 sample4
    """
    pass


def frag_size_samtools(infile, outfile=None):
    """
    Calculate Fragment length from BAM file
    BAM: column-9

    BAM is PE
    BAM is sorted
    BAM is indexed
    """
    d = {}
    for read in pysam.AlignmentFile(infile):
        if not read.is_proper_pair:
            continue
        if not read.is_read1:
            continue
        size = abs(read.tlen)
        d[size] = d.get(size, 0) + 1
    d = sorted(d.items(), key=lambda x: x[0])

    if not outfile is None:
        with open(outfile, 'wt') as w:
            for k, v in d:
                w.write(str(k) + '\t' + str(v) + '\n')

    return d

