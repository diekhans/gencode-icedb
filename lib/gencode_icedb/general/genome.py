"""
Functions around reading and classifying sequence
"""
from __future__ import print_function
import os
from pycbio.hgdata import dnaOps
from twobitreader import TwoBitFile

# twobitreader warning, pull request was submitted.
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module='twobitreader')


class GenomeReader(object):
    """
    Reads sequences from twoBitReader or perhaps other formats.
    """
    def __init__(self, twoBitReader):
        self.twoBitReader = twoBitReader
        self.chrom = self.size = None

    def close(self):
        "twobit package doesn't have close"
        self.twoBitReader = None
        self.chrom = self.size = None

    def __obtainChrom(self, chrom):
        if self.chrom != chrom:
            self.twoBitSeq = self.twoBitReader[chrom]
            self.chrom = chrom
            self.size = len(self.twoBitSeq)

    def get(self, chrom, start, end, strand=None):
        self.__obtainChrom(chrom)
        if strand == '-':
            start, end = dnaOps.reverseCoords(start, end, self.size)
        seq = self.twoBitSeq[start:end]
        if strand == '-':
            seq = dnaOps.reverseComplement(seq)
        return seq

    def haveChrom(self, chrom):
        return chrom in self.twoBitReader

    def getChroms(self):
        return sorted(self.twoBitReader.keys())

    def getChromSize(self, chrom):
        self.__obtainChrom(chrom)
        return self.size


class GenomeReaderFactory(object):
    """Obtain the appropriate genome reader."""
    def __init__(self, twoBitFile):
        self.twoBitFile = None
        if (twoBitFile is not None) and os.path.exists(twoBitFile):
            self.twoBitFile = twoBitFile
        if self.twoBitFile is None:
            raise Exception("no twoBitFile available")
        self.genomeReader = None

    def obtain(self):
        if self.genomeReader is None:
            self.genomeReader = GenomeReader(TwoBitFile(self.twoBitFile))
        return self.genomeReader

    def getOptionArgs(self, allowUpdate):
        "create a vector of options and values for passing to another program"
        args = []
        if self.twoBitFile is not None:
            args.append("--genomeSeqs={}".format(self.twoBitFile))
        return args

    @staticmethod
    def addCmdOptions(parser):
        """Add command options used to create a factory to the parser"""
        parser.add_argument('--genomeSeqs',
                            help="""Genome sequence twobit file to obtain splice sites""")

    @staticmethod
    def factoryFromCmdOptions(opts):
        """create a factory given options parse from command line"""
        # this will check sanity of options
        return GenomeReaderFactory(opts.genomeSeqs)
