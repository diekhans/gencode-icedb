"""
Functions around reading and classifying sequence
"""
from __future__ import print_function
import os
from pycbio.hgdata import dnaOps
from twobitreader import TwoBitFile

# FIXME: is factory really needed?  Derived classes might be better,
# interaction with command line and subprocess is part of this.

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

    def obtain(self):
        return GenomeReader(TwoBitFile(self.twoBitFile))

    def getOptionArgs(self):
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

    @staticmethod
    def factoryFromUcscDb(ucscDb):
        """create a factory from a file in ucsc locations, checking for cluster location"""
        tb = "/scratch/data/{ucscDb}/{ucscDb}.2bit".format(ucscDb=ucscDb)
        if not os.path.exists(tb):
            tb = "/hive/data/genomes/{ucscDb}/{ucscDb}.2bit".format(ucscDb=ucscDb)
        return GenomeReaderFactory(tb)
