"""
Functions around reading and classifying sequence
"""
from __future__ import print_function
from pycbio.hgdata import dnaOps


class GenomeReader(object):
    """
    Reads sequences from twoBitReader, which can be substituted with a mock
    version.
    """
    def __init__(self, twoBitReader):
        self.twoBitReader = twoBitReader
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

    def getChromSize(self, chrom):
        self.__obtainChrom(chrom)
        return self.size
