"""
Functions around reading and classifying sequence
"""
from __future__ import print_function
import os
from pycbio.sys import fileOps
from pycbio.tsv import TsvReader
from pycbio.hgdata import dnaOps
from twobitreader import TwoBitFile


class GenomeReader(object):
    """
    Reads sequences from twoBitReader, which can be substituted with a mock
    version for testing.
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


class MockGenomeReader(object):
    """Fake GenomeReader for testing when 2bit isn't there"""
    def __init__(self, mockDataTsv):
        self.mockDataSizes = dict()
        self.mockDataSeqs = dict()
        for row in TsvReader(mockDataTsv,
                             typeMap={"start": int,
                                      "end": int,
                                      "size": int}):
            self.mockDataSizes[row.chrom] = row.size
            self.mockDataSeqs[(row.chrom, row.start, row.end, row.strand)] = row.seq

    def get(self, chrom, start, end, strand=None):
        if strand is None:
            strand = '+'
        return self.mockDataSeqs[(chrom, start, end, strand)]

    def getChromSize(self, chrom):
        return self.mockDataSizes[chrom]


class MockSeqWriter(object):
    """Wrapper around GenomeReader that create mock data from twoBit """
    def __init__(self, genomeReader, mockDataTsv):
        self.genomeReader = genomeReader
        self.mockDataFh = open(mockDataTsv, "w")
        fileOps.prRowv(self.mockDataFh, "chrom", "start", "end", "size", "strand", "seq")

    def get(self, chrom, start, end, strand=None):
        if strand is None:
            strand = '+'
        assert strand == '+'
        size = self.genomeReader.getChromSize(chrom)
        seq = self.genomeReader.get(chrom, start, end, strand)
        fileOps.prRowv(self.mockDataFh, chrom, start, end, size, strand, seq)
        self.mockDataFh.flush()
        return seq

    def getChromSize(self, chrom):
        return self.genomeReader.getChromSize(chrom)


class GenomeReaderFactory(object):
    """Obtain the appropriate genome reader, either a reader of actually twobit,
    or a mock reader for tests."""
    def __init__(self, twoBitFile, mockTsv, updateMockReader, forceMockReader):
        self.twoBitFile = None
        if (twoBitFile is not None) and os.path.exists(twoBitFile):
            self.twoBitFile = twoBitFile
        self.mockTsv = None
        if (mockTsv is not None) and os.path.exists(mockTsv):
            self.mockTsv = mockTsv
        self.updateMockReader = updateMockReader
        self.forceMockReader = forceMockReader

        if (self.twoBitFile is None) and (self.mockTsv is None):
            raise Exception("need at least one of twoBitFile ({}) or mockTsv ({})".format(twoBitFile, mockTsv))
        if self.updateMockReader and ((self.twoBitFile is None) or (self.mockTsv is None)):
            raise Exception("need both of twoBitFile ({}) or mockTsv ({}) for updateMockReader".format(twoBitFile, mockTsv))
        if self.forceMockReader and (self.mockTsv is None):
            raise Exception("need mockTsv ({}) for forceMockReader".format(twoBitFile, mockTsv))
        self.genomeReader = None

    def __getReal(self):
        genomeReader = GenomeReader(TwoBitFile(twoBitHg19))
        if self.updateMockReader:
            genomeReader = MockSeqWriter(genomeReader, self.mockTsv)
        return genomeReader

    def __getMock(self):
        return MockGenomeReader(self.mockTsv)

    def obtain(self):
        if self.genomeReader is None:
            if (self.twoBitFile is not None) and not self.forceMockReader:
                self.genomeReader = self.__getReal()
            else:
                self.genomeReader = self.__getMock()
        return self.genomeReader
