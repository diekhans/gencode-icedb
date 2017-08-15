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

    def close(self):
        "doesn't do anything, but writer does use close"
        pass

    def get(self, chrom, start, end, strand=None):
        if strand is None:
            strand = '+'
        return self.mockDataSeqs[(chrom, start, end, strand)]

    def haveChrom(self, chrom):
        return chrom in self.mockDataSizes

    def getChroms(self):
        return sorted(self.mockDataSizes.keys())

    def getChromSize(self, chrom):
        return self.mockDataSizes[chrom]


class MockSeqWriter(object):
    """Wrapper around GenomeReader that create mock data from twoBit """
    def __init__(self, genomeReader, mockDataTsv):
        self.genomeReader = genomeReader
        self.mockDataFh = open(mockDataTsv, "w")
        fileOps.prRowv(self.mockDataFh, "chrom", "start", "end", "size", "strand", "seq")
        # track what we got so we can save chroms where only size was obtained
        # to mock file
        self.gotSizes = set()
        self.gotSeq = set()

    def __write(self, chrom, start, end, size, strand, seq):
        fileOps.prRowv(self.mockDataFh, chrom, start, end, size, strand, seq)

    def close(self):
        "make sure chromosomes where only the size is obtained are saved"
        for chrom in self.gotSizes - self.gotSeq:
            self.__write(chrom, 0, 0, self.getChromSize(chrom), "+", "")
        self.mockDataFh.close()
        self.mockDataFh = None
        self.genomeReader.close()

    def __del__(self):
        if self.mockDataFh is not None:
            self.close()

    def get(self, chrom, start, end, strand=None):
        self.gotSeq.add(chrom)
        if strand is None:
            strand = '+'
        assert strand == '+'
        size = self.genomeReader.getChromSize(chrom)
        seq = self.genomeReader.get(chrom, start, end, strand)
        self.__write(chrom, start, end, size, strand, seq)
        return seq

    def haveChrom(self, chrom):
        return self.genomeReader.haveChrom()

    def getChroms(self):
        return self.genomeReader.getChroms()

    def getChromSize(self, chrom):
        self.gotSizes.add(chrom)
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

        if self.updateMockReader and self.forceMockReader:
            raise Exception("can't specified both updateMockReader and forceMockReader")
        if (self.twoBitFile is None) and (self.mockTsv is None):
            raise Exception("need at least one of twoBitFile ({}) or mockTsv ({})".format(twoBitFile, mockTsv))
        if self.updateMockReader and ((self.twoBitFile is None) or (self.mockTsv is None)):
            raise Exception("need both of twoBitFile ({}) or mockTsv ({}) for updateMockReader".format(twoBitFile, mockTsv))
        if self.forceMockReader and (self.mockTsv is None):
            raise Exception("need mockTsv ({}) for forceMockReader".format(twoBitFile, mockTsv))
        self.genomeReader = None

    def __getReal(self):
        genomeReader = GenomeReader(TwoBitFile(self.twoBitFile))
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

    def getOptionArgs(self, allowUpdate):
        "create a vector of options and values for passing to another program"
        args = []
        if self.twoBitFile is not None:
            args.append("--genomeSeqs={}".format(self.twoBitFile))
        if self.mockTsv is not None:
            args.append("--mockGenomeSeqs={}".format(self.mockTsv))
        if self.updateMockReader and allowUpdate:
            args.append("--updateMockGenomeSeqs")
        if self.forceMockReader:
            args.append("--forceMockGenomeSeqs")
        return args

    @staticmethod
    def addCmdOptions(parser):
        """Add command options used to create a factory to the parser"""
        parser.add_argument('--genomeSeqs',
                            help="""Genome sequence twobit file to obtain splice sites""")
        parser.add_argument('--mockGenomeSeqs',
                            help="""Genome sequence TSV for mock reader for testing. This will be used if genomeSeqs is not specified or does not exist.""")
        parser.add_argument('--updateMockGenomeSeqs', action="store_true",
                            help="""update mock genome seqs.""")
        parser.add_argument('--forceMockGenomeSeqs', action="store_true",
                            help="""force using mock genome seqs, even if twobit exists""")

    @staticmethod
    def factoryFromCmdOptions(opts):
        """create a factory given options parse from command line"""
        # this will check sanity of options
        return GenomeReaderFactory(opts.genomeSeqs, opts.mockGenomeSeqs,
                                   opts.updateMockGenomeSeqs, opts.forceMockGenomeSeqs)
