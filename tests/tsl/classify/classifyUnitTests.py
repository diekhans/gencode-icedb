from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from pycbio.hgdata.psl import PslTbl
from pycbio.tsv import TsvReader
from pycbio.sys import fileOps
from gencode_icedb.tsl.evidFeatures import SeqReader, EvidFeatures
from twobitreader import TwoBitFile


twoBitHg19 = "/hive/data/genomes/hg19/hg19.2bit"
mockReaderHg19Tsv = "mockReaderHg19.tsv"
updateMockReader = False   # set this to update from a real run

debugResults = True   # print out results for updated expected

class MockSeqReader(object):
    """Fake SeqReader when 2bit isn't there"""
    def __init__(self, mockDataTsv):
        self.mockData = dict()
        for row in TsvReader(mockDataTsv,
                             typeMap={"start": int,
                                      "end": int}):
            self.mockData[(row.chrom, row.start, row.end, row.strand)] = row.seq

    def get(self, chrom, start, end, strand=None):
        if strand is None:
            strand = '+'
        return self.mockData[(chrom, start, end, strand)]


class MockSeqWriter(object):
    """Wrapper around SeqReader that create mock data from twoBit """
    def __init__(self, seqReader, mockDataTsv):
        self.seqReader = seqReader
        self.mockDataFh = open(mockDataTsv, "w")
        fileOps.prRowv(self.mockDataFh, "chrom", "start", "end", "strand", "seq")

    def get(self, chrom, start, end, strand=None):
        if strand is None:
            strand = '+'
        seq = self.seqReader.get(chrom, start, end, strand)
        fileOps.prRowv(self.mockDataFh, chrom, start, end, strand, seq)
        return seq


class EvidenceTests(TestCaseBase):
    seqReader = None
    set1Psls = None

    def __getRealSeqReader(self):
        seqReader = SeqReader(TwoBitFile(twoBitHg19))
        if updateMockReader:
            seqReader = MockSeqWriter(seqReader,
                                      self.getInputFile(mockReaderHg19Tsv))
        return seqReader

    def __getMockSeqReader(self):
        if updateMockReader:
            raise Exception("updateMockReader is True, however twobit {} is not available".format(twoBitHg19))
        return MockSeqReader(self.getInputFile(mockReaderHg19Tsv))

    def __obtainSeqReader(self):
        if EvidenceTests.seqReader is None:
            if os.path.exists(twoBitHg19):
                EvidenceTests.seqReader = self.__getRealSeqReader()
            else:
                EvidenceTests.seqReader = self.__getMockSeqReader()
        return EvidenceTests.seqReader

    def __requireSet1Psls(self):
        if EvidenceTests.set1Psls is None:
            EvidenceTests.set1Psls = PslTbl(self.getInputFile("set1.ucsc-mrna.psl"),
                                            qNameIdx=True)

    def __getSet1Psl(self, acc):
        self.__requireSet1Psls()
        return EvidenceTests.set1Psls.qNameMap[acc][0]

    def __assertFeatures(self, feats, expectFeatsStr, expectFeatStrs):
        if debugResults:
            print('"{}"'.format(feats))
            for f in feats:
                print('"{}",'.format(f))
        self.assertEqual(str(feats), expectFeatsStr)
        self.assertEqual([str(f) for f in feats], expectFeatStrs)

    def testAF010310(self):
        psl = self.__getSet1Psl("AF010310.1")
        feats = EvidFeatures(psl, self.__obtainSeqReader())
        self.__assertFeatures(feats, "t=chr22:18900294-18905926 -, q=AF010310.1:0=888 901",
                              ["exon 18900294-18900875 qIns=32 tIns=8",
                               "intron 18900875-18900950 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18900950-18901039 qIns=0 tIns=0",
                               "intron 18901039-18904402 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18904402-18904501 qIns=0 tIns=0",
                               "intron 18904501-18905828 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18905828-18905926 qIns=1 tIns=4"])

    def testX96484(self):
        psl = self.__getSet1Psl("X96484.1")
        feats = EvidFeatures(psl, self.__obtainSeqReader())
        self.__assertFeatures(feats, "t=chr22:18893922-18899592 +, q=X96484.1:48=1067 1080",
                              ["exon 18893922-18893997 qIns=0 tIns=0",
                               "intron 18893997-18894077 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18894077-18894238 qIns=1 tIns=0",
                               "intron 18894238-18897684 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18897684-18897785 qIns=0 tIns=0",
                               "intron 18897785-18898400 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18898400-18898541 qIns=0 tIns=0",
                               "intron 18898541-18899052 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18899052-18899592 qIns=0 tIns=0"])

    def testSet1(self):
        # just run through to see if they all convert
        self.__requireSet1Psls()
        cnt = 0
        for psl in self.set1Psls:
            EvidFeatures(psl, self.__obtainSeqReader())
            cnt += 1
        self.assertEqual(cnt, 81)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidenceTests))
    return ts


if __name__ == '__main__':
    unittest.main()
