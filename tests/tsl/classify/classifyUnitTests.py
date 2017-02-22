from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from pycbio.tsv import TsvReader
from pycbio.sys import fileOps
from gencode_icedb.tsl.evidFeatures import SeqReader, EvidFeatures, EvidFeaturesMap
from twobitreader import TwoBitFile
from pycbio.hgdata.hgLite import PslDbTable
import sqlite3


twoBitHg19 = "/hive/data/genomes/hg19/hg19.2bit"
mockReaderHg19Tsv = "mockReaderHg19.tsv"
updateMockReader = False   # set this to update from a real run

debugResults = False   # print out results for updated expected


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
    genomeReader = None
    set1PslDbTbl = None

    def __getRealSeqReader(self):
        genomeReader = SeqReader(TwoBitFile(twoBitHg19))
        if updateMockReader:
            genomeReader = MockSeqWriter(genomeReader,
                                         self.getInputFile(mockReaderHg19Tsv))
        return genomeReader

    def __getMockSeqReader(self):
        if updateMockReader:
            raise Exception("updateMockReader is True, however twobit {} is not available".format(twoBitHg19))
        return MockSeqReader(self.getInputFile(mockReaderHg19Tsv))

    def __obtainGenomeReader(self):
        if EvidenceTests.genomeReader is None:
            if os.path.exists(twoBitHg19):
                EvidenceTests.genomeReader = self.__getRealSeqReader()
            else:
                EvidenceTests.genomeReader = self.__getMockSeqReader()
        return EvidenceTests.genomeReader

    def __obtainSetPslDbTbl(self):
        if EvidenceTests.set1PslDbTbl is None:
            conn = sqlite3.connect(":memory")
            EvidenceTests.set1PslDbTbl = PslDbTable(conn, "set1", create=True)
            EvidenceTests.set1PslDbTbl.loadPslFile(self.getInputFile("set1.ucsc-mrna.psl"))
        return EvidenceTests.set1PslDbTbl

    def __getSet1Psl(self, acc):
        pslDbTbl = self.__obtainSetPslDbTbl()
        psls = pslDbTbl.getByQName(acc)
        if len(psls) == 0:
            raise Exception("psl not found: {}".format(acc))
        return psls[0]

    def __assertFeatures(self, feats, expectFeatsStr, expectFeatStrs):
        if debugResults:
            print('"{}"'.format(feats))
            for f in feats:
                print('"{}",'.format(f))
        self.assertEqual(str(feats), expectFeatsStr)
        self.assertEqual([str(f) for f in feats], expectFeatStrs)

    def testAF010310(self):
        psl = self.__getSet1Psl("AF010310.1")
        feats = EvidFeatures(psl, self.__obtainGenomeReader())
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
        feats = EvidFeatures(psl, self.__obtainGenomeReader())
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

    def testRangeMap1(self):
        # range is set1: chr22:18632931-19279166
        pslDbTbl = self.__obtainSetPslDbTbl()
        evidFeaturesMap = EvidFeaturesMap.dbFactory(pslDbTbl.conn, pslDbTbl.table,
                                                    "chr22", 18958026, 19109719,
                                                    self.__obtainGenomeReader())
        self.assertEqual(len(evidFeaturesMap.featuresList), 30)
        overFeats = list(evidFeaturesMap.overlapping("chr22", 18958026, 18982141))
        self.assertEqual(len(overFeats), 12)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidenceTests))
    return ts


if __name__ == '__main__':
    unittest.main()
