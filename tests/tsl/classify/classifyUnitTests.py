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
from gencode_icedb.genome import GenomeReader
from gencode_icedb.tsl.evidFeatures import EvidenceFeatureMap, EvidencePslFactory
from gencode_icedb.tsl.annotFeatures import AnnotationGenePredFactory
from twobitreader import TwoBitFile
from pycbio.hgdata.hgLite import PslDbTable, GenePredDbTable
import sqlite3


twoBitHg19 = "/hive/data/genomes/hg19/hg19.2bit"
mockReaderHg19PslTsv = "mockReaderHg19Psl.tsv"
mockReaderHg19GpTsv = "mockReaderHg19Gp.tsv"

updateMockReader = False   # set this to update from a real run
forceMockReader = False  # set this to check mock data
debugResults = False   # print out results for updated expected

if updateMockReader or forceMockReader or debugResults:
    print("Warning: debug variables set", file=sys.stderr)
if updateMockReader and forceMockReader:
    raise Exception("makes no sense to have both updateMockReader and forceMockReader set")


class MockGenomeReader(object):
    """Fake GenomeReader when 2bit isn't there"""
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
    """Wrapper around GenomeReader that create mock data from twoBit """
    def __init__(self, genomeReader, mockDataTsv):
        self.genomeReader = genomeReader
        self.mockDataFh = open(mockDataTsv, "w")
        fileOps.prRowv(self.mockDataFh, "chrom", "start", "end", "strand", "seq")

    def get(self, chrom, start, end, strand=None):
        if strand is None:
            strand = '+'
        assert strand == '+'
        seq = self.genomeReader.get(chrom, start, end, strand)
        fileOps.prRowv(self.mockDataFh, chrom, start, end, strand, seq)
        self.mockDataFh.flush()
        return seq


class GenomeReaderFactory(object):
    """obtain the appropriate genome reader"""
    def __init__(self, mockTsv):
        self.mockTsv = mockTsv
        self.genomeReader = None

    def __getReal(self, testCase):
        genomeReader = GenomeReader(TwoBitFile(twoBitHg19))
        if updateMockReader:
            genomeReader = MockSeqWriter(genomeReader,
                                         testCase.getInputFile(self.mockTsv))
        return genomeReader

    def __getMock(self, testCase):
        if updateMockReader:
            raise Exception("updateMockReader is True, however twobit {} is not available".format(twoBitHg19))
        return MockGenomeReader(testCase.getInputFile(self.mockTsv))

    def obtain(self, testCase):
        if self.genomeReader is None:
            if os.path.exists(twoBitHg19) and not forceMockReader:
                self.genomeReader = self.__getReal(testCase)
            else:
                self.genomeReader = self.__getMock(testCase)
        return self.genomeReader


class EvidenceTests(TestCaseBase):
    genomeReaderFactory = GenomeReaderFactory(mockReaderHg19PslTsv)
    set1PslDbTbl = None

    def __obtainGenomeReader(self):
        return EvidenceTests.genomeReaderFactory.obtain(self)

    def __obtainSetPslDbTbl(self):
        if EvidenceTests.set1PslDbTbl is None:
            conn = sqlite3.connect(":memory:")
            EvidenceTests.set1PslDbTbl = PslDbTable(conn, "set1", create=True)
            EvidenceTests.set1PslDbTbl.loadPslFile(self.getInputFile("set1.ucsc-mrna.psl"))
        return EvidenceTests.set1PslDbTbl

    def __getSet1Psl(self, acc):
        pslDbTbl = self.__obtainSetPslDbTbl()
        psls = pslDbTbl.getByQName(acc)
        if len(psls) == 0:
            raise Exception("psl not found: {}".format(acc))
        return psls[0]

    def __pslToEvidTranscript(self, psl):
        factory = EvidencePslFactory(self.__obtainGenomeReader())
        return factory.fromPsl(psl)

    def __assertFeatures(self, feats, expectFeatsStr, expectFeatStrs):
        if debugResults:
            print('"{}"'.format(feats))
            for f in feats:
                print('"{}",'.format(f))
        self.assertEqual(str(feats), expectFeatsStr)
        self.assertEqual([str(f) for f in feats], expectFeatStrs)

    def testAF010310(self):
        feats = self.__pslToEvidTranscript(self.__getSet1Psl("AF010310.1"))
        self.__assertFeatures(feats, "t=chr22:18900294-18905926 -, q=AF010310.1:0=888 901",
                              ["exon 18900294-18900875 qIns=32 tIns=8",
                               "intron 18900875-18900950 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18900950-18901039 qIns=0 tIns=0",
                               "intron 18901039-18904402 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18904402-18904501 qIns=0 tIns=0",
                               "intron 18904501-18905828 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18905828-18905926 qIns=1 tIns=4"])

    def testX96484(self):
        feats = self.__pslToEvidTranscript(self.__getSet1Psl("X96484.1"))
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

    def testX96484NoSJ(self):
        factory = EvidencePslFactory(None)
        feats = factory.fromPsl(self.__getSet1Psl("X96484.1"))
        self.__assertFeatures(feats, "t=chr22:18893922-18899592 +, q=X96484.1:48=1067 1080",
                              ["exon 18893922-18893997 qIns=0 tIns=0",
                               "intron 18893997-18894077 qDel=0 sjBases=None",
                               "exon 18894077-18894238 qIns=1 tIns=0",
                               "intron 18894238-18897684 qDel=0 sjBases=None",
                               "exon 18897684-18897785 qIns=0 tIns=0",
                               "intron 18897785-18898400 qDel=0 sjBases=None",
                               "exon 18898400-18898541 qIns=0 tIns=0",
                               "intron 18898541-18899052 qDel=0 sjBases=None",
                               "exon 18899052-18899592 qIns=0 tIns=0"])

    def testRangeMap1(self):
        # range is set1: chr22:18632931-19279166
        pslDbTbl = self.__obtainSetPslDbTbl()
        evidFeatureMap = EvidenceFeatureMap.dbFactory(pslDbTbl.conn, pslDbTbl.table,
                                                      "chr22", 18958026, 19109719,
                                                      self.__obtainGenomeReader())
        self.assertEqual(len(evidFeatureMap.transcripts), 30)
        overFeats = list(evidFeatureMap.overlapping("chr22", 18958026, 18982141))
        self.assertEqual(len(overFeats), 12)


class AnnotationTests(TestCaseBase):
    genomeReaderFactory = GenomeReaderFactory(mockReaderHg19GpTsv)
    set1GpDbTbl = None

    def __obtainGenomeReader(self):
        return AnnotationTests.genomeReaderFactory.obtain(self)

    def __obtainSetGpDbTbl(self):
        if AnnotationTests.set1GpDbTbl is None:
            conn = sqlite3.connect(":memory:")
            AnnotationTests.set1GpDbTbl = GenePredDbTable(conn, "set1", create=True)
            AnnotationTests.set1GpDbTbl.loadGenePredFile(self.getInputFile("set1.gencodeCompV19.gp"))
        return AnnotationTests.set1GpDbTbl

    def __getSet1Gp(self, acc):
        gpDbTbl = self.__obtainSetGpDbTbl()
        gps = gpDbTbl.getByName(acc)
        if len(gps) == 0:
            raise Exception("gp not found: {}".format(acc))
        return gps[0]

    def __gpToEvidTranscript(self, gp):
        factory = AnnotationGenePredFactory(self.__obtainGenomeReader())
        return factory.fromGenePred(gp)

    def __assertFeatures(self, feats, expectFeatsStr, expectFeatStrs):
        if debugResults:
            print('"{}"'.format(feats))
            for f in feats:
                print('"{}",'.format(f))
        self.assertEqual(str(feats), expectFeatsStr)
        self.assertEqual([str(f) for f in feats], expectFeatStrs)

    def testENST00000215794(self):
        # + strand
        feats = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000215794.7"))
        self.__assertFeatures(feats, "t=chr22:18632665-18660164 +, q=ENST00000215794.7:0=2129 2129",
                              ["exon 18632665-18632989 qIns=0 tIns=0",
                               "intron 18632989-18640324 qDel=0 sjBases=GC...AG (spliceGC_AG)",
                               "exon 18640324-18640587 qIns=0 tIns=0",
                               "intron 18640587-18642938 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18642938-18643035 qIns=0 tIns=0",
                               "intron 18643035-18644556 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18644556-18644702 qIns=0 tIns=0",
                               "intron 18644702-18650021 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18650021-18650101 qIns=0 tIns=0",
                               "intron 18650101-18650656 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18650656-18650803 qIns=0 tIns=0",
                               "intron 18650803-18652610 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18652610-18652706 qIns=0 tIns=0",
                               "intron 18652706-18653519 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18653519-18653687 qIns=0 tIns=0",
                               "intron 18653687-18655916 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18655916-18656048 qIns=0 tIns=0",
                               "intron 18656048-18656559 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18656559-18656609 qIns=0 tIns=0",
                               "intron 18656609-18659538 qDel=0 sjBases=GT...AG (spliceGT_AG)",
                               "exon 18659538-18660164 qIns=0 tIns=0"])

    def testENST00000334029(self):
        # - strand
        feats = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000334029.2"))
        self.__assertFeatures(feats, "t=chr22:18900294-18923964 -, q=ENST00000334029.2:0=1985 1985",
                              ["exon 18900294-18900875 qIns=0 tIns=0",
                               "intron 18900875-18900950 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18900950-18901039 qIns=0 tIns=0",
                               "intron 18901039-18904402 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18904402-18904501 qIns=0 tIns=0",
                               "intron 18904501-18905828 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18905828-18906004 qIns=0 tIns=0",
                               "intron 18906004-18906963 qDel=0 sjBases=CT...GC (spliceGC_AG)",
                               "exon 18906963-18907110 qIns=0 tIns=0",
                               "intron 18907110-18907218 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18907218-18907311 qIns=0 tIns=0",
                               "intron 18907311-18908854 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18908854-18908936 qIns=0 tIns=0",
                               "intron 18908936-18909837 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18909837-18909917 qIns=0 tIns=0",
                               "intron 18909917-18910329 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18910329-18910446 qIns=0 tIns=0",
                               "intron 18910446-18910627 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18910627-18910692 qIns=0 tIns=0",
                               "intron 18910692-18912563 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18912563-18912713 qIns=0 tIns=0",
                               "intron 18912713-18913200 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18913200-18913235 qIns=0 tIns=0",
                               "intron 18913235-18918502 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18918502-18918711 qIns=0 tIns=0",
                               "intron 18918711-18923902 qDel=0 sjBases=CT...AC (spliceGT_AG)",
                               "exon 18923902-18923964 qIns=0 tIns=0"])

    def testENST00000334029NoSJ(self):
        factory = AnnotationGenePredFactory(None)
        feats = factory.fromGenePred(self.__getSet1Gp("ENST00000334029.2"))
        self.__assertFeatures(feats, "t=chr22:18900294-18923964 -, q=ENST00000334029.2:0=1985 1985",
                              ["exon 18900294-18900875 qIns=0 tIns=0",
                               "intron 18900875-18900950 qDel=0 sjBases=None",
                               "exon 18900950-18901039 qIns=0 tIns=0",
                               "intron 18901039-18904402 qDel=0 sjBases=None",
                               "exon 18904402-18904501 qIns=0 tIns=0",
                               "intron 18904501-18905828 qDel=0 sjBases=None",
                               "exon 18905828-18906004 qIns=0 tIns=0",
                               "intron 18906004-18906963 qDel=0 sjBases=None",
                               "exon 18906963-18907110 qIns=0 tIns=0",
                               "intron 18907110-18907218 qDel=0 sjBases=None",
                               "exon 18907218-18907311 qIns=0 tIns=0",
                               "intron 18907311-18908854 qDel=0 sjBases=None",
                               "exon 18908854-18908936 qIns=0 tIns=0",
                               "intron 18908936-18909837 qDel=0 sjBases=None",
                               "exon 18909837-18909917 qIns=0 tIns=0",
                               "intron 18909917-18910329 qDel=0 sjBases=None",
                               "exon 18910329-18910446 qIns=0 tIns=0",
                               "intron 18910446-18910627 qDel=0 sjBases=None",
                               "exon 18910627-18910692 qIns=0 tIns=0",
                               "intron 18910692-18912563 qDel=0 sjBases=None",
                               "exon 18912563-18912713 qIns=0 tIns=0",
                               "intron 18912713-18913200 qDel=0 sjBases=None",
                               "exon 18913200-18913235 qIns=0 tIns=0",
                               "intron 18913235-18918502 qDel=0 sjBases=None",
                               "exon 18918502-18918711 qIns=0 tIns=0",
                               "intron 18918711-18923902 qDel=0 sjBases=None",
                               "exon 18923902-18923964 qIns=0 tIns=0"])


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidenceTests))
    return ts


if __name__ == '__main__':
    unittest.main()
