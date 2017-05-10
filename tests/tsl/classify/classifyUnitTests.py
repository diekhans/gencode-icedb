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


class FeatureTestBase(TestCaseBase):
    def _assertFeatures(self, trans, expectTransStr, expectFeatStrs):
        if debugResults:
            print(self.id())
            print('  "{}"'.format(trans))
            for f in trans.features:
                print('    "{}",'.format(f))
            print("")
        self.assertEqual(str(trans), expectTransStr)
        self.assertEqual([str(f) for f in trans.features], expectFeatStrs)


class EvidenceTests(FeatureTestBase):
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

    def testAF010310(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("AF010310.1"))
        self._assertFeatures(trans, "t=chr22:18900294-18905926/+, rna=AF010310.1:0-888/- 901",
                             ["exon 18900294-18900875 rna=13-618 rIns=32 cIns=8",
                              "intron 18900875-18900950 rna=618-618 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18900950-18901039 rna=618-707 rIns=0 cIns=0",
                              "intron 18901039-18904402 rna=707-707 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18904402-18904501 rna=707-806 rIns=0 cIns=0",
                              "intron 18904501-18905828 rna=806-806 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18905828-18905926 rna=806-901 rIns=1 cIns=4"])

    def testX96484(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("X96484.1"))
        self._assertFeatures(trans, "t=chr22:18893922-18899592/+, rna=X96484.1:48-1067/+ 1080",
                             ["exon 18893922-18893997 rna=48-123 rIns=0 cIns=0",
                              "intron 18893997-18894077 rna=123-123 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18894077-18894238 rna=123-285 rIns=1 cIns=0",
                              "intron 18894238-18897684 rna=285-285 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18897684-18897785 rna=285-386 rIns=0 cIns=0",
                              "intron 18897785-18898400 rna=386-386 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18898400-18898541 rna=386-527 rIns=0 cIns=0",
                              "intron 18898541-18899052 rna=527-527 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18899052-18899592 rna=527-1067 rIns=0 cIns=0"])

    def testX96484NoSJ(self):
        factory = EvidencePslFactory(None)
        trans = factory.fromPsl(self.__getSet1Psl("X96484.1"))
        self._assertFeatures(trans, "t=chr22:18893922-18899592/+, rna=X96484.1:48-1067/+ 1080",
                             ["exon 18893922-18893997 rna=48-123 rIns=0 cIns=0",
                              "intron 18893997-18894077 rna=123-123 rDel=0 sjBases=None",
                              "exon 18894077-18894238 rna=123-285 rIns=1 cIns=0",
                              "intron 18894238-18897684 rna=285-285 rDel=0 sjBases=None",
                              "exon 18897684-18897785 rna=285-386 rIns=0 cIns=0",
                              "intron 18897785-18898400 rna=386-386 rDel=0 sjBases=None",
                              "exon 18898400-18898541 rna=386-527 rIns=0 cIns=0",
                              "intron 18898541-18899052 rna=527-527 rDel=0 sjBases=None",
                              "exon 18899052-18899592 rna=527-1067 rIns=0 cIns=0"])

    def testAF010310Rc(self):
        trans = self.__pslToEvidTranscript(self.__getSet1Psl("AF010310.1"))
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self._assertFeatures(rcTrans,   "t=chr22:32398640-32404272/-, rna=AF010310.1:13-901/+ 901",
                             ["exon 32398640-32398738 rna=0-95 rIns=1 cIns=4",
                              "intron 32398738-32400065 rna=95-95 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32400065-32400164 rna=95-194 rIns=0 cIns=0",
                              "intron 32400164-32403527 rna=194-194 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32403527-32403616 rna=194-283 rIns=0 cIns=0",
                              "intron 32403616-32403691 rna=283-283 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32403691-32404272 rna=283-888 rIns=32 cIns=8"])

    def testRangeMap1(self):
        # range is set1: chr22:18632931-19279166
        pslDbTbl = self.__obtainSetPslDbTbl()
        evidFeatureMap = EvidenceFeatureMap.dbFactory(pslDbTbl.conn, pslDbTbl.table,
                                                      "chr22", 18958026, 19109719,
                                                      self.__obtainGenomeReader())
        self.assertEqual(len(evidFeatureMap.transcripts), 30)
        overFeats = list(evidFeatureMap.overlapping("chr22", 18958026, 18982141))
        self.assertEqual(len(overFeats), 12)


class AnnotationTests(FeatureTestBase):
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

    def testENST00000215794(self):
        # + strand
        trans = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000215794.7"))
        self._assertFeatures(trans, "t=chr22:18632665-18660164/+, rna=ENST00000215794.7:0-2129/+ 2129",
                             ["exon 18632665-18632989 rna=0-324 rIns=0 cIns=0",
                              "intron 18632989-18640324 rna=324-324 rDel=0 sjBases=GC...AG (spliceGC_AG)",
                              "exon 18640324-18640587 rna=324-587 rIns=0 cIns=0",
                              "intron 18640587-18642938 rna=587-587 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18642938-18643035 rna=587-684 rIns=0 cIns=0",
                              "intron 18643035-18644556 rna=684-684 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18644556-18644702 rna=684-830 rIns=0 cIns=0",
                              "intron 18644702-18650021 rna=830-830 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18650021-18650101 rna=830-910 rIns=0 cIns=0",
                              "intron 18650101-18650656 rna=910-910 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18650656-18650803 rna=910-1057 rIns=0 cIns=0",
                              "intron 18650803-18652610 rna=1057-1057 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18652610-18652706 rna=1057-1153 rIns=0 cIns=0",
                              "intron 18652706-18653519 rna=1153-1153 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18653519-18653687 rna=1153-1321 rIns=0 cIns=0",
                              "intron 18653687-18655916 rna=1321-1321 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18655916-18656048 rna=1321-1453 rIns=0 cIns=0",
                              "intron 18656048-18656559 rna=1453-1453 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18656559-18656609 rna=1453-1503 rIns=0 cIns=0",
                              "intron 18656609-18659538 rna=1503-1503 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 18659538-18660164 rna=1503-2129 rIns=0 cIns=0"])

    def testENST00000334029(self):
        # - strand
        trans = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000334029.2"))
        self._assertFeatures(trans, "t=chr22:18900294-18923964/+, rna=ENST00000334029.2:0-1985/- 1985",
                             ["exon 18900294-18900875 rna=0-581 rIns=0 cIns=0",
                              "intron 18900875-18900950 rna=581-581 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18900950-18901039 rna=581-670 rIns=0 cIns=0",
                              "intron 18901039-18904402 rna=670-670 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18904402-18904501 rna=670-769 rIns=0 cIns=0",
                              "intron 18904501-18905828 rna=769-769 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18905828-18906004 rna=769-945 rIns=0 cIns=0",
                              "intron 18906004-18906963 rna=945-945 rDel=0 sjBases=CT...GC (spliceGC_AG)",
                              "exon 18906963-18907110 rna=945-1092 rIns=0 cIns=0",
                              "intron 18907110-18907218 rna=1092-1092 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18907218-18907311 rna=1092-1185 rIns=0 cIns=0",
                              "intron 18907311-18908854 rna=1185-1185 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18908854-18908936 rna=1185-1267 rIns=0 cIns=0",
                              "intron 18908936-18909837 rna=1267-1267 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18909837-18909917 rna=1267-1347 rIns=0 cIns=0",
                              "intron 18909917-18910329 rna=1347-1347 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18910329-18910446 rna=1347-1464 rIns=0 cIns=0",
                              "intron 18910446-18910627 rna=1464-1464 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18910627-18910692 rna=1464-1529 rIns=0 cIns=0",
                              "intron 18910692-18912563 rna=1529-1529 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18912563-18912713 rna=1529-1679 rIns=0 cIns=0",
                              "intron 18912713-18913200 rna=1679-1679 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18913200-18913235 rna=1679-1714 rIns=0 cIns=0",
                              "intron 18913235-18918502 rna=1714-1714 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18918502-18918711 rna=1714-1923 rIns=0 cIns=0",
                              "intron 18918711-18923902 rna=1923-1923 rDel=0 sjBases=CT...AC (spliceGT_AG)",
                              "exon 18923902-18923964 rna=1923-1985 rIns=0 cIns=0"])

    def testENST00000334029NoSJ(self):
        factory = AnnotationGenePredFactory(None)
        trans = factory.fromGenePred(self.__getSet1Gp("ENST00000334029.2"))
        self._assertFeatures(trans, "t=chr22:18900294-18923964/+, rna=ENST00000334029.2:0-1985/- 1985",
                             ["exon 18900294-18900875 rna=0-581 rIns=0 cIns=0",
                              "intron 18900875-18900950 rna=581-581 rDel=0 sjBases=None",
                              "exon 18900950-18901039 rna=581-670 rIns=0 cIns=0",
                              "intron 18901039-18904402 rna=670-670 rDel=0 sjBases=None",
                              "exon 18904402-18904501 rna=670-769 rIns=0 cIns=0",
                              "intron 18904501-18905828 rna=769-769 rDel=0 sjBases=None",
                              "exon 18905828-18906004 rna=769-945 rIns=0 cIns=0",
                              "intron 18906004-18906963 rna=945-945 rDel=0 sjBases=None",
                              "exon 18906963-18907110 rna=945-1092 rIns=0 cIns=0",
                              "intron 18907110-18907218 rna=1092-1092 rDel=0 sjBases=None",
                              "exon 18907218-18907311 rna=1092-1185 rIns=0 cIns=0",
                              "intron 18907311-18908854 rna=1185-1185 rDel=0 sjBases=None",
                              "exon 18908854-18908936 rna=1185-1267 rIns=0 cIns=0",
                              "intron 18908936-18909837 rna=1267-1267 rDel=0 sjBases=None",
                              "exon 18909837-18909917 rna=1267-1347 rIns=0 cIns=0",
                              "intron 18909917-18910329 rna=1347-1347 rDel=0 sjBases=None",
                              "exon 18910329-18910446 rna=1347-1464 rIns=0 cIns=0",
                              "intron 18910446-18910627 rna=1464-1464 rDel=0 sjBases=None",
                              "exon 18910627-18910692 rna=1464-1529 rIns=0 cIns=0",
                              "intron 18910692-18912563 rna=1529-1529 rDel=0 sjBases=None",
                              "exon 18912563-18912713 rna=1529-1679 rIns=0 cIns=0",
                              "intron 18912713-18913200 rna=1679-1679 rDel=0 sjBases=None",
                              "exon 18913200-18913235 rna=1679-1714 rIns=0 cIns=0",
                              "intron 18913235-18918502 rna=1714-1714 rDel=0 sjBases=None",
                              "exon 18918502-18918711 rna=1714-1923 rIns=0 cIns=0",
                              "intron 18918711-18923902 rna=1923-1923 rDel=0 sjBases=None",
                              "exon 18923902-18923964 rna=1923-1985 rIns=0 cIns=0"])


    def testENST00000334029Rc(self):
        # - strand
        trans = self.__gpToEvidTranscript(self.__getSet1Gp("ENST00000334029.2"))
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self._assertFeatures(rcTrans, "t=chr22:32380602-32404272/-, rna=ENST00000334029.2:0-1985/+ 1985",
                             ["exon 32380602-32380664 rna=0-62 rIns=0 cIns=0",
                              "intron 32380664-32385855 rna=62-62 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32385855-32386064 rna=62-271 rIns=0 cIns=0",
                              "intron 32386064-32391331 rna=271-271 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32391331-32391366 rna=271-306 rIns=0 cIns=0",
                              "intron 32391366-32391853 rna=306-306 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32391853-32392003 rna=306-456 rIns=0 cIns=0",
                              "intron 32392003-32393874 rna=456-456 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32393874-32393939 rna=456-521 rIns=0 cIns=0",
                              "intron 32393939-32394120 rna=521-521 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32394120-32394237 rna=521-638 rIns=0 cIns=0",
                              "intron 32394237-32394649 rna=638-638 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32394649-32394729 rna=638-718 rIns=0 cIns=0",
                              "intron 32394729-32395630 rna=718-718 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32395630-32395712 rna=718-800 rIns=0 cIns=0",
                              "intron 32395712-32397255 rna=800-800 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32397255-32397348 rna=800-893 rIns=0 cIns=0",
                              "intron 32397348-32397456 rna=893-893 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32397456-32397603 rna=893-1040 rIns=0 cIns=0",
                              "intron 32397603-32398562 rna=1040-1040 rDel=0 sjBases=GC...AG (spliceGC_AG)",
                              "exon 32398562-32398738 rna=1040-1216 rIns=0 cIns=0",
                              "intron 32398738-32400065 rna=1216-1216 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32400065-32400164 rna=1216-1315 rIns=0 cIns=0",
                              "intron 32400164-32403527 rna=1315-1315 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32403527-32403616 rna=1315-1404 rIns=0 cIns=0",
                              "intron 32403616-32403691 rna=1404-1404 rDel=0 sjBases=GT...AG (spliceGT_AG)",
                              "exon 32403691-32404272 rna=1404-1985 rIns=0 cIns=0"])

    def testENST00000334029NoSJRc(self):
        # no splice sites, just sizes
        factory = AnnotationGenePredFactory(chromSizeFunc=self.__obtainGenomeReader().getChromSize)
        trans = factory.fromGenePred(self.__getSet1Gp("ENST00000334029.2"))
        rcTrans = trans.reverseComplement()
        self.assertEqual(len(rcTrans.features), len(trans.features))
        self._assertFeatures(rcTrans, "t=chr22:32380602-32404272/-, rna=ENST00000334029.2:0-1985/+ 1985",
                             ["exon 32380602-32380664 rna=0-62 rIns=0 cIns=0",
                              "intron 32380664-32385855 rna=62-62 rDel=0 sjBases=None",
                              "exon 32385855-32386064 rna=62-271 rIns=0 cIns=0",
                              "intron 32386064-32391331 rna=271-271 rDel=0 sjBases=None",
                              "exon 32391331-32391366 rna=271-306 rIns=0 cIns=0",
                              "intron 32391366-32391853 rna=306-306 rDel=0 sjBases=None",
                              "exon 32391853-32392003 rna=306-456 rIns=0 cIns=0",
                              "intron 32392003-32393874 rna=456-456 rDel=0 sjBases=None",
                              "exon 32393874-32393939 rna=456-521 rIns=0 cIns=0",
                              "intron 32393939-32394120 rna=521-521 rDel=0 sjBases=None",
                              "exon 32394120-32394237 rna=521-638 rIns=0 cIns=0",
                              "intron 32394237-32394649 rna=638-638 rDel=0 sjBases=None",
                              "exon 32394649-32394729 rna=638-718 rIns=0 cIns=0",
                              "intron 32394729-32395630 rna=718-718 rDel=0 sjBases=None",
                              "exon 32395630-32395712 rna=718-800 rIns=0 cIns=0",
                              "intron 32395712-32397255 rna=800-800 rDel=0 sjBases=None",
                              "exon 32397255-32397348 rna=800-893 rIns=0 cIns=0",
                              "intron 32397348-32397456 rna=893-893 rDel=0 sjBases=None",
                              "exon 32397456-32397603 rna=893-1040 rIns=0 cIns=0",
                              "intron 32397603-32398562 rna=1040-1040 rDel=0 sjBases=None",
                              "exon 32398562-32398738 rna=1040-1216 rIns=0 cIns=0",
                              "intron 32398738-32400065 rna=1216-1216 rDel=0 sjBases=None",
                              "exon 32400065-32400164 rna=1216-1315 rIns=0 cIns=0",
                              "intron 32400164-32403527 rna=1315-1315 rDel=0 sjBases=None",
                              "exon 32403527-32403616 rna=1315-1404 rIns=0 cIns=0",
                              "intron 32403616-32403691 rna=1404-1404 rDel=0 sjBases=None",
                              "exon 32403691-32404272 rna=1404-1985 rIns=0 cIns=0"])
def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidenceTests))
    return ts


if __name__ == '__main__':
    unittest.main()
