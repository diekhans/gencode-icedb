from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from pycbio.hgdata.psl import Psl
from gencode_icedb.general.genome import GenomeReaderFactory
from gencode_icedb.general.gencodeDb import UcscGencodeReader
from gencode_icedb.tsl.evidenceDb import EvidenceSource, EvidenceReader
from gencode_icedb.tsl.supportDefs import EvidenceSupport
from gencode_icedb.tsl.supportClassify import compareMegWithEvidence, writeTsvHeaders, classifyGeneTranscripts


class EvidCompareTest(TestCaseBase):
    """low-level comparison tests"""
    UCSC_DB = "hg38"
    EVIDENCE_DB = "output/evid.db"
    GENCODE_DB = "output/gencode.db"

    @classmethod
    def setUpClass(cls):
        cls.genomeReader = GenomeReaderFactory.factoryFromUcscDb(cls.UCSC_DB).obtain()
        cls.gencodeReader = UcscGencodeReader(cls.GENCODE_DB, cls.genomeReader)
        cls.evidenceReader = EvidenceReader(cls.EVIDENCE_DB, cls.genomeReader)

    @classmethod
    def tearDownClass(cls):
        cls.genomeReader.close()
        cls.gencodeReader.close()
        cls.evidenceReader.close()

    def _getAnnot(self, transId):
        return self.gencodeReader.getByTranscriptIds(transId)[0]

    def _pslToTrans(self, psl):
        return self.evidenceReader.evidFactory.fromPsl(psl)

    def _evalAnnotTransEvid(self, annotTrans, evidSrc, evidNames=None, allowExtension=False):
        "low-level testing of specific cases"
        if evidNames is None:
            evidTranses = self.evidenceReader.genOverlapping(evidSrc, annotTrans.chrom.name, annotTrans.chrom.start, annotTrans.chrom.end, annotTrans.rna.strand)
        else:
            evidTranses = self.evidenceReader.genByNames(evidSrc, evidNames)
        for evidTrans in evidTranses:
            yield (annotTrans.rna.name, evidTrans.rna.name,
                   compareMegWithEvidence(annotTrans, evidTrans, allowExtension=allowExtension))

    def _classifyTest(self, annotTranses, noDiff=False):
        "TSL testing with TSV"
        outTslTsv = self.getOutputFile(".tsl.tsv")
        outDetailsTsv = self.getOutputFile(".details.tsv")
        with open(outTslTsv, 'w') as tslTsvFh, open(outDetailsTsv, 'w') as detailsTsvFh:
            writeTsvHeaders(tslTsvFh, detailsTsvFh)
            classifyGeneTranscripts(self.evidenceReader, annotTranses, tslTsvFh, detailsTsvFh)
        if not noDiff:
            self.diffFiles(self.getExpectedFile(".tsl.tsv"), outTslTsv)
            self.diffFiles(self.getExpectedFile(".details.tsv"), outDetailsTsv)

    def _classifyGeneTest(self, geneId):
        geneTranses = self.gencodeReader.getByGeneId(geneId)
        if len(geneTranses) == 0:
            raise Exception("no transcripts found for {}".format(geneId))
        self._classifyTest(geneTranses)

    def _classifyTransTest(self, transId, noDiff=False):
        transes = self.gencodeReader.getByTranscriptId(transId)
        if len(transes) == 0:
            raise Exception("no transcripts found for {}".format(transId))
        self._classifyTest(transes, noDiff)

    def testGAB4(self):
        self._classifyGeneTest("ENSG00000215568.7")

    def testBCR(self):
        self._classifyGeneTest("ENSG00000186716.20")

    def testIL17RA(self):
        self._classifyGeneTest("ENSG00000177663.13")

    def testSHOX(self):
        # PAR gene
        self._classifyGeneTest("ENSG00000185960.13")

    def testExtendWithTwoExonsOverInitial(self):
        # EST AA227241.1 has two 5' exons overlapping 5' exon, caused failure with allowExtension
        annotName = "ENST00000489867.1"
        evidName = "AA227241.1"
        results = tuple(self._evalAnnotTransEvid(self._getAnnot(annotName),
                                                 EvidenceSource.UCSC_EST,
                                                 evidNames=evidName, allowExtension=True))
        self.assertEqual(results,
                         (('ENST00000489867.1', 'AA227241.1', EvidenceSupport.feat_mismatch),))

    def _rawPslCmpr(self, annotTrans, evidRawPsl, allowExtension):
        evidTrans = self._pslToTrans(Psl(evidRawPsl))
        return compareMegWithEvidence(annotTrans, evidTrans, allowExtension=allowExtension)

    def testExact(self):
        # test allowExtension with fake psl that exactly matches
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["651", "0", "0", "0", "0", "0", "3", "14865", "+", "ENST00000477874.1_aln", "651", "0", "651", "chr22", "50818468", "17084953", "17100469", "4", "276,147,113,115,", "0,276,423,536,", "17084953,17097796,17098774,17100354,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport, EvidenceSupport.good)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport, EvidenceSupport.good)

    def testExtend5(self):
        # test allowExtension with fake psl on extended on 5 side
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["751", "0", "0", "0", "0", "0", "4", "15218", "+", "ENST00000477874.1_aln", "751", "0", "751", "chr22", "50818468", "17084500", "17100469", "5", "100,276,147,113,115,", "0,100,376,523,636,", "17084500,17084953,17097796,17098774,17100354,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport, EvidenceSupport.good)

    def testExtend5Inexact(self):
        # test allowExtension with fake psl on extended on 5 side and 5'
        # annotated exon being shorter than evidence, which is not allowed
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["704", "0", "0", "0", "0", "0", "4", "15265", "+", "ENST00000477874.1_aln", "704", "0", "704", "chr22", "50818468", "17084500", "17100469", "5", "100,229,147,113,115,", "0,100,329,476,589,", "17084500,17085000,17097796,17098774,17100354,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport, EvidenceSupport.feat_mismatch)

    def testExtend3(self):
        # test allowExtension with fake psl on extended on 3' side
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["741", "0", "0", "0", "0", "0", "4", "14896", "+", "ENST00000477874.1_aln", "741", "0", "741", "chr22", "50818468", "17084953", "17100590", "5", "276,147,113,115,90,", "0,276,423,536,651,", "17084953,17097796,17098774,17100354,17100500,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport, EvidenceSupport.good)

    def testExtend3Inexact(self):
        # test allowExtension with fake psl on extended on 3' side
        # and 3' annotated exon being shorter than evidence, which is
        # not allowed
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["682", "0", "0", "0", "0", "0", "4", "15165", "+", "ENST00000477874.1_aln", "682", "0", "682", "chr22", "50818468", "17084953", "17100800", "5", "276,147,113,46,100,", "0,276,423,536,582,", "17084953,17097796,17098774,17100354,17100700,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport, EvidenceSupport.feat_mismatch)

    def testExtend3Regression1(self):
        # Important test: Final two annotation exons mostly covered by an a
        # single exon of BC071657.1, however the overlapping exons ends
        # differently.  This was incorrectly called as an extension.  Both
        # happen to have the same numbner of features.  This is with the
        # Ensembl alignment, UCSC aligned differently.
        annotTrans = self._getAnnot("ENST00000315091.7")
        evidRawPsl = ["2763", "0", "0", "0", "0", "0", "8", "107096", "+", "BC071657.1", "2833", "0", "2763", "chr1", "248956422", "11012670", "11122529", "9", "73,250,164,141,171,1919,17,9,19,", "0,73,323,487,628,799,2718,2735,2744,", "11012670,11013715,11016843,11018732,11020428,11022123,11034700,11034718,11122510,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport, EvidenceSupport.feat_mismatch)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidCompareTest))
    return ts


if __name__ == '__main__':
    unittest.main()
