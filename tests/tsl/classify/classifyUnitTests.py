import sys
import os
if __name__ == '__main__':
    rootDir = "../../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from pycbio.hgdata.psl import Psl
from gencode_icedb.general.genome import GenomeReader
from gencode_icedb.general.ucscGencodeSource import UcscGencodeReader
from gencode_icedb.general.geneAnnot import geneAnnotGroup
from gencode_icedb.general.evidFeatures import EvidencePslFactory
from gencode_icedb.tsl.evidenceDataDb import evidenceAlignsReaderFactory
from gencode_icedb.tsl.supportDefs import EvidenceType, EvidenceSupport
from gencode_icedb.tsl.supportEval import tightExonPolymorphicSizeLimit, tightExonPolymorphicFactionLimit, EvidenceQualityEval, MegSupportEvaluator, FullLengthSupportEvaluator
from gencode_icedb.tsl.supportEvalDb import SupportEvidEvalResult


genbankUuids = {
    EvidenceType.RNA: "02c995e3-372c-4cde-b216-5d3376c51988",
    EvidenceType.EST: "fcfc5121-7e07-42f6-93a2-331a513eeb2c",
}


class EvidCompareTest(TestCaseBase):
    """low-level comparison tests"""
    UCSC_DB = "hg38"
    EVIDENCE_DB_DIR = "output/db/evidDb"
    GENCODE_DB = "output/db/gencode.db"

    @classmethod
    def setUpClass(cls):
        cls.genomeReader = GenomeReader.getFromUcscDbName(cls.UCSC_DB)
        cls.gencodeReader = UcscGencodeReader(cls.GENCODE_DB, cls.genomeReader)
        cls.evidenceFactory = EvidencePslFactory(cls.genomeReader)
        cls.qualEval = EvidenceQualityEval(tightExonPolymorphicSizeLimit, tightExonPolymorphicFactionLimit)
        cls.evidenceReaders = {}
        cls.evaluators = {}
        for evidType in (EvidenceType.RNA, EvidenceType.EST):
            evidSetUuid = genbankUuids[evidType]
            cls.evidenceReaders[evidType] = evidenceAlignsReaderFactory(evidSetUuid, os.path.join(cls.EVIDENCE_DB_DIR, str(evidType) + ".psl.gz"), cls.genomeReader)
            for allowExtension in (True, False):
                cls.evaluators[(evidType, allowExtension)] = MegSupportEvaluator(evidSetUuid, cls.qualEval, allowExtension=allowExtension)

    @classmethod
    def tearDownClass(cls):
        cls.genomeReader.close()
        cls.gencodeReader.close()
        for rdr in cls.evidenceReaders.values():
            rdr.close()

    @classmethod
    def _getEvaluator(cls, evidType, allowExtension):
        return cls.evaluators[(evidType, allowExtension)]

    def _getAnnot(self, transId):
        annots = self.gencodeReader.getByTranscriptIds(transId)
        if len(annots) == 0:
            raise Exception("annotation {} not found".format(transId))
        return annots[0]

    def _pslToTrans(self, psl):
        return self.evidenceFactory.fromPsl(psl)

    def _evalAnnotTransEvid(self, annotTrans, evidType, evidNames=None, allowExtension=False):
        "low-level testing of specific cases"
        evaluator = self._getEvaluator(evidType, allowExtension)
        evidReader = self.evidenceReaders[evidType]
        try:
            evidReader.setNameSubset(evidNames)  # could be None
            evidTranses = evidReader.genOverlapping(annotTrans.chrom, annotTrans.transcriptionStrand)
            for evidTrans in evidTranses:
                yield evaluator.compare(annotTrans, evidTrans)
        finally:
            # must finished reading for resetting, as this is a generator
            evidReader.setNameSubset(None)

    def _evalTest(self, geneAnnots, noDiff=False):
        "Support evaluation testing with TSV"
        outSupportTsv = self.getOutputFile(".support.tsv")
        outDetailsTsv = self.getOutputFile(".details.tsv")
        evaluatorRna = FullLengthSupportEvaluator(self.evidenceReaders[EvidenceType.RNA], self.qualEval)
        evaluatorEst = FullLengthSupportEvaluator(self.evidenceReaders[EvidenceType.EST], self.qualEval)
        with open(outSupportTsv, 'w') as evalTsvFh, open(outDetailsTsv, 'w') as detailsTsvFh:
            evaluatorRna.writeTsvHeaders(evalTsvFh, detailsTsvFh)
            for geneAnnot in geneAnnots:
                evaluatorRna.evaluateGeneTranscripts(geneAnnot, evalTsvFh, detailsTsvFh)
                evaluatorEst.evaluateGeneTranscripts(geneAnnot, evalTsvFh, detailsTsvFh)
        if not noDiff:
            self.diffFiles(self.getExpectedFile(".support.tsv"), outSupportTsv)
            self.diffFiles(self.getExpectedFile(".details.tsv"), outDetailsTsv)

    def _evalGeneTest(self, geneId):
        transAnnots = self.gencodeReader.getByGeneId(geneId)
        if len(transAnnots) == 0:
            raise Exception("no transcripts found for {}".format(geneId))
        self._evalTest(geneAnnotGroup(transAnnots))

    def _evalTransTest(self, transId, noDiff=False):
        transAnnots = self.gencodeReader.getByTranscriptId(transId)
        if len(transAnnots) == 0:
            raise Exception("no transcripts found for {}".format(transId))
        self._evalTest(geneAnnotGroup(transAnnots), noDiff)

    def testGAB4(self):
        self._evalGeneTest("ENSG00000215568.8")

    def testBCR(self):
        self._evalGeneTest("ENSG00000186716.20")

    def testIL17RA(self):
        self._evalGeneTest("ENSG00000177663.13")

    def testSHOX(self):
        # PAR gene
        self._evalGeneTest("ENSG00000185960.14")

    def testExtendWithTwoExonsOverInitial(self):
        # EST AA227241.1 has two 5' exons overlapping 5' exon, caused failure with allowExtension
        annotName = "ENST00000489867.2"
        evidName = "AA227241.1"
        results = tuple(self._evalAnnotTransEvid(self._getAnnot(annotName), EvidenceType.EST,
                                                 evidNames=evidName, allowExtension=True))
        self.assertEqual((SupportEvidEvalResult('ENST00000489867.2', genbankUuids[EvidenceType.EST], 'AA227241.1', EvidenceSupport.feat_count_mismatch),),
                         results)

    def _rawPslCmpr(self, annotTrans, evidRawPsl, allowExtension):
        evidTrans = self._pslToTrans(Psl.fromRow(evidRawPsl))
        evaluator = self._getEvaluator(EvidenceType.EST, allowExtension)
        return evaluator.compare(annotTrans, evidTrans)

    def testExact(self):
        # test allowExtension with fake psl that exactly matches
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["651", "0", "0", "0", "0", "0", "3", "14865", "+", "ENST00000477874.1_aln", "651", "0", "651", "chr22", "50818468", "17084953", "17100469", "4", "276,147,113,115,", "0,276,423,536,", "17084953,17097796,17098774,17100354,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport.support, EvidenceSupport.good)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport.support, EvidenceSupport.good)

    def testExtend5(self):
        # test allowExtension with fake psl on extended on 5 side
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["751", "0", "0", "0", "0", "0", "4", "15218", "+", "ENST00000477874.1_aln", "751", "0", "751", "chr22", "50818468", "17084500", "17100469", "5", "100,276,147,113,115,", "0,100,376,523,636,", "17084500,17084953,17097796,17098774,17100354,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport.support, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport.support, EvidenceSupport.extends_exons)

    def testExtend5Inexact(self):
        # test allowExtension with fake psl on extended on 5 side and 5'
        # annotated exon being shorter than evidence, which is not allowed
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["704", "0", "0", "0", "0", "0", "4", "15265", "+", "ENST00000477874.1_aln", "704", "0", "704", "chr22", "50818468", "17084500", "17100469", "5", "100,229,147,113,115,", "0,100,329,476,589,", "17084500,17085000,17097796,17098774,17100354,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport.support, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport.support, EvidenceSupport.feat_mismatch)

    def testExtend3(self):
        # test allowExtension with fake psl on extended on 3' side
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["741", "0", "0", "0", "0", "0", "4", "14896", "+", "ENST00000477874.1_aln", "741", "0", "741", "chr22", "50818468", "17084953", "17100590", "5", "276,147,113,115,90,", "0,276,423,536,651,", "17084953,17097796,17098774,17100354,17100500,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport.support, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport.support, EvidenceSupport.extends_exons)

    def testExtend3Inexact(self):
        # test allowExtension with fake psl on extended on 3' side
        # and 3' annotated exon being shorter than evidence, which is
        # not allowed
        annotTrans = self._getAnnot("ENST00000477874.1")
        evidRawPsl = ["682", "0", "0", "0", "0", "0", "4", "15165", "+", "ENST00000477874.1_aln", "682", "0", "682", "chr22", "50818468", "17084953", "17100800", "5", "276,147,113,46,100,", "0,276,423,536,582,", "17084953,17097796,17098774,17100354,17100700,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport.support, EvidenceSupport.feat_count_mismatch)
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport.support, EvidenceSupport.feat_mismatch)

    def testExtend3Regression1(self):
        # The evidence has an additional 5' exon.  The second to last annotation exon
        # mostly a single exon of BC071657.1, ending differently by one base.
        # This is followed by an annotations exon that is in an intron of the evidence.
        # This was incorrectly called as an extension.  This is with the Ensembl
        # alignment, UCSC aligned differently.
        #    annotation          evidence
        #    ...                 11012670-11012743   5' extension
        #    11013715-11013965   11013715-11013965
        #    11016843-11017007   11016843-11017007
        #    11018732-11018873   11018732-11018873
        #    11020428-11020599   11020428-11020599
        #    11022123-11022241   11022123-11024042   different end
        #    11023192-11024183   ...                 in intron
        #    ...                 11034700-11034727
        #    ...                 11122510-11122529
        annotTrans = self._getAnnot("ENST00000315091.7")
        evidRawPsl = ["2763", "0", "0", "0", "0", "0", "8", "107096", "+", "BC071657.1", "2833", "0", "2763", "chr1", "248956422", "11012670", "11122529", "9", "73,250,164,141,171,1919,17,9,19,", "0,73,323,487,628,799,2718,2735,2744,", "11012670,11013715,11016843,11018732,11020428,11022123,11034700,11034718,11122510,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport.support, EvidenceSupport.feat_mismatch)

    def testEst3Minus(self):
        annotTrans = self._getAnnot("ENST00000191063.8")
        evidRawPsl = ["857", "11", "0", "8", "3", "17", "8", "10236", "--", "AL558145.3", "1038", "51", "944", "chr10", "133797422", "5878172", "5889284", "11", "5,12,31,114,80,221,43,109,162,87,12,", "94,101,113,145,259,339,560,603,712,874,975,", "127908138,127908143,127908156,127908187,127908302,127909355,127911657,127913345,127914255,127919135,127919238,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, False)
        self.assertEqual(evidSupport.support, EvidenceSupport.large_indel_size)

    def testPolymophic(self):
        # test polymorphic with fake psl
        annotTrans = self._getAnnot("ENST00000400588.5")
        evidRawPsl = ["2627", "0", "0", "0", "1", "3", "10", "43660", "-", "ENST00000400588.5_poly", "2630", "0", "2630", "chr22", "50818468", "16961935", "17008222", "12", "941,105,97,91,146,119,86,251,208,304,161,118,", "0,941,1046,1143,1234,1383,1502,1588,1839,2047,2351,2512,", "16961935,16963724,16964765,16965177,16966099,16966245,16968297,16969942,16987959,16991872,17007940,17008104,"]
        evidSupport = self._rawPslCmpr(annotTrans, evidRawPsl, True)
        self.assertEqual(evidSupport.support, EvidenceSupport.polymorphic)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidCompareTest))
    return ts


if __name__ == '__main__':
    unittest.main()
