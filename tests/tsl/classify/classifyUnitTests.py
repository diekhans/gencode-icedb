from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from gencode_icedb.general.genome import GenomeReaderFactory
from gencode_icedb.general.gencodeDb import UcscGencodeReader
from gencode_icedb.tsl.evidenceDb import EvidenceSource, EvidenceReader
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

    def __evalWithEvid(self, annotTrans, evidSrc, evidTrans):
        sup = compareMegWithEvidence(annotTrans, evidTrans)
        print("annot: {} {}: {} => {}".format(annotTrans.rna.name, str(evidSrc), evidTrans.rna.name, str(sup)))

    def __evalAnnotTransEvidSrc(self, annotTrans, evidSrc):
        for evidTrans in self.evidenceReader.getOverlapping(evidSrc, annotTrans.chrom.name, annotTrans.chrom.start, annotTrans.chrom.end, annotTrans.rna.strand):
            self.__evalWithEvid(annotTrans, evidSrc, evidTrans)

    def __evalAnnotTrans(self, annotTrans):
        for evidSrc in EvidenceSource:
            self.__evalAnnotTransEvidSrc(annotTrans, evidSrc)

    def __classifyTest(self, annotTranses, noDiff=False):
        outTslTsv = self.getOutputFile(".tsl.tsv")
        outDetailsTsv = self.getOutputFile(".details.tsv")
        with open(outTslTsv, 'w') as tslTsvFh, open(outDetailsTsv, 'w') as detailsTsvFh:
            writeTsvHeaders(tslTsvFh, detailsTsvFh)
            classifyGeneTranscripts(self.evidenceReader, annotTranses, tslTsvFh, detailsTsvFh)
        if not noDiff:
            self.diffFiles(self.getExpectedFile(".tsl.tsv"), outTslTsv)
            self.diffFiles(self.getExpectedFile(".details.tsv"), outDetailsTsv)

    def __classifyGeneTest(self, geneId):
        geneTranses = self.gencodeReader.getByGeneId(geneId)
        if len(geneTranses) == 0:
            raise Exception("no transcripts found for {}".format(geneId))
        self.__classifyTest(geneTranses)

    def __classifyTransTest(self, transId, noDiff=False):
        transes = self.gencodeReader.getByTranscriptId(transId)
        if len(transes) == 0:
            raise Exception("no transcripts found for {}".format(transId))
        self.__classifyTest(transes, noDiff)

    def testGAB4(self):
        self.__classifyGeneTest("ENSG00000215568.7")

    def testBCR(self):
        self.__classifyGeneTest("ENSG00000186716.20")

    def testIL17RA(self):
        self.__classifyGeneTest("ENSG00000177663.13")

    def testSHOX(self):
        # PAR gene
        self.__classifyGeneTest("ENSG00000185960.13")

    def skip_testDebug(self):
        "Debug a transcript"
        self.__classifyTransTest("ENST00000359540.7", noDiff=True)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidCompareTest))
    return ts


if __name__ == '__main__':
    unittest.main()
