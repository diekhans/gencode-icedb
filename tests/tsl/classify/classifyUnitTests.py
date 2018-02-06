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
from gencode_icedb.general.evidenceDb import EvidenceSource, EvidenceReader
from gencode_icedb.tsl.supportClassify import compareMegWithEvidence, SupportClassifier
from gencode_icedb.tsl.supportDefs import EvidenceEval


class EvidCompareTest(TestCaseBase):
    """low-level comparison tests"""
    UCSC_DB = "hg38"
    EVIDENCE_DB = "output/hsEvid.db"
    GENCODE_DB = "output/hsGencode.db"

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
        for evidTrans in self.evidenceReader.overlappingGen(evidSrc, annotTrans.chrom.name, annotTrans.chrom.start, annotTrans.chrom.end, annotTrans.rna.strand):
            self.__evalWithEvid(annotTrans, evidSrc, evidTrans)

    def __evalAnnotTrans(self, annotTrans):
        for evidSrc in EvidenceSource:
            self.__evalAnnotTransEvidSrc(annotTrans, evidSrc)

    def XtestShit(self):
        for annotTrans in self.gencodeReader.allGen():
            if len(annotTrans.features) >= 3:  # multi-exon only
                self.__evalAnnotTrans(annotTrans)

    def testGAB4(self):
        geneId = "ENSG00000215568.7"
        geneAnnotTranses = list(self.gencodeReader.getByGeneId(geneId))
        classifier = SupportClassifier(self.evidenceReader, self.genomeReader)
        outTslTsv = self.getOutputFile(".tsl.tsv")
        outDetailsTsv = self.getOutputFile(".details.tsv")
        with open(outTslTsv, 'w') as tslTsvFh, \
              open(outDetailsTsv, 'w') as detailsTsvFh:
            classifier.writeTsvHeaders(tslTsvFh, detailsTsvFh)
            classifier.classifyGeneTranscripts(geneAnnotTranses, tslTsvFh, detailsTsvFh)
        self.diffFiles(self.getExpectedFile(".tsl.tsv"), outTslTsv)
        self.diffFiles(self.getExpectedFile(".details.tsv"), outDetailsTsv)

def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidCompareTest))
    return ts


if __name__ == '__main__':
    unittest.main()
