from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from pycbio.hgdata.hgLite import sqliteConnect, PslDbTable, GenePredDbTable
from pycbio.hgdata.rangeFinder import RangeFinder
from gencode_icedb.general.genome import GenomeReaderFactory
from gencode_icedb.general.evidFeatures import EvidencePslDbFactory
from gencode_icedb.general.annotFeatures import AnnotationGenePredDbFactory
from gencode_icedb.tsl.supportClassify import compareMegWithEvidence
from gencode_icedb.tsl.supportDefs import UCSC_RNA_ALN_TBL, UCSC_EST_ALN_TBL, ENSEMBL_RNA_ALN_TBL
from gencode_icedb.tsl.supportDefs import GENCODE_ANN_TBL
from gencode_icedb.tsl.supportDefs import EvidenceComparison


class TestData(object):
    """Pre-loads all test data from a database to speed up tests"""

    def __init__(self, ucscDb, evidDb, annotDb):
        self.genomeReader = GenomeReaderFactory.factoryFromUcscDb(ucscDb).obtain()
        self.evidDbConn = sqliteConnect(evidDb)
        self.ucscRnas = EvidencePslDbFactory(self.evidDbConn, UCSC_RNA_ALN_TBL, self.genomeReader)
        self.ucscEsts = EvidencePslDbFactory(self.evidDbConn, UCSC_EST_ALN_TBL, self.genomeReader)
        self.ensemblRnas = EvidencePslDbFactory(self.evidDbConn, ENSEMBL_RNA_ALN_TBL, self.genomeReader)
        self.annotDbConn = sqliteConnect(annotDb)
        self.annots = AnnotationGenePredDbFactory(self.annotDbConn, GENCODE_ANN_TBL, self.genomeReader)


class EvidCompareTest(TestCaseBase):
    """low-level comparison tests"""
    data = None

    @classmethod
    def setUpClass(cls):
        cls.data = TestData("hg38", "output/hsEvid.db", "output/hsGencode.db")

    def __evalWithEvid(self, annotTrans, etype, evidTrans):
        sup = compareMegWithEvidence(annotTrans, evidTrans)
        print("annot: {} {}: {} => {}".format(annotTrans.rna.name, etype, evidTrans.rna.name, str(sup)))

    def __evalAnnotTransEvidSrc(self, annotTrans, etype, evidReader):
        for evidTrans in evidReader.overlappingGen(annotTrans.chrom.name, annotTrans.chrom.start, annotTrans.chrom.end, annotTrans.rna.strand):
            self.__evalWithEvid(annotTrans, etype, evidTrans)

    def __evalAnnotTrans(self, annotTrans):
        self.__evalAnnotTransEvidSrc(annotTrans, "ucscRna", self.data.ucscRnas)
        self.__evalAnnotTransEvidSrc(annotTrans, "ucscEst", self.data.ucscEsts)
        self.__evalAnnotTransEvidSrc(annotTrans, "ensemblRna", self.data.ensemblRnas)

    def testShit(self):
        for annotTrans in self.data.annots.allGen():
            if len(annotTrans.features) >= 3:  # multi-exon only
                self.__evalAnnotTrans(annotTrans)


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidCompareTest))
    return ts


if __name__ == '__main__':
    unittest.main()
