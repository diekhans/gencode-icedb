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
from gencode_icedb.general.evidFeatures import EvidenceMap
from gencode_icedb.general.annotFeatures import AnnotationMap
from gencode_icedb.tsl.supportClassify import compareMegWithEvidence
from gencode_icedb.tsl.supportDefs import UCSC_RNA_ALN_TBL, UCSC_EST_ALN_TBL, ENSEMBL_RNA_ALN_TBL
from gencode_icedb.tsl.supportDefs import GENCODE_ANN_TBL
from gencode_icedb.tsl.supportDefs import EvidenceComparison


class TestData(object):
    """Pre-loads all test data from a database to speed up tests"""

    def __init__(self, ucscDb, evidDb, annotDb):
        self.genomeReader = GenomeReaderFactory.factoryFromUcscDb(ucscDb).obtain()
        conn = sqliteConnect(evidDb)
        try:
            self.ucscRnas = EvidenceMap.dbFactory(conn, UCSC_RNA_ALN_TBL, self.genomeReader)
            self.ucscEsts = EvidenceMap.dbFactory(conn, UCSC_EST_ALN_TBL, self.genomeReader)
            self.ensemblRna = EvidenceMap.dbFactory(conn, ENSEMBL_RNA_ALN_TBL, self.genomeReader)
        finally:
            conn.close()
        conn = sqliteConnect(annotDb)
        try:
            self.annon = AnnotationMap.dbFactory(conn, GENCODE_ANN_TBL, self.genomeReader)
        finally:
            conn.close()


class EvidCompareTest(TestCaseBase):
    """low-level comparison tests"""
    data = None

    @classmethod
    def setUpClass(cls):
        cls.evidenceCache = TestData("hg38", "output/hsEvid.db", "output/hsGencode.db")

    def testShit(self):
        pass


def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidCompareTest))
    return ts


if __name__ == '__main__':
    unittest.main()
