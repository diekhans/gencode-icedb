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
from gencode_icedb.tsl.evidFeatures import EvidFeatures


class EvidenceTests(TestCaseBase):
    set1Psls = None

    def __requireSet1Psls(self):
        if self.set1Psls is None:
            self.set1Psls = PslTbl(self.getInputFile("set1.ucsc-mrna.psl"),
                                   qNameIdx=True)

    def __getSet1Psl(self, acc):
        self.__requireSet1Psls()
        return self.set1Psls.qNameMap[acc][0]

    def __assertFeatures(self, feats, expectFeatsStr, expectFeatStrs):
        self.assertEqual(str(feats), expectFeatsStr)
        self.assertEqual([str(f) for f in feats], expectFeatStrs)
    
    def testAF010310(self):
        psl = self.__getSet1Psl("AF010310.1")
        feats = EvidFeatures(psl, '-', None)
        self.__assertFeatures(feats, "t=18900294 -, q=AF010310.1:0=888 901",
                              ['exon 32398640-32398738 qIns=1 tIns=4',
                               'intron 32398738-32400065 qDel=0',
                               'exon 32400065-32400164 qIns=0 tIns=0',
                               'intron 32400164-32403527 qDel=0',
                               'exon 32403527-32403616 qIns=0 tIns=0',
                               'intron 32403616-32403691 qDel=0',
                               'exon 32403691-32403763 qIns=1 tIns=0',
                               'intron 32403763-32403763 qDel=1',
                               'exon 32403763-32403897 qIns=2 tIns=0',
                               'intron 32403897-32403897 qDel=1',
                               'exon 32403897-32403932 qIns=3 tIns=0',
                               'intron 32403932-32403938 qDel=10',
                               'exon 32403938-32403987 qIns=6 tIns=1',
                               'intron 32403987-32403987 qDel=1',
                               'exon 32403987-32404031 qIns=3 tIns=0',
                               'intron 32404031-32404031 qDel=1',
                               'exon 32404031-32404272 qIns=3 tIns=1'])


    def testX96484(self):
        psl = self.__getSet1Psl("X96484.1")
        feats = EvidFeatures(psl, '+', None)
        self.__assertFeatures(feats, "t=18893922 +, q=X96484.1:48=1067 1080",
                              ['exon 18893922-18893997 qIns=0 tIns=0',
                               'intron 18893997-18894077 qDel=0',
                               'exon 18894077-18894238 qIns=1 tIns=0',
                               'intron 18894238-18897684 qDel=0',
                               'exon 18897684-18897785 qIns=0 tIns=0',
                               'intron 18897785-18898400 qDel=0',
                               'exon 18898400-18898541 qIns=0 tIns=0',
                               'intron 18898541-18899052 qDel=0',
                               'exon 18899052-18899592 qIns=0 tIns=0'])
        
    def testSet1(self):
        # just run through to see if they all convert
        self.__requireSet1Psls()
        cnt = 0
        for psl in self.set1Psls:
            EvidFeatures(psl, psl.getTStrand(), None)
            cnt += 1
        self.assertEqual(cnt, 81)

def suite():
    ts = unittest.TestSuite()
    ts.addTest(unittest.makeSuite(EvidenceTests))
    return ts

if __name__ == '__main__':
    unittest.main()
