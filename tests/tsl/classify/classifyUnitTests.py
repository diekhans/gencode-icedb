from __future__ import print_function
import sys
import os
if __name__ == '__main__':
    rootDir = "../../.."
    sys.path = [os.path.join(rootDir, "lib"),
                os.path.join(rootDir, "extern/pycbio/lib")] + sys.path
import unittest
from pycbio.sys.testCaseBase import TestCaseBase
from pycbio.hgdata.hgLite import PslDbTable, GenePredDbTable
from gencode_icedb.tsl.supportClassify import compareMegWithEvidence
from gencode_icedb.tsl.supportDefs import EvidenceComparison


class TestEvidence(object):
    """pre-loads all test evidence from a database"""
    def __init__(self, evidDb):


class EvidCompareTest(TestCaseBase):
    """low-level comparison tests"""
