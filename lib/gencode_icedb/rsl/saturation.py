"""Support for estimating saturation for intron support"""
from __future__ import print_function
from gencode_icedb.rsl.gencodeIntronEvid import IntronSupportLevel, intronEvidSupportLevel
from pycbio.sys import fileOps


class SaturationCounts(object):
    """counts of introns a different levels of support"""

    tsvAnnonCols = {l: "{}_annot_cnt".format(l) for l in IntronSupportLevel}
    tsvNovelCols = {l: "{}_novel_cnt".format(l) for l in IntronSupportLevel}
    tsvHeader = tuple(["numExprs", "annotCnt", "novelCnt"]
                      + [tsvAnnonCols[l] for l in IntronSupportLevel]
                      + [tsvNovelCols[l] for l in IntronSupportLevel])

    tsvTypeMap = {h: int for h in tsvHeader}

    def __init__(self, numExprs):
        self.numExprs = numExprs
        # indexed by IntronSupportLevel
        self.annotCnts = {l: 0 for l in IntronSupportLevel}
        self.novelCnts = {l: 0 for l in IntronSupportLevel}

    def count(self, isNovel, sjcnts):
        level = intronEvidSupportLevel(sjcnts.numUniqueMapReads)
        if isNovel:
            self.novelCnts[level] += 1
        else:
            self.annotCnts[level] += 1

    def sumTsvRow(self, row):
        for l in IntronSupportLevel:
            self.annotCnts[l] += getattr(row, self.tsvAnnonCols[l])
            self.novelCnts[l] += getattr(row, self.tsvNovelCols[l])

    @staticmethod
    def writeHeader(fh):
        fileOps.prRow(fh, SaturationCounts.tsvHeader)

    def write(self, fh):
        fileOps.prRow(fh, [self.numExprs, sum(self.annotCnts.values()), sum(self.novelCnts.values())]
                      + [self.annotCnts[l] for l in IntronSupportLevel]
                      + [self.novelCnts[l] for l in IntronSupportLevel])
