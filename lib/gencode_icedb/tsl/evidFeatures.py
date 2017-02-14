"""
Convert PSL to features, closing (and tracking gaps)
"""
from __future__ import print_function


class SeqRegion(object):
    def __init__(self, twoBitReader, chrom, start, end):
        self.chrom, self.start, self.end = chrom, start, end
        self.twoBitSeq = twoBitReader[chrom]


class EvidFeature(object):
    __slots__ = ("start", "end")

    def __init__(self, start, end):
        self.start, self.end = start, end


class EvidExon(EvidFeature):
    __slots__ = ("qInsertBases", "tInsertBases")

    def __init__(self, start, end, qInsertBases, tInsertBases):
        super(EvidExon, self).__init__(start, end)
        self.qInsertBases, self.tInsertBases = qInsertBases, tInsertBases

    def __str__(self):
        return "exon {}-{} qIns={} tIns={}".format(self.start, self.end, self.qInsertBases, self.tInsertBases)


class EvidIntron(EvidFeature):
    __slots__ = ("qDeleteBases", )

    def __init__(self, start, end, qDeleteBases):
        super(EvidIntron, self).__init__(start, end)
        self.qDeleteBases = qDeleteBases

    def __str__(self):
        return "intron {}-{} qDel={}".format(self.start, self.end, self.qDeleteBases)


class EvidFeatures(list):
    """
    Set of features derived from a PSL alignment, features are in genomic order.
    """
    minIntronSize = 30

    def __init__(self, psl, tStrand, seqRegion):
        if psl.getTStrand() != tStrand:
            psl = psl.reverseComplement()
        self.chrom = psl.tStart
        self.strand = tStrand
        self.qName = psl.qName
        self.qStart = psl.qStart
        self.qEnd = psl.qEnd
        self.qSize = psl.qSize
        self.__buildFeatures(psl)

    def __str__(self):
        return "t={} {}, q={}:{}={} {}".format(self.chrom, self.strand,
                                               self.qName, self.qStart,
                                               self.qEnd, self.qSize)

    def __buildFeatures(self, psl):
        iBlkStart = 0
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self.__findExonEnd(psl, iBlkStart)
            self.__makeExon(psl, iBlkStart, iBlkEnd)
            if iBlkEnd < len(psl.blocks):
                self.__makeIntron(psl, iBlkEnd)
            iBlkStart = iBlkEnd

    def __tGapSize(self, psl, iBlkStart, iBlkEnd):
        return psl.blocks[iBlkEnd].tStart - psl.blocks[iBlkStart].tEnd

    def __findExonEnd(self, psl, iBlkStart):
        "finds half-open end of blocks covering current exon"
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(psl.blocks)) and (self.__tGapSize(psl, iBlkStart, iBlkEnd) < self.minIntronSize):
            iBlkEnd += 1
        return iBlkEnd

    def __makeExon(self, psl, iBlkStart, iBlkEnd):
        qInsertBases = tInsertBases = 0
        for iBlk in xrange(iBlkStart + 1, iBlkEnd):
            qInsertBases += psl.blocks[iBlk].qStart - psl.blocks[iBlk - 1].qEnd
            tInsertBases += psl.blocks[iBlk].tStart - psl.blocks[iBlk - 1].tEnd
        self.append(EvidExon(psl.blocks[iBlkStart].tStart,
                             psl.blocks[iBlkEnd - 1].tEnd,
                             qInsertBases, tInsertBases))

    def __makeIntron(self, psl, iBlkNext):
        qDeleteBases = psl.blocks[iBlkNext].qStart - psl.blocks[iBlkNext - 1].qEnd
        self.append(EvidIntron(psl.blocks[iBlkNext - 1].tEnd,
                               psl.blocks[iBlkNext].tStart,
                               qDeleteBases))
