"""
Convert PSL to features, closing (and tracking gaps)
"""
from __future__ import print_function
from pycbio.hgdata import dnaOps
from pycbio.sys.symEnum import SymEnum


class SeqReader(object):
    """
    Reads sequences from twoBitReader, which can be substituted with a mock
    version.
    """
    def __init__(self, twoBitReader):
        self.twoBitReader = twoBitReader
        self.chrom = self.size = None

    def __obtainChrom(self, chrom):
        if self.chrom != chrom:
            self.twoBitSeq = self.twoBitReader[chrom]
            self.chrom = chrom
            self.size = len(self.twoBitSeq)

    def get(self, chrom, start, end, strand=None):
        self.__obtainChrom(chrom)
        if strand == '-':
            start, end = self.size - end, self.size - start
        seq = self.twoBitSeq[start:end]
        if strand == '-':
            seq = dnaOps.reverseComplement(seq)
        return seq

SpliceSites = SymEnum("SpliceSite",
                     ("spliceGT_AG", "spliceGC_AG", "spliceAT_AC", "spliceOther"))

spliceSitesMap = {
    ("gt", "ag"): SpliceSites.spliceGT_AG,
    ("gc", "ag"): SpliceSites.spliceGC_AG,
    ("at", "ac"): SpliceSites.spliceAT_AC,
}

def spliceSitesClassify(donor, acceptor):
    return spliceSitesMap.get((donor.lower(), acceptor.lower()), SpliceSites.spliceOther)

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
        return "exon {}-{} qIns={} tIns={} ".format(self.start, self.end,
                                                    self.qInsertBases, self.tInsertBases)


class EvidIntron(EvidFeature):
    __slots__ = ("qDeleteBases", "startBases", "endBases", "spliceSites")

    def __init__(self, start, end, qDeleteBases, startBases, endBases, strand):
        super(EvidIntron, self).__init__(start, end)
        self.qDeleteBases, self.startBases, self.endBases = qDeleteBases, startBases, endBases
        if strand == '+':
            self.spliceSites = spliceSitesClassify(startBases, endBases)
        else:
            self.spliceSites = spliceSitesClassify(dnaOps.reverseComplement(endBases), dnaOps.reverseComplement(startBases))
        
    def __str__(self):
        return "intron {}-{} qDel={} sjBases={}...{} ({})".format(self.start,
                                                                  self.end, self.qDeleteBases, self.startBases, self.endBases,
                                                                  self.spliceSites)


class EvidFeatures(list):
    """
    Set of features derived from a PSL alignment, features are in genomic order.
    """
    minIntronSize = 30

    def __init__(self, psl, seqRegion):
        if psl.getTStrand() != '+':
            raise Exception("expected PSL target strand of +")
        self.chrom = psl.tName
        self.strand = psl.getQStrand()
        self.qName = psl.qName
        self.qStart = psl.qStart
        self.qEnd = psl.qEnd
        self.qSize = psl.qSize
        self.__buildFeatures(psl, seqRegion)

    def __str__(self):
        return "t={} {}, q={}:{}={} {}".format(self.chrom, self.strand,
                                               self.qName, self.qStart,
                                               self.qEnd, self.qSize)

    def __buildFeatures(self, psl, seqRegion):
        iBlkStart = 0
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self.__findExonEnd(psl, iBlkStart)
            self.__makeExon(psl, iBlkStart, iBlkEnd)
            if iBlkEnd < len(psl.blocks):
                self.__makeIntron(psl, iBlkEnd, seqRegion)
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

    def __makeIntron(self, psl, iBlkNext, seqRegion):
        qDeleteBases = psl.blocks[iBlkNext].qStart - psl.blocks[iBlkNext - 1].qEnd
        startBases = seqRegion.get(psl.tName, psl.blocks[iBlkNext].tStart,
                                   psl.blocks[iBlkNext].tStart + 2, psl.getTStrand())
        endBases = seqRegion.get(psl.tName, psl.blocks[iBlkNext].tEnd - 2,
                                 psl.blocks[iBlkNext].tEnd, psl.getTStrand())
        self.append(EvidIntron(psl.blocks[iBlkNext - 1].tEnd,
                               psl.blocks[iBlkNext].tStart,
                               qDeleteBases, startBases, endBases, self.strand))
