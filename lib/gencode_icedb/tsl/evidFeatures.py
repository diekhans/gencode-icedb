"""
Convert PSL to features, closing (and tracking gaps)
"""
from __future__ import print_function
from pycbio.hgdata import dnaOps
from pycbio.hgdata.rangeFinder import RangeFinder
from pycbio.sys.symEnum import SymEnum
from pycbio.hgdata.hgLite import PslDbTable


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
            start, end = dnaOps.reverseCoords(start, end, self.size)
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
    __slots__ = ("parent", "start", "end")

    def __init__(self, parent, start, end):
        self.parent, self.start, self.end = parent, start, end


class EvidExon(EvidFeature):
    __slots__ = ("qInsertBases", "tInsertBases")

    def __init__(self, parent, start, end, qInsertBases, tInsertBases):
        super(EvidExon, self).__init__(parent, start, end)
        self.qInsertBases, self.tInsertBases = qInsertBases, tInsertBases

    def __str__(self):
        return "exon {}-{} qIns={} tIns={}".format(self.start, self.end,
                                                   self.qInsertBases, self.tInsertBases)


class EvidIntron(EvidFeature):
    __slots__ = ("qDeleteBases", "startBases", "endBases", "spliceSites")

    def __init__(self, parent, start, end, qDeleteBases, startBases, endBases, spliceSites):
        super(EvidIntron, self).__init__(parent, start, end)
        self.qDeleteBases, self.startBases, self.endBases, self.spliceSites = qDeleteBases, startBases, endBases, spliceSites

    def __str__(self):
        return "intron {}-{} qDel={} sjBases={}...{} ({})".format(self.start, self.end, self.qDeleteBases,
                                                                  self.startBases, self.endBases, self.spliceSites)


class EvidFeatures(list):
    """
    Set of features derived from a PSL alignment, features are kept in genomic
    order.
    """
    minIntronSize = 30

    def __init__(self, psl, genomeReader):
        if psl.getTStrand() != '+':
            raise Exception("unexpected target `-' strand on {} ".format(psl.qName))
        self.chrom = psl.tName
        self.strand = psl.getQStrand()
        self.start, self.end = psl.tStart, psl.tEnd
        self.size = psl.tSize
        self.qName = psl.qName
        self.qStart, self.qEnd = psl.qStart, psl.qEnd
        self.qSize = psl.qSize
        self.__buildFeatures(psl, genomeReader)

    def __str__(self):
        return "t={}:{}-{} {}, q={}:{}={} {}".format(self.chrom, self.start, self.end, self.strand,
                                                     self.qName, self.qStart, self.qEnd, self.qSize)

    def __buildFeatures(self, psl, genomeReader):
        iBlkStart = 0
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self.__findExonEnd(psl, iBlkStart)
            self.__makeExon(psl, iBlkStart, iBlkEnd)
            if iBlkEnd < len(psl.blocks):
                self.__makeIntron(psl, iBlkEnd, genomeReader)
            iBlkStart = iBlkEnd

    def __tGapSize(self, psl, iBlk):
        "size of gap before the block"
        return psl.blocks[iBlk].tStart - psl.blocks[iBlk - 1].tEnd

    def __findExonEnd(self, psl, iBlkStart):
        "finds half-open end of blocks covering current exon"
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(psl.blocks)) and (self.__tGapSize(psl, iBlkEnd) < self.minIntronSize):
            iBlkEnd += 1
        return iBlkEnd

    def __makeExon(self, psl, iBlkStart, iBlkEnd):
        qInsertBases = tInsertBases = 0
        for iBlk in xrange(iBlkStart + 1, iBlkEnd):
            qInsertBases += psl.blocks[iBlk].qStart - psl.blocks[iBlk - 1].qEnd
            tInsertBases += psl.blocks[iBlk].tStart - psl.blocks[iBlk - 1].tEnd
        self.append(EvidExon(self, psl.blocks[iBlkStart].tStart, psl.blocks[iBlkEnd - 1].tEnd,
                             qInsertBases, tInsertBases))

    def __makeIntron(self, psl, iBlkNext, genomeReader):
        qDeleteBases = psl.blocks[iBlkNext].qStart - psl.blocks[iBlkNext - 1].qEnd
        startBases = genomeReader.get(psl.tName, psl.blocks[iBlkNext - 1].tEnd,
                                      psl.blocks[iBlkNext - 1].tEnd + 2)
        endBases = genomeReader.get(psl.tName, psl.blocks[iBlkNext].tStart - 2,
                                    psl.blocks[iBlkNext].tStart)
        if self.strand == '+':
            spliceSites = spliceSitesClassify(startBases, endBases)
        else:
            spliceSites = spliceSitesClassify(dnaOps.reverseComplement(endBases), dnaOps.reverseComplement(startBases))
        self.append(EvidIntron(self, psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart,
                               qDeleteBases, startBases, endBases, spliceSites))


class EvidFeaturesMap(list):
    "map by exon coordinates of EvidFeatures objects"

    def __init__(self, psls, genomeReader):
        self.featuresList = []
        self.exonRangeMap = RangeFinder()
        for psl in psls:
            self.__addPsl(psl, genomeReader)

    @staticmethod
    def dbFactory(conn, table, chrom, start, end, genomeReader):
        "constructor from a sqlite3 databases"
        pslDbTable = PslDbTable(conn, table)
        psls = pslDbTable.getTRangeOverlap(chrom, start, end)
        return EvidFeaturesMap(psls, genomeReader)

    def __addPsl(self, psl, genomeReader):
        evidFeatures = EvidFeatures(psl, genomeReader)
        self.featuresList.append(evidFeatures)
        for evidFeat in evidFeatures:
            if isinstance(evidFeat, EvidExon):
                self.exonRangeMap.add(evidFeatures.chrom, evidFeat.start, evidFeat.end, evidFeatures, evidFeatures.strand)

    def overlapping(self, chrom, start, end, strand=None):
        "generator over ExonFeatures overlaping the specified range"
        return self.exonRangeMap.overlapping(chrom, start, end, strand)
