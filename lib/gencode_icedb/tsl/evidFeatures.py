"""
Convert PSL to features, closing (and tracking gaps)
"""
from __future__ import print_function
from pycbio.hgdata.rangeFinder import RangeFinder
from pycbio.hgdata.hgLite import PslDbTable
from gencode_icedb.sequence import spliceSitesClassifyStrand
from gencode_icedb.tsl import minIntronSize


class EvidFeature(object):
    "a feature of evidence aligned to the genome"
    __slots__ = ("parent", "start", "end")

    def __init__(self, parent, start, end):
        self.parent, self.start, self.end = parent, start, end


class EvidExon(EvidFeature):
    "exon from alignment, with target deletions closed"
    __slots__ = ("qInsertBases", "tInsertBases")

    def __init__(self, parent, start, end, qInsertBases, tInsertBases):
        super(EvidExon, self).__init__(parent, start, end)
        self.qInsertBases, self.tInsertBases = qInsertBases, tInsertBases

    def __str__(self):
        return "exon {}-{} qIns={} tIns={}".format(self.start, self.end,
                                                   self.qInsertBases, self.tInsertBases)


class EvidIntron(EvidFeature):
    "intron from alignment"
    __slots__ = ("qDeleteBases", "startBases", "endBases", "spliceSites")

    def __init__(self, parent, start, end, qDeleteBases, startBases, endBases, spliceSites):
        super(EvidIntron, self).__init__(parent, start, end)
        self.qDeleteBases, self.startBases, self.endBases, self.spliceSites = qDeleteBases, startBases, endBases, spliceSites

    def __str__(self):
        if self.startBases is None:
            sjDesc = None
        else:
            sjDesc = "{}...{} ({})".format(self.startBases, self.endBases, self.spliceSites)
        return "intron {}-{} qDel={} sjBases={}".format(self.start, self.end, self.qDeleteBases, sjDesc)


class EvidTranscript(list):
    """
    Set of features for a transcript derived from a PSL alignment, features are
    kept in genomic order.
    """

    def __init__(self, chrom, strand, start, end, size, qName, qStart, qEnd, qSize,
                 features):
        self.chrom, self.strand, self.start, self.end, self.size = chrom, strand, start, end, size
        self.qName, self.qStart, self.qEnd, self.qSize = qName, qStart, qEnd, qSize
        self.extend(features)

    def __str__(self):
        return "t={}:{}-{} {}, q={}:{}={} {}".format(self.chrom, self.start, self.end, self.strand,
                                                     self.qName, self.qStart, self.qEnd, self.qSize)


class EvidTranscriptPslFactory(object):
    """
    factory to create evidence features from PSLs
    """

    def __init__(self, genomeReader):
        """genomeReader maybe None if splice sites are not desired """
        self.genomeReader = genomeReader

    def __buildFeatures(self, psl):
        iBlkStart = 0
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self.__findExonEnd(psl, iBlkStart)
            yield self.__makeExon(psl, iBlkStart, iBlkEnd)
            if iBlkEnd < len(psl.blocks):
                yield self.__makeIntron(psl, iBlkEnd)
            iBlkStart = iBlkEnd

    def __tGapSize(self, psl, iBlk):
        "size of gap before the block"
        return psl.blocks[iBlk].tStart - psl.blocks[iBlk - 1].tEnd

    def __findExonEnd(self, psl, iBlkStart):
        "finds half-open end of blocks covering current exon"
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(psl.blocks)) and (self.__tGapSize(psl, iBlkEnd) < minIntronSize):
            iBlkEnd += 1
        return iBlkEnd

    def __makeExon(self, psl, iBlkStart, iBlkEnd):
        qInsertBases = tInsertBases = 0
        for iBlk in xrange(iBlkStart + 1, iBlkEnd):
            qInsertBases += psl.blocks[iBlk].qStart - psl.blocks[iBlk - 1].qEnd
            tInsertBases += psl.blocks[iBlk].tStart - psl.blocks[iBlk - 1].tEnd
        return EvidExon(self, psl.blocks[iBlkStart].tStart, psl.blocks[iBlkEnd - 1].tEnd,
                        qInsertBases, tInsertBases)

    def __getSpliceSites(self, psl, iBlkNext):
        startBases = self.genomeReader.get(psl.tName, psl.blocks[iBlkNext - 1].tEnd,
                                           psl.blocks[iBlkNext - 1].tEnd + 2)
        endBases = self.genomeReader.get(psl.tName, psl.blocks[iBlkNext].tStart - 2,
                                         psl.blocks[iBlkNext].tStart)
        spliceSites = spliceSitesClassifyStrand(psl.getQStrand(), startBases, endBases)
        return startBases, endBases, spliceSites
        
    def __makeIntron(self, psl, iBlkNext):
        qDeleteBases = psl.blocks[iBlkNext].qStart - psl.blocks[iBlkNext - 1].qEnd
        if self.genomeReader is None:
            startBases = endBases = spliceSites = None
        else:
            startBases, endBases, spliceSites = self.__getSpliceSites(psl, iBlkNext)
        return EvidIntron(self, psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart,
                          qDeleteBases, startBases, endBases, spliceSites)

    def fromPsl(self, psl):
        "convert a psl to an EvidTranscript"
        if psl.getTStrand() != '+':
            raise Exception("unexpected target `-' strand on {} ".format(psl.qName))
        return EvidTranscript(psl.tName, psl.getQStrand(), psl.tStart, psl.tEnd, psl.tSize,
                              psl.qName, psl.qStart, psl.qEnd, psl.qSize,
                              self.__buildFeatures(psl))


class EvidFeatureMap(list):
    "map by exon coordinates of EvidFeatures objects"

    def __init__(self):
        self.transcripts = []
        self.exonRangeMap = RangeFinder()

    def overlapping(self, chrom, start, end, strand=None):
        "generator over ExonFeatures overlaping the specified range"
        return self.exonRangeMap.overlapping(chrom, start, end, strand)

    def addTranscript(self, evidTrans):
        self.transcripts.append(evidTrans)
        for evidFeat in evidTrans:
            if isinstance(evidFeat, EvidExon):
                self.exonRangeMap.add(evidTrans.chrom, evidFeat.start, evidFeat.end, evidTrans, evidTrans.strand)

    @staticmethod
    def dbFactory(conn, table, chrom, start, end, genomeReader):
        "constructor from a sqlite3 databases"
        pslDbTable = PslDbTable(conn, table)
        evidFactory = EvidTranscriptPslFactory(genomeReader)
        evidFeatureMap = EvidFeatureMap()
        for psl in pslDbTable.getTRangeOverlap(chrom, start, end):
            evidFeatureMap.addTranscript(evidFactory.fromPsl(psl))
        return evidFeatureMap
