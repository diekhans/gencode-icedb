"""
Convert PSL to features, closing (and tracking gaps)
"""
from __future__ import print_function
from pycbio.hgdata.rangeFinder import RangeFinder
from pycbio.hgdata.hgLite import PslDbTable
from gencode_icedb.genome import spliceSitesClassifyStrand
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.tsl.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures


class EvidencePslFactory(object):
    """
    Factory to create evidence features from PSLs
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
        return ExonFeature(self, psl.blocks[iBlkStart].tStart, psl.blocks[iBlkEnd - 1].tEnd,
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
        return IntronFeature(self, psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart,
                             qDeleteBases, startBases, endBases, spliceSites)

    def fromPsl(self, psl):
        "convert a psl to an TranscriptFeatures object"
        if psl.getTStrand() != '+':
            raise Exception("unexpected target `-' strand on {} ".format(psl.qName))
        return TranscriptFeatures(psl.tName, psl.getQStrand(), psl.tStart, psl.tEnd, psl.tSize,
                                  psl.qName, psl.qStart, psl.qEnd, psl.qSize,
                                  self.__buildFeatures(psl))


class EvidenceFeatureMap(list):
    "map by exon coordinates of TranscriptFeatures objects"

    def __init__(self):
        self.transcripts = []
        self.exonRangeMap = RangeFinder()

    def overlapping(self, chrom, start, end, strand=None):
        "generator over ExonFeatures overlaping the specified range"
        return self.exonRangeMap.overlapping(chrom, start, end, strand)

    def addTranscript(self, transFeatures):
        self.transcripts.append(transFeatures)
        for evidFeat in transFeatures:
            if isinstance(evidFeat, ExonFeature):
                self.exonRangeMap.add(transFeatures.chrom, evidFeat.start, evidFeat.end, transFeatures, transFeatures.strand)

    @staticmethod
    def dbFactory(conn, table, chrom, start, end, genomeReader):
        "constructor from a sqlite3 databases"
        pslDbTable = PslDbTable(conn, table)
        evidFactory = EvidencePslFactory(genomeReader)
        evidFeatureMap = EvidenceFeatureMap()
        for psl in pslDbTable.getTRangeOverlap(chrom, start, end):
            evidFeatureMap.addTranscript(evidFactory.fromPsl(psl))
        return evidFeatureMap
