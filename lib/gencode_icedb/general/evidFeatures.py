"""
Convert alignments (PSL) to features, closing and tracking gaps.
"""
from __future__ import print_function
from pycbio.hgdata.rangeFinder import RangeFinder
from pycbio.hgdata.hgLite import PslDbTable
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, AlignedFeature, ChromInsertFeature, RnaInsertFeature


class EvidencePslFactory(object):
    """
    Factory to create evidence features from PSLs
    """

    def __init__(self, genomeReader):
        """genomeReader maybe None if splice sites are not desired """
        self.genomeReader = genomeReader

    def __buildFeatures(self, psl, trans):
        iBlkStart = 0
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self.__findExonEnd(psl, iBlkStart)
            yield self.__makeExon(psl, iBlkStart, iBlkEnd, trans)
            if iBlkEnd < len(psl.blocks):
                yield self.__makeIntron(psl, iBlkEnd, trans)
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

    def __addAlignedFeature(self, psl, iBlk, exon, alignFeatures):
        blk = psl.blocks[iBlk]
        alignFeatures.append(AlignedFeature(exon, blk.tStart, blk.tEnd, blk.qStart, blk.qEnd))

    def __addUnalignedFeatures(self, psl, iBlk, exon, alignFeatures):
        prevBlk = psl.blocks[iBlk - 1]
        blk = psl.blocks[iBlk]
        if blk.qStart > prevBlk.qEnd:
            alignFeatures.append(RnaInsertFeature(exon, prevBlk.qEnd, blk.qStart))
        if blk.tStart > prevBlk.tEnd:
            alignFeatures.append(ChromInsertFeature(exon, prevBlk.tEnd, blk.tStart))

    def __addAlignFeatures(self, psl, iBlkStart, iBlkEnd, exon):
        alignFeatures = []
        for iBlk in range(iBlkStart, iBlkEnd):
            if iBlk > iBlkStart:
                self.__addUnalignedFeatures(psl, iBlk, exon, alignFeatures)
            self.__addAlignedFeature(psl, iBlk, exon, alignFeatures)  # after since unaligned is before block
        exon.alignFeatures = tuple(alignFeatures)

    def __makeExon(self, psl, iBlkStart, iBlkEnd, trans):
        exon = ExonFeature(trans, psl.blocks[iBlkStart].tStart, psl.blocks[iBlkEnd - 1].tEnd,
                           psl.blocks[iBlkStart].qStart, psl.blocks[iBlkEnd - 1].qEnd)
        self.__addAlignFeatures(psl, iBlkStart, iBlkEnd, exon)
        return exon

    def __getSpliceSites(self, psl, iBlkNext):
        if self.genomeReader is None:
            return (None, None)
        else:
            return spliceJuncsGetSeqs(self.genomeReader, psl.tName, psl.blocks[iBlkNext - 1].tEnd,
                                      psl.blocks[iBlkNext].tStart, psl.getQStrand())

    def __makeIntron(self, psl, iBlkNext, trans):
        donorSeq, acceptorSeq = self.__getSpliceSites(psl, iBlkNext)
        return IntronFeature(trans, psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart,
                             psl.blocks[iBlkNext - 1].qEnd, psl.blocks[iBlkNext].qStart,
                             donorSeq, acceptorSeq)

    def fromPsl(self, psl):
        "convert a psl to an TranscriptFeatures object"
        if psl.getTStrand() != '+':
            raise Exception("unexpected target `-' strand on {} ".format(psl.qName))
        trans = TranscriptFeatures(psl.tName, psl.getTStrand(), psl.tStart, psl.tEnd, psl.tSize,
                                   psl.qName, psl.getQStrand(), psl.qStart, psl.qEnd, psl.qSize)
        trans.features = tuple(self.__buildFeatures(psl, trans))
        return trans


class EvidenceFeatureMap(list):
    "map by exon coordinates of TranscriptFeatures objects"

    def __init__(self):
        self.transcripts = []
        self.exonRangeMap = RangeFinder()

    def overlapping(self, chrom, start, end, strand=None):
        "generator over ExonFeatures overlaping the specified range"
        return self.exonRangeMap.overlapping(chrom, start, end, strand)

    def addTranscript(self, trans):
        self.transcripts.append(trans)
        for feat in trans.features:
            if isinstance(feat, ExonFeature):
                self.exonRangeMap.add(trans.chrom, feat.chromStart, feat.chromEnd, trans, trans.chromStrand)

    @staticmethod
    def dbFactory(conn, table, chrom, start, end, genomeReader):
        "constructor from a sqlite3 databases"
        pslDbTable = PslDbTable(conn, table)
        evidFactory = EvidencePslFactory(genomeReader)
        evidFeatureMap = EvidenceFeatureMap()
        for psl in pslDbTable.getTRangeOverlap(chrom, start, end):
            evidFeatureMap.addTranscript(evidFactory.fromPsl(psl))
        return evidFeatureMap
