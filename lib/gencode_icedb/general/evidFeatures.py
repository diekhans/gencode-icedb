"""
Convert alignments (PSL) to features, closing and tracking gaps.
"""
from __future__ import print_function
from pycbio.hgdata.rangeFinder import RangeFinder
from pycbio.hgdata.hgLite import PslDbTable
from pycbio.hgdata.coords import Coords
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
        alignFeatures.append(AlignedFeature(exon,
                                            exon.chrom.subrange(blk.tStart, blk.tEnd),
                                            exon.rna.subrange(blk.qStart, blk.qEnd)))

    def __addUnalignedFeatures(self, psl, iBlk, feat, alignFeatures):
        prevBlk = psl.blocks[iBlk - 1]
        blk = psl.blocks[iBlk]
        if blk.qStart > prevBlk.qEnd:
            alignFeatures.append(RnaInsertFeature(feat,
                                                  feat.rna.subrange(prevBlk.qEnd, blk.qStart)))
        if blk.tStart > prevBlk.tEnd:
            alignFeatures.append(ChromInsertFeature(feat,
                                                    feat.chrom.subrange(prevBlk.tEnd, blk.tStart)))

    def __addAlignFeatures(self, psl, iBlkStart, iBlkEnd, exon):
        alignFeatures = []
        for iBlk in range(iBlkStart, iBlkEnd):
            if iBlk > iBlkStart:
                self.__addUnalignedFeatures(psl, iBlk, exon, alignFeatures)
            self.__addAlignedFeature(psl, iBlk, exon, alignFeatures)  # after since unaligned is before block
        exon.alignFeatures = tuple(alignFeatures)

    def __makeExon(self, psl, iBlkStart, iBlkEnd, trans):
        exon = ExonFeature(trans,
                           trans.chrom.subrange(psl.blocks[iBlkStart].tStart, psl.blocks[iBlkEnd - 1].tEnd),
                           trans.rna.subrange(psl.blocks[iBlkStart].qStart, psl.blocks[iBlkEnd - 1].qEnd))
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
        alignFeatures = []
        intron = IntronFeature(trans,
                               trans.chrom.subrange(psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart),
                               trans.rna.subrange(psl.blocks[iBlkNext - 1].qEnd, psl.blocks[iBlkNext].qStart),
                               donorSeq, acceptorSeq)
        self.__addUnalignedFeatures(psl, iBlkNext, intron, alignFeatures)
        intron.alignFeatures = tuple(alignFeatures)
        return intron

    def fromPsl(self, psl):
        "convert a psl to an TranscriptFeatures object"
        chrom = Coords(psl.tName, psl.tStart, psl.tEnd, '+', psl.tSize)
        if psl.getTStrand() == '-':
            chrom = chrom.reverse()
        rna = Coords(psl.qName, psl.qStart, psl.qEnd, '+', psl.qSize)
        if psl.getQStrand() == '-':
            rna = rna.reverse()
        trans = TranscriptFeatures(chrom, rna)
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
                self.exonRangeMap.add(trans.chrom.name, feat.chrom.start, feat.chrom.end, trans, trans.chrom.strand)

    @staticmethod
    def dbFactory(conn, table, chrom, start, end, genomeReader):
        "constructor from a sqlite3 databases"
        pslDbTable = PslDbTable(conn, table)
        evidFactory = EvidencePslFactory(genomeReader)
        evidFeatureMap = EvidenceFeatureMap()
        for psl in pslDbTable.getTRangeOverlap(chrom, start, end):
            evidFeatureMap.addTranscript(evidFactory.fromPsl(psl))
        return evidFeatureMap
