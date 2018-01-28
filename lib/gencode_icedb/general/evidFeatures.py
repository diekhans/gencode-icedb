"""
Convert alignments (PSL) to features, closing and tracking gaps.
"""
from __future__ import print_function
from collections import defaultdict
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
        features = []
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self.__findExonEnd(psl, iBlkStart)
            features.append(self.__makeExon(psl, iBlkStart, iBlkEnd, trans, len(features)))
            if iBlkEnd < len(psl.blocks):
                features.append(self.__makeIntron(psl, iBlkEnd, trans, len(features)))
            iBlkStart = iBlkEnd
        return features

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
        alignFeatures.append(AlignedFeature(exon, len(alignFeatures),
                                            exon.chrom.subrange(blk.tStart, blk.tEnd),
                                            exon.rna.subrange(blk.qStart, blk.qEnd)))

    def __addUnalignedFeatures(self, psl, iBlk, feat, alignFeatures):
        prevBlk = psl.blocks[iBlk - 1]
        blk = psl.blocks[iBlk]
        if blk.qStart > prevBlk.qEnd:
            alignFeatures.append(RnaInsertFeature(feat, len(alignFeatures),
                                                  feat.rna.subrange(prevBlk.qEnd, blk.qStart)))
        if blk.tStart > prevBlk.tEnd:
            alignFeatures.append(ChromInsertFeature(feat, len(alignFeatures),
                                                    feat.chrom.subrange(prevBlk.tEnd, blk.tStart)))

    def __addAlignFeatures(self, psl, iBlkStart, iBlkEnd, exon):
        alignFeatures = []
        for iBlk in range(iBlkStart, iBlkEnd):
            if iBlk > iBlkStart:
                self.__addUnalignedFeatures(psl, iBlk, exon, alignFeatures)
            self.__addAlignedFeature(psl, iBlk, exon, alignFeatures)  # after since unaligned is before block
        exon.alignFeatures = tuple(alignFeatures)

    def __makeExon(self, psl, iBlkStart, iBlkEnd, trans, iFeat):
        exon = ExonFeature(trans, iFeat,
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

    def __makeIntron(self, psl, iBlkNext, trans, iFeat):
        donorSeq, acceptorSeq = self.__getSpliceSites(psl, iBlkNext)
        alignFeatures = []
        intron = IntronFeature(trans, iFeat,
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


class EvidenceMap(list):
    "map by overall coordinates of TranscriptFeatures objects"

    def __init__(self):
        self.byName = defaultdict(list)
        self.byRange = RangeFinder()

    def overlapping(self, chrom, start, end, strand=None):
        "generator over Features overlaping the specified range"
        return self.byRange.overlapping(chrom, start, end, strand)

    def add(self, trans):
        self.append(trans)
        self.byRange.add(trans.chrom.name, trans.chrom.start, trans.chrom.end, trans, trans.chrom.strand)

    @staticmethod
    def dbFactory(conn, table, genomeReader, chrom=None, start=None, end=None):
        """Factory from a sqlite3 databases.  This can load a range or an entire
        table."""
        assert (chrom is None) == (start is None) == (end is None), "chrom, start, end must all be None or all specified"
        pslDbTable = PslDbTable(conn, table)
        evidFactory = EvidencePslFactory(genomeReader)
        evidMap = EvidenceMap()
        if chrom is None:
            pslGen = pslDbTable.getAll()
        else:
            pslGen = pslDbTable.getTRangeOverlap(chrom, start, end)
        for psl in pslGen:
            evidMap.add(evidFactory.fromPsl(psl))
        return evidMap
