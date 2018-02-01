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
    Factory to create evidence features from PSLs.
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

class EvidencePslDbFactory(EvidencePslFactory):
    """
    Factory to create evidence features from a PSLs in an sqlite3 database.
    """
    def __init__(self, conn, table, genomeReader):
        """genomeReader maybe None if splice sites are not desired """
        super(EvidencePslDbFactory, self).__init__(genomeReader)
        self.pslDbTable = PslDbTable(conn, table)

    def overlappingGen(self, chrom, start, end, qStrand=None):
        """Generator get overlapping alignments as TranscriptFeatures.
        """
        for psl in self.pslDbTable.getTRangeOverlap(chrom, start, end):
            if psl.getQStrand() == qStrand:
                yield self.fromPsl(psl)
