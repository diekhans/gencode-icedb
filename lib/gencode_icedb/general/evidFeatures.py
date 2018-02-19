"""
Convert alignments (PSL) to features, closing and tracking gaps.
"""
from __future__ import print_function
from pycbio.hgdata.coords import Coords
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, AlignedFeature, ChromInsertFeature, RnaInsertFeature

# FIXME: rename alignfeaturs.


class EvidencePslFactory(object):
    """
    Factory to create evidence features from PSLs.
    """
    def __init__(self, genomeReader):
        """genomeReader maybe None if splice sites are not desired """
        self.genomeReader = genomeReader

    def _buildFeatures(self, psl, trans):
        iBlkStart = 0
        features = []
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self._findExonEnd(psl, iBlkStart)
            self._makeExon(psl, iBlkStart, iBlkEnd, trans, features)
            if iBlkEnd < len(psl.blocks):
                self._makeIntron(psl, iBlkEnd, trans, features)
            iBlkStart = iBlkEnd
        return features

    def _tGapSize(self, psl, iBlk):
        "size of gap before the block"
        return psl.blocks[iBlk].tStart - psl.blocks[iBlk - 1].tEnd

    def _findExonEnd(self, psl, iBlkStart):
        "finds half-open end of blocks covering current exon"
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(psl.blocks)) and (self._tGapSize(psl, iBlkEnd) < minIntronSize):
            iBlkEnd += 1
        return iBlkEnd

    def _addFirstUnaligned(self, psl, exon, alignFeatures):
        return  # FIXME: ignored, see issues.org
        if exon.rna.strand == '+':
            qStart, qEnd = 0, psl.qStart
        else:
            qStart, qEnd = 0, psl.qSize - psl.qEnd
        if qStart < qEnd:
            alignFeatures.append(RnaInsertFeature(exon, len(alignFeatures),
                                                  exon.rna.subrange(qStart, qEnd)))

    def _addLastUnaligned(self, psl, exon, alignFeatures):
        return  # FIXME: ignored, see issues.org
        if exon.rna.strand == '+':
            qStart, qEnd = psl.qEnd, psl.qSize
        else:
            qStart, qEnd = 0, psl.qSize - psl.qStart
        if qStart < qEnd:
            alignFeatures.append(RnaInsertFeature(exon, len(alignFeatures),
                                                  exon.rna.subrange(qStart, qEnd)))

    def _addAlignedFeature(self, psl, iBlk, exon, alignFeatures):
        blk = psl.blocks[iBlk]
        alignFeatures.append(AlignedFeature(exon, len(alignFeatures),
                                            exon.chrom.subrange(blk.tStart, blk.tEnd),
                                            exon.rna.subrange(blk.qStart, blk.qEnd)))

    def _addUnalignedFeatures(self, psl, iBlk, feat, alignFeatures):
        prevBlk = psl.blocks[iBlk - 1]
        blk = psl.blocks[iBlk]
        if blk.qStart > prevBlk.qEnd:
            alignFeatures.append(RnaInsertFeature(feat, len(alignFeatures),
                                                  feat.rna.subrange(prevBlk.qEnd, blk.qStart)))
        if blk.tStart > prevBlk.tEnd:
            alignFeatures.append(ChromInsertFeature(feat, len(alignFeatures),
                                                    feat.chrom.subrange(prevBlk.tEnd, blk.tStart)))

    def _addAlignFeatures(self, psl, iBlkStart, iBlkEnd, exon):
        alignFeatures = []
        if iBlkStart == 0:
            self._addFirstUnaligned(psl, exon, alignFeatures)
        for iBlk in range(iBlkStart, iBlkEnd):
            if iBlk > iBlkStart:
                self._addUnalignedFeatures(psl, iBlk, exon, alignFeatures)
            self._addAlignedFeature(psl, iBlk, exon, alignFeatures)  # after since unaligned is before block
        if iBlkEnd == psl.blockCount - 1:
            self._addLastUnaligned(psl, exon, alignFeatures)
        exon.alignFeatures = tuple(alignFeatures)

    def _makeExon(self, psl, iBlkStart, iBlkEnd, trans, features):
        exon = ExonFeature(trans, len(features),
                           trans.chrom.subrange(psl.blocks[iBlkStart].tStart, psl.blocks[iBlkEnd - 1].tEnd),
                           trans.rna.subrange(psl.blocks[iBlkStart].qStart, psl.blocks[iBlkEnd - 1].qEnd))
        features.append(exon)
        self._addAlignFeatures(psl, iBlkStart, iBlkEnd, exon)

    def _getSpliceSites(self, psl, iBlkNext):
        if self.genomeReader is None:
            return (None, None)
        else:
            return spliceJuncsGetSeqs(self.genomeReader, psl.tName, psl.blocks[iBlkNext - 1].tEnd,
                                      psl.blocks[iBlkNext].tStart, psl.getQStrand())

    def _makeIntron(self, psl, iBlkNext, trans, features):
        donorSeq, acceptorSeq = self._getSpliceSites(psl, iBlkNext)
        intron = IntronFeature(trans, len(features),
                               trans.chrom.subrange(psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart),
                               trans.rna.subrange(psl.blocks[iBlkNext - 1].qEnd, psl.blocks[iBlkNext].qStart),
                               donorSeq, acceptorSeq)
        features.append(intron)
        alignFeatures = []
        self._addUnalignedFeatures(psl, iBlkNext, intron, alignFeatures)
        intron.alignFeatures = tuple(alignFeatures)

    def fromPsl(self, psl, attrs=None):
        "convert a psl to an TranscriptFeatures object"
        chrom = Coords(psl.tName, psl.tStart, psl.tEnd, '+', psl.tSize)
        if psl.getTStrand() == '-':
            chrom = chrom.reverse()
        rna = Coords(psl.qName, psl.qStart, psl.qEnd, '+', psl.qSize)
        if psl.getQStrand() == '-':
            rna = rna.reverse()
        trans = TranscriptFeatures(chrom, rna, attrs=attrs)
        trans.features = tuple(self._buildFeatures(psl, trans))
        return trans
