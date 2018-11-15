"""
Convert alignments (PSL or BAM) to features, closing and tracking gaps.
"""
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

    def _buildFeatures(self, psl, trans):
        iBlkStart = 0
        features = []
        while iBlkStart < len(psl.blocks):
            iBlkEnd = self._findExonEnd(psl, iBlkStart)
            self._addExon(psl, iBlkStart, iBlkEnd, trans, features)
            if iBlkEnd < len(psl.blocks):
                self._addIntron(psl, iBlkEnd, trans, features)
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
        blk = psl.blocks[0]
        if blk.qStart > 0:
            alignFeatures.append(RnaInsertFeature(exon, len(alignFeatures),
                                                  exon.rna.subrange(0, blk.qStart, psl.qStrand)))

    def _addLastUnaligned(self, psl, exon, alignFeatures):
        blk = psl.blocks[-1]
        if blk.qEnd < psl.qSize:
            alignFeatures.append(RnaInsertFeature(exon, len(alignFeatures),
                                                  exon.rna.subrange(blk.qEnd, psl.qSize, psl.qStrand)))

    def _addAlignedFeature(self, psl, iBlk, exon, alignFeatures):
        blk = psl.blocks[iBlk]
        alignFeatures.append(AlignedFeature(exon, len(alignFeatures),
                                            exon.chrom.subrange(blk.tStart, blk.tEnd, psl.tStrand),
                                            exon.rna.subrange(blk.qStart, blk.qEnd, psl.qStrand)))

    def _addUnalignedFeatures(self, psl, iBlk, feat, alignFeatures):
        prevBlk = psl.blocks[iBlk - 1]
        blk = psl.blocks[iBlk]
        if blk.qStart > prevBlk.qEnd:
            alignFeatures.append(RnaInsertFeature(feat, len(alignFeatures),
                                                  feat.rna.subrange(prevBlk.qEnd, blk.qStart, psl.qStrand)))
        if blk.tStart > prevBlk.tEnd:
            alignFeatures.append(ChromInsertFeature(feat, len(alignFeatures),
                                                    feat.chrom.subrange(prevBlk.tEnd, blk.tStart, psl.tStrand)))

    def _addAlignFeatures(self, psl, iBlkStart, iBlkEnd, exon):
        alignFeatures = []
        if iBlkStart == 0:
            self._addFirstUnaligned(psl, exon, alignFeatures)
        for iBlk in range(iBlkStart, iBlkEnd):
            if iBlk > iBlkStart:
                self._addUnalignedFeatures(psl, iBlk, exon, alignFeatures)
            self._addAlignedFeature(psl, iBlk, exon, alignFeatures)  # after since unaligned is before block
        if iBlkEnd == psl.blockCount:
            self._addLastUnaligned(psl, exon, alignFeatures)
        exon.alignFeatures = tuple(alignFeatures)

    def _addExon(self, psl, iBlkStart, iBlkEnd, trans, features):
        # include either initial or terminal unaligned parts of RNA if at start or end
        qStart = psl.blocks[iBlkStart].qStart if iBlkStart > 0 else 0
        qEnd = psl.blocks[iBlkEnd - 1].qEnd if iBlkEnd < psl.blockCount else psl.qSize

        exon = ExonFeature(trans, len(features),
                           trans.chrom.subrange(psl.blocks[iBlkStart].tStart, psl.blocks[iBlkEnd - 1].tEnd, psl.tStrand),
                           trans.rna.subrange(qStart, qEnd, psl.qStrand))
        features.append(exon)
        self._addAlignFeatures(psl, iBlkStart, iBlkEnd, exon)

    def _getSpliceSites(self, psl, iBlkNext, trans):
        if self.genomeReader is None:
            return (None, None)
        else:
            # this handles 3' ESTs
            coords = Coords(psl.tName, psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart, psl.tStrand, psl.tSize)
            if coords.strand == '-':
                coords = coords.reverse()
            return spliceJuncsGetSeqs(self.genomeReader, coords.name, coords.start, coords.end, trans.transcriptionStrand)

    def _addIntron(self, psl, iBlkNext, trans, features):
        donorSeq, acceptorSeq = self._getSpliceSites(psl, iBlkNext, trans)
        intron = IntronFeature(trans, len(features),
                               trans.chrom.subrange(psl.blocks[iBlkNext - 1].tEnd, psl.blocks[iBlkNext].tStart, psl.tStrand),
                               trans.rna.subrange(psl.blocks[iBlkNext - 1].qEnd, psl.blocks[iBlkNext].qStart, psl.qStrand),
                               donorSeq, acceptorSeq)
        features.append(intron)
        alignFeatures = []
        self._addUnalignedFeatures(psl, iBlkNext, intron, alignFeatures)
        intron.alignFeatures = tuple(alignFeatures)

    def fromPsl(self, psl, attrs=None, orientChrom=True):
        """Convert a psl to a TranscriptFeatures object.  If orientChrom is
        True, then always ensure chrom strand is positive.
        """
        # After tslGetUcscRnaAligns adjustments, 3' EST are `+-'
        # for positive stand genes and `--' for negative strand genes.
        # Thus transcriptionStrand is always the qStrand regardless of tStrand.
        transcriptionStrand = psl.qStrand

        if orientChrom and psl.tStrand == '-':
            psl = psl.reverseComplement()
        chrom = Coords(psl.tName, psl.tStart, psl.tEnd, '+', psl.tSize)
        if psl.tStrand == '-':
            chrom = chrom.reverse()
        rna = Coords(psl.qName, 0, psl.qSize, '+', psl.qSize)
        if psl.qStrand == '-':
            rna = rna.reverse()

        trans = TranscriptFeatures(chrom, rna, transcriptionStrand=transcriptionStrand, attrs=attrs)
        trans.features = tuple(self._buildFeatures(psl, trans))
        return trans


class EvidenceSamFactory(object):
    """
    Factory to create evidence features from SAM/BAM/CRAM records.
    """
    pass
