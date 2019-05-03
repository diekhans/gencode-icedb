"""
Convert alignments (PSL or BAM) to features, closing and tracking gaps.
"""
from collections import namedtuple
from pycbio.hgdata.coords import Coords
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, AlignedFeature, ChromInsertFeature, RnaInsertFeature

class _EvidenceFactoryCreator(object):
    """Class to generate TranscriptFeatures from a polymorphic representation
    of the alignment.  This representation uses a duck-typed subset of the
    pycbio Psl and PslBlock fields."""
    def __init__(self, genomeReader):
        """genomeReader maybe None if splice sites are not desired """
        self.genomeReader = genomeReader

    def _buildFeatures(self, aln, trans):
        iBlkStart = 0
        features = []
        while iBlkStart < len(aln.blocks):
            iBlkEnd = self._findExonEnd(aln, iBlkStart)
            self._addExon(aln, iBlkStart, iBlkEnd, trans, features)
            if iBlkEnd < len(aln.blocks):
                self._addIntron(aln, iBlkEnd, trans, features)
            iBlkStart = iBlkEnd
        return features

    def _tGapSize(self, aln, iBlk):
        "size of gap before the block"
        return aln.blocks[iBlk].tStart - aln.blocks[iBlk - 1].tEnd

    def _findExonEnd(self, aln, iBlkStart):
        "finds half-open end of blocks covering current exon"
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(aln.blocks)) and (self._tGapSize(aln, iBlkEnd) < minIntronSize):
            iBlkEnd += 1
        return iBlkEnd

    def _addFirstUnaligned(self, aln, exon, alignFeatures):
        blk = aln.blocks[0]
        if blk.qStart > 0:
            alignFeatures.append(RnaInsertFeature(exon, len(alignFeatures),
                                                  exon.rna.subrange(0, blk.qStart, aln.qStrand)))

    def _addLastUnaligned(self, aln, exon, alignFeatures):
        blk = aln.blocks[-1]
        if blk.qEnd < aln.qSize:
            alignFeatures.append(RnaInsertFeature(exon, len(alignFeatures),
                                                  exon.rna.subrange(blk.qEnd, aln.qSize, aln.qStrand)))

    def _addAlignedFeature(self, aln, iBlk, exon, alignFeatures):
        blk = aln.blocks[iBlk]
        alignFeatures.append(AlignedFeature(exon, len(alignFeatures),
                                            exon.chrom.subrange(blk.tStart, blk.tEnd, aln.tStrand),
                                            exon.rna.subrange(blk.qStart, blk.qEnd, aln.qStrand)))

    def _addUnalignedFeatures(self, aln, iBlk, feat, alignFeatures):
        prevBlk = aln.blocks[iBlk - 1]
        blk = aln.blocks[iBlk]
        if blk.qStart > prevBlk.qEnd:
            alignFeatures.append(RnaInsertFeature(feat, len(alignFeatures),
                                                  feat.rna.subrange(prevBlk.qEnd, blk.qStart, aln.qStrand)))
        if blk.tStart > prevBlk.tEnd:
            alignFeatures.append(ChromInsertFeature(feat, len(alignFeatures),
                                                    feat.chrom.subrange(prevBlk.tEnd, blk.tStart, aln.tStrand)))

    def _addAlignFeatures(self, aln, iBlkStart, iBlkEnd, exon):
        alignFeatures = []
        if iBlkStart == 0:
            self._addFirstUnaligned(aln, exon, alignFeatures)
        for iBlk in range(iBlkStart, iBlkEnd):
            if iBlk > iBlkStart:
                self._addUnalignedFeatures(aln, iBlk, exon, alignFeatures)
            self._addAlignedFeature(aln, iBlk, exon, alignFeatures)  # after since unaligned is before block
        if iBlkEnd == len(aln.blocks):
            self._addLastUnaligned(aln, exon, alignFeatures)
        exon.alignFeatures = tuple(alignFeatures)

    def _addExon(self, aln, iBlkStart, iBlkEnd, trans, features):
        # include either initial or terminal unaligned parts of RNA if at start or end
        qStart = aln.blocks[iBlkStart].qStart if iBlkStart > 0 else 0
        qEnd = aln.blocks[iBlkEnd - 1].qEnd if iBlkEnd < len(aln.blocks) else aln.qSize

        exon = ExonFeature(trans, len(features),
                           trans.chrom.subrange(aln.blocks[iBlkStart].tStart, aln.blocks[iBlkEnd - 1].tEnd, aln.tStrand),
                           trans.rna.subrange(qStart, qEnd, aln.qStrand))
        features.append(exon)
        self._addAlignFeatures(aln, iBlkStart, iBlkEnd, exon)

    def _getSpliceSites(self, aln, iBlkNext, trans):
        if self.genomeReader is None:
            return (None, None)
        else:
            # this handles 3' ESTs
            coords = Coords(aln.tName, aln.blocks[iBlkNext - 1].tEnd, aln.blocks[iBlkNext].tStart, aln.tStrand, aln.tSize)
            if coords.strand == '-':
                coords = coords.reverse()
            return spliceJuncsGetSeqs(self.genomeReader, coords.name, coords.start, coords.end, trans.transcriptionStrand)

    def _addIntron(self, aln, iBlkNext, trans, features):
        donorSeq, acceptorSeq = self._getSpliceSites(aln, iBlkNext, trans)
        intron = IntronFeature(trans, len(features),
                               trans.chrom.subrange(aln.blocks[iBlkNext - 1].tEnd, aln.blocks[iBlkNext].tStart, aln.tStrand),
                               trans.rna.subrange(aln.blocks[iBlkNext - 1].qEnd, aln.blocks[iBlkNext].qStart, aln.qStrand),
                               donorSeq, acceptorSeq)
        features.append(intron)
        alignFeatures = []
        self._addUnalignedFeatures(aln, iBlkNext, intron, alignFeatures)
        intron.alignFeatures = tuple(alignFeatures)

    def fromAlignment(self, aln, attrs, orientChrom):
        """Convert a abstract alginment to a TranscriptFeatures object.  If
        orientChrom is True, then always ensure chrom strand is positive.
        """
        # After tslGetUcscRnaAligns adjustments, 3' EST are `+-'
        # for positive stand genes and `--' for negative strand genes.
        # Thus transcriptionStrand is always the qStrand regardless of tStrand.
        transcriptionStrand = aln.qStrand

        if orientChrom and aln.tStrand == '-':
            aln = aln.reverseComplement()
        chrom = Coords(aln.tName, aln.tStart, aln.tEnd, '+', aln.tSize)
        if aln.tStrand == '-':
            chrom = chrom.reverse()
        rna = Coords(aln.qName, 0, aln.qSize, '+', aln.qSize)
        if aln.qStrand == '-':
            rna = rna.reverse()

        trans = TranscriptFeatures(chrom, rna, transcriptionStrand=transcriptionStrand, attrs=attrs)
        trans.features = tuple(self._buildFeatures(aln, trans))
        return trans


class EvidencePslFactory(object):
    """
    Factory to create evidence features from PSLs.
    """
    def __init__(self, genomeReader=None):
        """genomeReader maybe None if splice sites are not desired """
        self.creator = _EvidenceFactoryCreator(genomeReader)

    def fromPsl(self, psl, attrs=None, orientChrom=True):
        """Convert a psl to a TranscriptFeatures object.  If orientChrom is
        True, then always ensure chrom strand is positive.
        """
        return self.creator.fromAlignment(psl, attrs, orientChrom)


BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CONSUMES_QUERY_OPS = frozenset([BAM_CMATCH, BAM_CINS, BAM_CSOFT_CLIP, BAM_CEQUAL, BAM_CDIFF])
BAM_CONSUMES_REF_OPS = frozenset([BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP, BAM_CEQUAL, BAM_CDIFF])
BAM_CLIP_OPS = frozenset([BAM_CSOFT_CLIP, BAM_CHARD_CLIP])


class SamAlign(object):
    """Object that implements a PSL-like interface for a pysam AlignedSegment object"""

    class Block(namedtuple("Block",
                           ("qStart", "tStart", "size"))):
        "Ungapped block gap extracted from the align,ent"
        __slots__ = ()

        @property
        def qEnd(self):
            return self.qStart + self.size

        @property
        def tEnd(self):
            return self.tStart + self.size

        def __str__(self):
            return "{}..{} {}..{} [{}]".format(self.qStart, self.qEnd, self.tStart, self.tEnd, self.size)

    def __init__(self, samfh, alnseg):
        self.samfh = samfh
        self.alnseg = alnseg
        self.blocks = self._buildBlocks(alnseg)

    def _buildBlocks(self, alnseg):
        """Intermediate transform to list of ungapped blocks.  This transform
        eliminates a lot of passing around book-keeping variables"""
        qNext = 0
        tNext = alnseg.reference_start
        blocks = []
        for op, size in alnseg.cigartuples:
            qStart = qNext if op in BAM_CONSUMES_QUERY_OPS else None
            tStart = tNext if op in BAM_CONSUMES_REF_OPS else None
            if (qStart is not None) and (tStart is not None) and (op not in BAM_CLIP_OPS):
                blocks.append(self.Block(qStart, tStart, size))
            if qStart is not None:
                qNext += size
            if tStart is not None:
                tNext += size
        return blocks

    @property
    def qName(self):
        return self.alnseg.query_name

    @property
    def qStart(self):
        return self.alnseg. query_alignment_start

    @property
    def qEnd(self):
        return self.alnseg.query_alignment_end

    @property
    def qSize(self):
        size = self.alnseg.query_length
        if size == 0:
            size = self.alnseg.infer_query_length()
        return size

    @property
    def qStrand(self):
        return '-' if self.alnseg.is_reverse else '+'

    @property
    def tName(self):
        return self.alnseg.reference_name

    @property
    def tStart(self):
        return self.blocks[0].tStart

    @property
    def tEnd(self):
        return self.blocks[-1].tEnd

    @property
    def tSize(self):
        return self.samfh.get_reference_length(self.alnseg.reference_name)

    @property
    def tStrand(self):
        return '+'


class EvidenceSamFactory(object):
    """
    Factory to create evidence features from SAM/BAM/CRAM records read by pysam.
    """

    def __init__(self, genomeReader=None):
        """genomeReader maybe None if splice sites are not desired """
        self.creator = _EvidenceFactoryCreator(genomeReader)

    def fromSam(self, samfh, alnseg, attrs=None, orientChrom=True):
        """Convert a BAM/SAM/CRAM record to a TranscriptFeatures object.  If
        orientChrom is True, then always ensure chrom strand is positive.
        """

        return self.creator.fromAlignment(SamAlign(samfh, alnseg), attrs, orientChrom)
