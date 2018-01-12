"""
Features of a transcript annotation or alignment.
"""
from __future__ import print_function
from pycbio.hgdata import dnaOps
from pycbio.hgdata.frame import Frame
from pycbio.hgdata.bed import Bed
from gencode_icedb.general.spliceJuncs import SpliceJuncs, spliceJuncsClassify


def _reverseComplementChildren(rcParent, features):
    "reverse complement a list of child TransFeatures, return None if features is None"
    if features is None:
        return None
    rcFeatures = []
    for i in range(len(features) - 1, -1, -1):
        rcFeatures.append(features[i].reverseComplement(rcParent))
    return tuple(rcFeatures)


class TransFeature(object):
    """A feature of annotation or alignment to the genome.  The RNA range is
    in genomic coordinates and it's length may not match the length of the
    feature in gapped feature types, due to gap closing.  For introns, are the
    interbase coordinates of the intron, which might not be equal if there are
    unaligned bases in the intron.  It is possible for either the chrom or RNA
    range to be None."""
    __slots__ = ("parent", "chrom", "rna")
    name = None  # overridden by derived case

    def __init__(self, parent, chrom, rna):
        assert (chrom is None) or ((chrom.start is not None) and (chrom.end is not None))
        assert (rna is None) or ((rna.start is not None) and (rna.end is not None))
        self.parent = parent
        self.chrom = chrom
        self.rna = rna

    def __str__(self):
        "default str(), uses name"
        return "{} {}".format(self.name, self.coordsStr())

    def coordsStr(self):
        "get coordinates string"
        cBounds = (self.chrom.start, self.chrom.end) if self.chrom is not None else (None, None)
        rBounds = (self.rna.start, self.rna.end) if self.rna is not None else (None, None)
        return "{}-{} rna={}-{}".format(cBounds[0], cBounds[1], rBounds[0], rBounds[1])

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        return (str(self),)

    def _getChildrenStrTree(self, features):
        "build tuple for a list of children"
        r = []
        for feat in features:
            r.append(feat.toStrTree())
        return tuple(r)

    @property
    def transcript(self):
        "get the transcript feature by walking the parents"
        if isinstance(self, TranscriptFeatures):
            return self
        elif self.parent is None:
            raise TypeError("no TranscriptFeatures ancestor found: {}".format(self))
        else:
            return self.parent.transcript

    def reverseComplement(self, rcParent):
        "default reverseComplement when __init__ takes just coordinates"
        return self.__class__(rcParent, self.chrom.reverse(), self.rna.reverse())


class AlignedFeature(TransFeature):
    """Ungapped alignment block."""
    name = "aln"
    __slots__ = ()


class ChromInsertFeature(TransFeature):
    """Chromosome insert region (or RNA deletion)."""
    name = "cins"
    __slots__ = ()

    def __init__(self, parent, chrom):
        super(ChromInsertFeature, self).__init__(parent, chrom, None)

    def reverseComplement(self, rcParent):
        return ChromInsertFeature(rcParent, self.chrom.reverse())


class RnaInsertFeature(TransFeature):
    """RNA insert region (or chrom deletion)."""
    name = "rins"
    __slots__ = ()

    def __init__(self, parent, rna):
        super(RnaInsertFeature, self).__init__(parent, None, rna)

    def reverseComplement(self, rcParent):
        return RnaInsertFeature(rcParent, self.rna.reverse())


class Utr5RegionFeature(TransFeature):
    "A un-gapped 5'UTR region in an exon"
    name = "5'UTR"
    __slots__ = ()


class CdsRegionFeature(TransFeature):
    "A un-gapped CDS region in an exon"
    name = "CDS"
    __slots__ = ("frame")

    def __init__(self, parent, chrom, rna, frame):
        assert(isinstance(frame, Frame))
        super(CdsRegionFeature, self).__init__(parent, chrom, rna)
        self.frame = frame

    def __str__(self):
        return "{} {} {}".format(self.name, self.coordsStr(), self.frame)

    def reverseComplement(self, rcParent):
        rcRna = self.rna.reverse()
        rcFrame = self.frame + len(rcRna)
        return CdsRegionFeature(rcParent, self.chrom.reverse(), rcRna, rcFrame)


class Utr3RegionFeature(TransFeature):
    "A un-gapped 3'UTR region in an exon"
    name = "3'UTR"
    __slots__ = ()


class NonCodingRegionFeature(TransFeature):
    "An un-gapped non-coding region in an exon"
    name = "NC"
    __slots__ = ()


class ExonFeature(TransFeature):
    """exon with target gaps closed."""
    name = "exon"
    __slots__ = ("rnaFeatures", "alignFeatures")

    def __init__(self, parent, chrom, rna):
        super(ExonFeature, self).__init__(parent, chrom, rna)
        self.rnaFeatures = self.alignFeatures = None

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.rnaFeatures is not None:
            r.append(self._getChildrenStrTree(self.rnaFeatures))
        if self.alignFeatures is not None:
            r.append(self._getChildrenStrTree(self.alignFeatures))
        return tuple(r)

    def reverseComplement(self, rcParent):
        rcExon = ExonFeature(rcParent, self.chrom.reverse(), self.rna.reverse())
        rcExon.rnaFeatures = _reverseComplementChildren(rcExon, self.rnaFeatures)
        rcExon.alignFeatures = _reverseComplementChildren(rcExon, self.alignFeatures)
        return rcExon

    def rnaOverlaps(self, exon2):
        "does RNA range overlap another exon?"
        return (self.rna.start < exon2.rna.end) and (self.rna.end > exon2.rna.start)

    def __countAlignedBases(self):
        alignedCnt = 0
        for blk in self.alignFeatures:
            if isinstance(blk, AlignedFeature):
                alignedCnt += len(blk.rna)
        return alignedCnt

    @property
    def alignedBases(self):
        """if there are alignment subfeatures, return the actual number of
        aligned bases, otherwise, the RNA for an annotation"""
        if len(self.alignFeatures) > 0:
            return self.__countAlignedBases()
        else:
            return len(self.rna)


class IntronFeature(TransFeature):
    """Intron from annotation or alignment, splice junction information maybe
    None. The alignFeatures field are for query insertions in intron, and maybe None.
    The donorSeq and acceptorSeq values are in the direction of transcription and
    are lower-case if splice junction motif is unknown and upper case if known.
    The spliceJuncs field is a SpliceJuncs object or None.
    """
    name = "intron"
    __slots__ = ("donorSeq", "acceptorSeq", "spliceJuncs", "alignFeatures")

    def __init__(self, parent, chrom, rna, donorSeq, acceptorSeq):
        super(IntronFeature, self).__init__(parent, chrom, rna)
        self.donorSeq = self.acceptorSeq = self.spliceJuncs = None
        if donorSeq is not None:
            self.spliceJuncs = spliceJuncsClassify(donorSeq, acceptorSeq)
            if self.spliceJuncs == SpliceJuncs.unknown:
                self.donorSeq, self.acceptorSeq = donorSeq.lower(), acceptorSeq.lower()
            else:
                self.donorSeq, self.acceptorSeq = donorSeq.upper(), acceptorSeq.upper()
        self.alignFeatures = None

    def __str__(self):
        if self.donorSeq is None:
            sjDesc = None
        else:
            sjDesc = "{}...{} ({})".format(self.donorSeq, self.acceptorSeq, str(self.spliceJuncs))
        return "intron {} sjBases={}".format(self.coordsStr(), sjDesc)

    @property
    def sjBases(self):
        """Get splice junction patterns, or None if not available"""
        if self.donorSeq is None:
            return None
        else:
            return "{}/{}".format(self.donorSeq, self.acceptorSeq)

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.alignFeatures is not None:
            r.append(self._getChildrenStrTree(self.alignFeatures))
        return tuple(r)

    def reverseComplement(self, rcParent):
        rcIntron = IntronFeature(rcParent, self.chrom.reverse(), self.rna.reverse(),
                                 self.donorSeq, self.acceptorSeq)
        rcIntron.alignFeatures = _reverseComplementChildren(rcParent, self.alignFeatures)
        return rcIntron

    def rnaIntersect(self, exon2):
        "does RNA interbase range intersect another intron? (overlap allowing zero length)"
        return (self.rna.start <= exon2.rna.end) and (self.rna.end >= exon2.rna.start)


class TranscriptFeatures(TransFeature):
    """
    Set of features for a transcript derived from an alignment or annotation,
    features are kept in chromosome order (positive strand).
    """
    name = "trans"
    __slots__ = ("chrom", "rna", "cdsChromStart", "cdsChromEnd", "features")

    def __init__(self, chrom, rna, cdsChromStart=None, cdsChromEnd=None):
        super(TranscriptFeatures, self).__init__(None, chrom, rna)
        self.chrom = chrom
        self.rna = rna
        self.cdsChromStart = cdsChromStart
        self.cdsChromEnd = cdsChromEnd
        self.features = None

    def __str__(self):
        return "t={}/{}, rna={}/{} {}".format(str(self.chrom), self.chrom.strand,
                                              str(self.rna), self.rna.strand, self.rna.size)

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.features is not None:
            r.append(self._getChildrenStrTree(self.features))
        return tuple(r)

    def toBed(self, itemRgb=""):
        """convert transcript and CDS to a Bed object"""
        blocks = []
        for feat in self.features:
            if isinstance(feat, ExonFeature):
                blocks.append(Bed.Block(feat.chrom.start, feat.chrom.end))
        cdsStart = self.cdsChromStart if self.cdsChromStart is not None else self.chrom.end
        cdsEnd = self.cdsChromEnd if self.cdsChromEnd is not None else self.chrom.end
        return Bed(self.chrom.name, self.chrom.start, self.chrom.end,
                   self.rna.name, 0, self.rna.strand, cdsStart, cdsEnd,
                   itemRgb, blocks)

    @property
    def alignedBases(self):
        alignedCnt = 0
        for feature in self.features:
            if isinstance(feature, ExonFeature):
                alignedCnt += feature.alignedBases
        return alignedCnt

    def reverseComplement(self):
        "return a new TranscriptFeatures object that is reverse complemented"
        if self.chrom.size is None:
            raise Exception("can't reverse-complement transcript without chrom.size: {}".format(self))
        rcCdsChromStart, rcCdsChromEnd = None, None if self.cdsChromStart is None else dnaOps.reverseCoords(self.cdsChromStart, self.cdsChromEnd, self.chrom.size)

        rcTrans = TranscriptFeatures(self.chrom.reverse(), self.rna.reverse(),
                                     rcCdsChromStart, rcCdsChromEnd)
        rcTrans.features = _reverseComplementChildren(rcTrans, self.features)
        return rcTrans
