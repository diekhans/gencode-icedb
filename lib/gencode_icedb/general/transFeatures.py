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
    for i in xrange(len(features) - 1, -1, -1):
        rcFeatures.append(features[i].reverseComplement(rcParent))
    return tuple(rcFeatures)


class TransFeature(object):
    """A feature of annotation or alignment to the genome.  The RNA range is
    in genomic coordinates and it's length may not match the length of the
    feature in gapped feature types, due to gap closing.  For introns, are the
    interbase coordinates of the intron, which might not be equal if there are
    unaligned bases in the intron.  It is possible for either the chrom or RNA
    range to be None."""
    __slots__ = ("parent", "chromStart", "chromEnd", "rnaStart", "rnaEnd")

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd):
        self.parent = parent
        self.chromStart, self.chromEnd = chromStart, chromEnd
        self.rnaStart, self.rnaEnd = rnaStart, rnaEnd

    def __str__(self):
        "default str(), uses name"
        return "{} {}".format(self.name, self.coordsStr())

    def coordsStr(self):
        "get coordinates string"
        return "{}-{} rna={}-{}".format(self.chromStart, self.chromEnd,
                                        self.rnaStart, self.rnaEnd)

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

    def reverseCoords(self):
        "get tuple of reverse coordinates"
        trans = self.transcript
        if self.chromStart is not None:
            rcChromStart, rcChromEnd = dnaOps.reverseCoords(self.chromStart, self.chromEnd, trans.chromSize)
        else:
            rcChromStart, rcChromEnd = None, None
        if self.rnaStart is not None:
            rcRnaStart, rcRnaEnd = dnaOps.reverseCoords(self.rnaStart, self.rnaEnd, trans.rnaSize)
        else:
            rcRnaStart, rcRnaEnd = None, None
        return (rcChromStart, rcChromEnd, rcRnaStart, rcRnaEnd)

    def reverseComplement(self, rcParent):
        "default reverseComplement when __init__ takes just coordinates"
        rcCoords = self.reverseCoords()
        return self.__class__(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3])

    @property
    def chromLength(self):
        if self.chromStart is None:
            return 0
        else:
            return self.chromEnd - self.chromStart

    @property
    def rnaLength(self):
        if self.rnaStart is None:
            return 0
        else:
            return self.rnaEnd - self.rnaStart


class AlignedFeature(TransFeature):
    """Ungapped alignment block."""
    name = "aln"
    __slots__ = ()


class ChromInsertFeature(TransFeature):
    """Chromosome insert region (or RNA deletion)."""
    name = "cins"
    __slots__ = ()

    def __init__(self, parent, chromStart, chromEnd):
        super(ChromInsertFeature, self).__init__(parent, chromStart, chromEnd, None, None)

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        return ChromInsertFeature(rcParent, rcCoords[0], rcCoords[1])


class RnaInsertFeature(TransFeature):
    """RNA insert region (or chrom deletion)."""
    name = "rins"
    __slots__ = ()

    def __init__(self, parent, rnaStart, rnaEnd):
        super(RnaInsertFeature, self).__init__(parent, None, None, rnaStart, rnaEnd)

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        return RnaInsertFeature(rcParent, rcCoords[2], rcCoords[3])


class Utr5RegionFeature(TransFeature):
    "A un-gapped 5'UTR region in an exon"
    name = "5'UTR"
    __slots__ = ()


class CdsRegionFeature(TransFeature):
    "A un-gapped CDS region in an exon"
    name = "CDS"
    __slots__ = ("frame")

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd, frame):
        assert(isinstance(frame, Frame))
        super(CdsRegionFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)
        self.frame = frame

    def __str__(self):
        return "{} {} {}".format(self.name, self.coordsStr(), self.frame)

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        rcFrame = self.frame + (rcCoords[3] - rcCoords[2])
        return CdsRegionFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3], rcFrame)


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

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd):
        super(ExonFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)
        self.rnaFeatures, self.alignFeatures = None, None

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
        rcCoords = self.reverseCoords()
        rcExon = ExonFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3])
        rcExon.rnaFeatures = _reverseComplementChildren(rcExon, self.rnaFeatures)
        rcExon.alignFeatures = _reverseComplementChildren(rcExon, self.alignFeatures)
        return rcExon

    def rnaOverlaps(self, exon2):
        "does RNA range overlap another exon?"
        return (self.rnaStart < exon2.rnaEnd) and (self.rnaEnd > exon2.rnaStart)

    def __countAlignedBases(self):
        alignedCnt = 0
        for blk in self.alignFeatures:
            if isinstance(blk, AlignedFeature):
                alignedCnt += blk.rnaLength
        return alignedCnt

    @property
    def alignedBases(self):
        """if there are alignment subfeatures, return the actual number of
        aligned bases, otherwise, the RNA for an annotation"""
        if len(self.alignFeatures) > 0:
            return self.__countAlignedBases()
        else:
            return self.rnaLength


class IntronFeature(TransFeature):
    """Intron from annotation or alignment, splice junction information maybe
    None. The alignFeatures field are for query insertions in intron, and maybe None.
    The donorSeq and acceptorSeq values are in the direction of transcription and
    are lower-case if splice junction motif is unknown and upper case if known.
    The spliceJuncs field is a SpliceJuncs object or None.
    """
    name = "intron"
    __slots__ = ("donorSeq", "acceptorSeq", "spliceJuncs", "alignFeatures")

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd, donorSeq, acceptorSeq):
        super(IntronFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)
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
            sjDesc = "{}...{} ({})".format(self.donorSeq, self.acceptorSeq, self.spliceJuncs)
        return "intron {} sjBases={}".format(self.coordsStr(), sjDesc)

    def sjBases(self):
        """Get splice junction patterns"""
        if self.donorSeq is None:
            return "??/??"
        else:
            sjDesc = "{}/{}".format(self.donorSeq, self.acceptorSeq)
        
    
    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.alignFeatures is not None:
            r.append(self._getChildrenStrTree(self.alignFeatures))
        return tuple(r)

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        rcIntron = IntronFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3],
                                 self.donorSeq, self.acceptorSeq)
        rcIntron.alignFeatures = _reverseComplementChildren(rcParent, self.alignFeatures)
        return rcIntron

    def rnaIntersect(self, exon2):
        "does RNA interbase range intersect another intron? (overlap allowing zero length)"
        return (self.rnaStart <= exon2.rnaEnd) and (self.rnaEnd >= exon2.rnaStart)


class TranscriptFeatures(TransFeature):
    """
    Set of features for a transcript derived from an alignment or annotation,
    features are kept in chromosome order (positive strand).
    """
    name = "trans"
    __slots__ = ("chrom", "chromStrand", "chromSize", "rnaName", "rnaStrand", "rnaSize",
                 "cdsChromStart", "cdsChromEnd", "features")

    def __init__(self, chrom, chromStrand, chromStart, chromEnd, chromSize,
                 rnaName, rnaStrand, rnaStart, rnaEnd, rnaSize,
                 cdsChromStart=None, cdsChromEnd=None):
        super(TranscriptFeatures, self).__init__(None, chromStart, chromEnd, rnaStart, rnaEnd)
        self.chrom, self.chromStrand, self.chromSize = chrom, chromStrand, chromSize
        self.rnaName, self.rnaStrand, self.rnaSize = rnaName, rnaStrand, rnaSize
        self.cdsChromStart, self.cdsChromEnd = cdsChromStart, cdsChromEnd
        self.features = None

    def __str__(self):
        return "t={}:{}-{}/{}, rna={}:{}-{}/{} {}".format(self.chrom, self.chromStart, self.chromEnd, self.chromStrand,
                                                          self.rnaName, self.rnaStart, self.rnaEnd, self.rnaStrand, self.rnaSize)

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.features is not None:
            r.append(self._getChildrenStrTree(self.features))
        return tuple(r)

    def toBed(self, itemRgb):
        """convert transcript and CDS to a bed vector"""
        blocks = []
        for feat in self.features:
            if isinstance(feat, ExonFeature):
                blocks.append(Bed.Block(feat.chromStart, feat.chromEnd))
        cdsStart = self.cdsChromStart if self.cdsChromStart is not None else self.chromEnd
        cdsEnd = self.cdsChromEnd if self.cdsChromEnd is not None else self.chromEnd
        return Bed(self.chrom, self.chromStart, self.chromEnd,
                   self.rnaName, 0, self.rnaStrand, cdsStart, cdsEnd,
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
        if self.chromSize is None:
            raise Exception("can't reverse-complement transcript without chromSize: {}".format(self))
        rcCoords = self.reverseCoords()
        rcCdsChromStart, rcCdsChromEnd = None, None if self.cdsChromStart is None else dnaOps.reverseCoords(self.cdsChromStart, self.cdsChromEnd, self.chromSize)

        rcTrans = TranscriptFeatures(self.chrom, dnaOps.reverseStrand(self.chromStrand), rcCoords[0], rcCoords[1], self.chromSize,
                                     self.rnaName, dnaOps.reverseStrand(self.rnaStrand), rcCoords[2], rcCoords[3], self.rnaSize,
                                     rcCdsChromStart, rcCdsChromEnd)
        rcTrans.features = _reverseComplementChildren(rcTrans, self.features)
        return rcTrans
