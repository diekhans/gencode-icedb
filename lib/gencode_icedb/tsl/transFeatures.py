"""
Features of a transcript annotation or alignment.
"""
from __future__ import print_function
from pycbio.hgdata import dnaOps
from pycbio.hgdata.frame import Frame


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
        return "{}-{} rna={}-{}".format(self.chromStart, self.chromEnd,
                                        self.rnaStart, self.rnaEnd)

    def _baseStr(self):
        "get base object string, easier than using super"
        return TransFeature.__str__(self)

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

    @property
    def chromSize(self):
        if self.chromStart is None:
            return 0
        else:
            return self.chromEnd - self.chromStart


class AlignBlockFeature(TransFeature):
    """Ungapped alignment block, or unaligned region.  Insertions will have None for the
    other side of the pairwise alignment"""
    __slots__ = ()

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        return AlignBlockFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3])


class Utr5RegionFeature(TransFeature):
    "A un-gapped 5'UTR region in an exon"
    __slots__ = ()

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd):
        super(Utr5RegionFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)

    def __str__(self):
        return "5'UTR {}".format(self._baseStr())

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        return Utr5RegionFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3])


class CdsRegionFeature(TransFeature):
    "A un-gapped CDS region in an exon"
    __slots__ = ("frame")

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd, frame):
        assert(isinstance(frame, Frame))
        super(CdsRegionFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)
        self.frame = frame

    def __str__(self):
        return "CDS {} {}".format(self._baseStr(), self.frame)

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        rcFrame = self.frame + (rcCoords[3] - rcCoords[2])
        return CdsRegionFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3], rcFrame)


class Utr3RegionFeature(TransFeature):
    "A un-gapped 3'UTR region in an exon"
    __slots__ = ()

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd):
        super(Utr3RegionFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)

    def __str__(self):
        return "3'UTR {}".format(self._baseStr())

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        return Utr3RegionFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3])


class ExonFeature(TransFeature):
    """exon with target gaps closed."""
    __slots__ = ("codingFeatures", "alignFeatures")

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd):
        super(ExonFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)
        self.codingFeatures, self.alignFeatures = None, None

    def __str__(self):
        return "exon {}".format(self._baseStr())

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        rcExon = ExonFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3])
        rcExon.codingFeatures = _reverseComplementChildren(rcExon, self.codingFeatures)
        rcExon.alignFeatures = _reverseComplementChildren(rcExon, self.alignFeatures)
        return rcExon

    def rnaOverlaps(self, exon2):
        "does RNA range overlap another exon?"
        return (self.rnaStart < exon2.rnaEnd) and (self.rnaEnd > exon2.rnaStart)


class IntronFeature(TransFeature):
    """intron from annotation or alignment, splice junction information maybe
    None. alignFeatures are for query insertions in intron, and maybe None"""
    __slots__ = ("donorSeq", "acceptorSeq", "spliceSites", "alignFeatures")

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd, donorSeq, acceptorSeq, spliceSites):
        super(IntronFeature, self).__init__(parent, chromStart, chromEnd, rnaStart, rnaEnd)
        self.donorSeq, self.acceptorSeq, self.spliceSites = donorSeq, acceptorSeq, spliceSites
        self.alignFeatures = None

    def __str__(self):
        if self.donorSeq is None:
            sjDesc = None
        else:
            sjDesc = "{}...{} ({})".format(self.donorSeq, self.acceptorSeq, self.spliceSites)
        return "intron {} sjBases={}".format(self._baseStr(), sjDesc)

    def reverseComplement(self, rcParent):
        rcCoords = self.reverseCoords()
        rcDonorSeq, rcAcceptorSeq = (None, None) if self.donorSeq is None else (dnaOps.reverseComplement(self.acceptorSeq), dnaOps.reverseComplement(self.donorSeq))
        rcIntron = IntronFeature(rcParent, rcCoords[0], rcCoords[1], rcCoords[2], rcCoords[3],
                                 rcDonorSeq, rcAcceptorSeq, self.spliceSites)
        rcIntron.alignFeatures = _reverseComplementChildren(rcParent, self.alignFeatures)
        return rcIntron

    def rnaIntersect(self, exon2):
        "does RNA interbase range intersect another intron? (overlap allowing zero length)"
        return (self.rnaStart <= exon2.rnaEnd) and (self.rnaEnd >= exon2.rnaStart)


class TranscriptFeatures(TransFeature):
    """
    Set of features for a transcript derived from an alignment or annotation,
    features are kept in chrom strand order.
    """
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

    @property
    def alignedBases(self):
        alignedCnt = 0
        for feature in self:
            if isinstance(feature, ExonFeature):
                alignedCnt += (feature.rnaEnd - feature.rnaStart) - feature.rnaInsertCnt
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
