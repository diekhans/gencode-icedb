"""
Features of a transcript annotation or alignment.
"""


class TransFeature(object):
    "a feature of annotation or alignment to the genome"
    __slots__ = ("parent", "chromStart", "chromEnd")

    def __init__(self, parent, chromStart, chromEnd):
        self.parent, self.chromStart, self.chromEnd = parent, chromStart, chromEnd

    @property
    def size(self):
        return self.chromEnd - self.chromStart


class ExonFeature(TransFeature):
    """exon with target gaps closed.  The query range is in genomic coordinates and
    it's length may not match the length of the feature, due to gap closing."""
    __slots__ = ("rnaStart", "rnaEnd", "rnaInsertCnt", "chromInsertCnt")

    def __init__(self, parent, chromStart, chromEnd, rnaStart, rnaEnd, rnaInsertCnt, chromInsertCnt):
        super(ExonFeature, self).__init__(parent, chromStart, chromEnd)
        self.rnaStart, self.rnaEnd = rnaStart, rnaEnd
        self.rnaInsertCnt, self.chromInsertCnt = rnaInsertCnt, chromInsertCnt

    def __str__(self):
        return "exon {}-{} q=({}-{}) qIns={} tIns={}".format(self.chromStart, self.chromEnd,
                                                             self.rnaStart, self.rnaEnd,
                                                             self.rnaInsertCnt, self.chromInsertCnt)

    def rnaOverlaps(self, exon2):
        "does RNA rnage overlap?"
        return (self.rnaStart < exon2.rnaEnd) and (self.rnaEnd > exon2.rnaStart)


class IntronFeature(TransFeature):
    """intron from annotation or alignment, splice junction information maybe
    None"""
    __slots__ = ("rnaDeleteCnt", "donorSeq", "acceptorSeq", "spliceSites")

    def __init__(self, parent, chromStart, chromEnd, rnaDeleteCnt, donorSeq, acceptorSeq, spliceSites):
        super(IntronFeature, self).__init__(parent, chromStart, chromEnd)
        self.rnaDeleteCnt, self.donorSeq, self.acceptorSeq, self.spliceSites = rnaDeleteCnt, donorSeq, acceptorSeq, spliceSites

    def __str__(self):
        if self.donorSeq is None:
            sjDesc = None
        else:
            sjDesc = "{}...{} ({})".format(self.donorSeq, self.acceptorSeq, self.spliceSites)
        return "intron {}-{} qDel={} sjBases={}".format(self.chromStart, self.chromEnd, self.rnaDeleteCnt, sjDesc)


class TranscriptFeatures(list):
    """
    Set of features for a transcript derived from an alignment or annotation,
    features are kept in genomic order.
    """

    def __init__(self, chrom, chromStrand, chromStart, chromEnd, chromSize,
                 rnaName, rnaStrand, rnaStart, rnaEnd, rnaSize,
                 features, cdsChromStart=None, cdsChromEnd=None):
        self.chrom, self.chromStrand, self.chromStart, self.chromEnd, self.chromSize = chrom, chromStrand, chromStart, chromEnd, chromSize
        self.rnaName, self.rnaStrand, self.rnaStart, self.rnaEnd, self.rnaSize = rnaName, rnaStrand, rnaStart, rnaEnd, rnaSize
        self.cdsChromStart, self.cdsChromEnd = cdsChromStart, cdsChromEnd
        self.extend(features)

    def __str__(self):
        return "t={}:{}-{}/{}, q={}:{}-{}/{} {}".format(self.chrom, self.chromStart, self.chromEnd, self.chromStrand,
                                                      self.rnaName, self.rnaStart, self.rnaEnd, self.rnaStrand, self.rnaSize)

    @property
    def alignedBases(self):
        alignedCnt = 0
        for feature in self:
            if isinstance(feature, ExonFeature):
                alignedCnt += (feature.rnaEnd - feature.rnaStart) - feature.rnaInsertCnt
        return alignedCnt
