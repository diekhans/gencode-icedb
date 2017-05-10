"""
Features of a transcript annotation or alignment.
"""
from __future__ import print_function
from pycbio.hgdata import dnaOps


class TransFeature(object):
    """A feature of annotation or alignment to the genome.  The RNA range is
    in genomic coordinates and it's length may not match the length of the
    feature, due to gap closing.  For introns, are the interbase coordinates
    of the intron, which might not be equal if there are unaligned bases in
    the intron"""
    __slots__ = ("trans", "chromStart", "chromEnd", "rnaStart", "rnaEnd", )

    def __init__(self, trans, chromStart, chromEnd, rnaStart, rnaEnd):
        self.trans = trans
        self.chromStart, self.chromEnd = chromStart, chromEnd
        self.rnaStart, self.rnaEnd = rnaStart, rnaEnd

    @property
    def size(self):
        return self.chromEnd - self.chromStart


class ExonFeature(TransFeature):
    """exon with target gaps closed.."""
    __slots__ = ("rnaInsertCnt", "chromInsertCnt")

    def __init__(self, trans, chromStart, chromEnd, rnaStart, rnaEnd, rnaInsertCnt, chromInsertCnt):
        super(ExonFeature, self).__init__(trans, chromStart, chromEnd, rnaStart, rnaEnd)
        self.rnaInsertCnt, self.chromInsertCnt = rnaInsertCnt, chromInsertCnt

    def __str__(self):
        return "exon {}-{} rna={}-{} rIns={} cIns={}".format(self.chromStart, self.chromEnd,
                                                             self.rnaStart, self.rnaEnd,
                                                             self.rnaInsertCnt, self.chromInsertCnt)

    def reverseComplement(self, rcTrans):
        rcChromStart, rcChromEnd = dnaOps.reverseCoords(self.chromStart, self.chromEnd, self.trans.chromSize)
        rcRnaStart, rcRnaEnd = dnaOps.reverseCoords(self.rnaStart, self.rnaEnd, self.trans.rnaSize)
        return ExonFeature(rcTrans, rcChromStart, rcChromEnd, rcRnaStart, rcRnaEnd, self.rnaInsertCnt, self.chromInsertCnt)

    def rnaOverlaps(self, exon2):
        "does RNA rnage overlap?"
        return (self.rnaStart < exon2.rnaEnd) and (self.rnaEnd > exon2.rnaStart)


class IntronFeature(TransFeature):
    """intron from annotation or alignment, splice junction information maybe
    None"""
    __slots__ = ("rnaDeleteCnt", "donorSeq", "acceptorSeq", "spliceSites")

    def __init__(self, trans, chromStart, chromEnd, rnaStart, rnaEnd, rnaDeleteCnt, donorSeq, acceptorSeq, spliceSites):
        assert (rnaEnd - rnaStart) == rnaDeleteCnt
        super(IntronFeature, self).__init__(trans, chromStart, chromEnd, rnaStart, rnaEnd)
        self.rnaDeleteCnt, self.donorSeq, self.acceptorSeq, self.spliceSites = rnaDeleteCnt, donorSeq, acceptorSeq, spliceSites

    def __str__(self):
        if self.donorSeq is None:
            sjDesc = None
        else:
            sjDesc = "{}...{} ({})".format(self.donorSeq, self.acceptorSeq, self.spliceSites)
        return "intron {}-{} rna={}-{} rDel={} sjBases={}".format(self.chromStart, self.chromEnd,
                                                                  self.rnaStart, self.rnaEnd,
                                                                  self.rnaDeleteCnt, sjDesc)

    def reverseComplement(self, rcTrans):
        rcChromStart, rcChromEnd = dnaOps.reverseCoords(self.chromStart, self.chromEnd, self.trans.chromSize)
        rcRnaStart, rcRnaEnd = dnaOps.reverseCoords(self.rnaStart, self.rnaEnd, self.trans.rnaSize)
        rcDonorSeq, rcAcceptorSeq = (None, None) if self.donorSeq is None else (dnaOps.reverseComplement(self.acceptorSeq), dnaOps.reverseComplement(self.donorSeq))
        return IntronFeature(rcTrans, rcChromStart, rcChromEnd, rcRnaStart, rcRnaEnd,
                             self.rnaDeleteCnt, rcDonorSeq, rcAcceptorSeq, self.spliceSites)


class TranscriptFeatures(object):
    """
    Set of features for a transcript derived from an alignment or annotation,
    features are kept in chrom strand order.
    """

    def __init__(self, chrom, chromStrand, chromStart, chromEnd, chromSize,
                 rnaName, rnaStrand, rnaStart, rnaEnd, rnaSize,
                 cdsChromStart=None, cdsChromEnd=None):
        self.chrom, self.chromStrand, self.chromStart, self.chromEnd, self.chromSize = chrom, chromStrand, chromStart, chromEnd, chromSize
        self.rnaName, self.rnaStrand, self.rnaStart, self.rnaEnd, self.rnaSize = rnaName, rnaStrand, rnaStart, rnaEnd, rnaSize
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
        rcChromStart, rcChromEnd = dnaOps.reverseCoords(self.chromStart, self.chromEnd, self.chromSize)
        rcRnaStart, rcRnaEnd = dnaOps.reverseCoords(self.rnaStart, self.rnaEnd, self.rnaSize)
        rcCdsChromStart, rcCdsChromEnd = None, None if self.cdsChromStart is None else dnaOps.reverseCoords(self.cdsChromStart, self.cdsChromEnd, self.chromSize)

        rcTrans = TranscriptFeatures(self.chrom, dnaOps.reverseStrand(self.chromStrand), rcChromStart, rcChromEnd, self.chromSize,
                                     self.rnaName, dnaOps.reverseStrand(self.rnaStrand), rcRnaStart, rcRnaEnd, self.rnaSize,
                                     rcCdsChromStart, rcCdsChromEnd)
        features = []
        for i in xrange(len(self.features)-1, -1, -1):
            features.append(self.features[i].reverseComplement(rcTrans))
        rcTrans.features = tuple(features)
        return rcTrans
