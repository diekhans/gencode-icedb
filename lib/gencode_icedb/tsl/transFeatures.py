"""
Features of a transcript annotation or alignment.
"""


class TransFeature(object):
    "a feature of annotation or alignment to the genome"
    __slots__ = ("parent", "start", "end")

    def __init__(self, parent, start, end):
        self.parent, self.start, self.end = parent, start, end


class ExonFeature(TransFeature):
    "exon with target gaps closed"
    __slots__ = ("qInsertBases", "tInsertBases")

    def __init__(self, parent, start, end, qInsertBases, tInsertBases):
        super(ExonFeature, self).__init__(parent, start, end)
        self.qInsertBases, self.tInsertBases = qInsertBases, tInsertBases

    def __str__(self):
        return "exon {}-{} qIns={} tIns={}".format(self.start, self.end,
                                                   self.qInsertBases, self.tInsertBases)


class IntronFeature(TransFeature):
    """intron from annotation or alignment, splice junction information maybe
    None"""
    __slots__ = ("qDeleteBases", "startBases", "endBases", "spliceSites")

    def __init__(self, parent, start, end, qDeleteBases, startBases, endBases, spliceSites):
        super(IntronFeature, self).__init__(parent, start, end)
        self.qDeleteBases, self.startBases, self.endBases, self.spliceSites = qDeleteBases, startBases, endBases, spliceSites

    def __str__(self):
        if self.startBases is None:
            sjDesc = None
        else:
            sjDesc = "{}...{} ({})".format(self.startBases, self.endBases, self.spliceSites)
        return "intron {}-{} qDel={} sjBases={}".format(self.start, self.end, self.qDeleteBases, sjDesc)


class TranscriptFeatures(list):
    """
    Set of features for a transcript derived from an alignment or annotation,
    features are kept in genomic order.
    """

    def __init__(self, chrom, strand, start, end, size, qName, qStart, qEnd, qSize,
                 features):
        self.chrom, self.strand, self.start, self.end, self.size = chrom, strand, start, end, size
        self.qName, self.qStart, self.qEnd, self.qSize = qName, qStart, qEnd, qSize
        self.extend(features)

    def __str__(self):
        return "t={}:{}-{} {}, q={}:{}={} {}".format(self.chrom, self.start, self.end, self.strand,
                                                     self.qName, self.qStart, self.qEnd, self.qSize)
