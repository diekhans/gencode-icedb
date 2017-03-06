from __future__ import print_function
from pycbio.hgdata.rangeFinder import RangeFinder
from pycbio.hgdata.hgLite import PslDbTable
from gencode_icedb.sequence import spliceSitesClassifyStrand
from gencode_icedb import minIntronSize


class AnnotFeature(object):
    "an annotation feature on the genome"
    __slots__ = ("parent", "start", "end")

    def __init__(self, parent, start, end):
        self.parent, self.start, self.end = parent, start, end


class AnnotExon(AnnotFeature):
    "exon from annotation, with target deletions closed and counted"
    __slots__ = ("tGapBases",)

    def __init__(self, parent, start, end, tGapBases):
        super(AnnotExon, self).__init__(parent, start, end)
        self.tGapBases = tGapBases

    def __str__(self):
        return "exon {}-{} tGap={}".format(self.start, self.end, self.tGapBases)


class AnnotIntron(AnnotFeature):
    "intron from annotation"
    __slots__ = ("startBases", "endBases", "spliceSites")

    def __init__(self, parent, start, end, startBases, endBases, spliceSites):
        super(AnnotIntron, self).__init__(parent, start, end)
        self.startBases, self.endBases, self.spliceSites = startBases, endBases, spliceSites

    def __str__(self):
        return "intron {}-{} sjBases={}...{} ({})".format(self.start, self.end,
                                                          self.startBases, self.endBases, self.spliceSites)


class AnnotTranscript(list):
    """
    Set of features for a transcript annotation. Features are
    kept in genomic order.
    """

    def __init__(self, name, chrom, strand, start, end, features):
        self.name, self.chrom, self.strand, self.start, self.end = name, chrom, strand, start, end,
        self.extend(features)

    def __str__(self):
        return "annot={} {}:{}-{} {}".format(self.name, self.chrom, self.start, self.end, self.strand)


class AnnotTranscriptGenePredFactory(object):
    """
    factory to create annotation features from genePreds
    """

    def __init__(self, genomeReader):
        self.genomeReader = genomeReader

    def __buildFeatures(self, gp):
        iBlkStart = 0
        while iBlkStart < len(gp.exons):
            iBlkEnd = self.__findExonEnd(gp, iBlkStart)
            yield self.__makeExon(gp, iBlkStart, iBlkEnd)
            if iBlkEnd < len(gp.exons):
                yield self.__makeIntron(gp, iBlkEnd)
            iBlkStart = iBlkEnd

    def __tGapSize(self, gp, iBlk):
        "size of gap before the block"
        return gp.exons[iBlk].start - gp.exons[iBlk - 1].end

    def __findExonEnd(self, gp, iBlkStart):
        "finds half-open end of blocks covering current exon"
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(gp.exons)) and (self.__tGapSize(gp, iBlkEnd) < minIntronSize):
            iBlkEnd += 1
        return iBlkEnd

    def __makeExon(self, gp, iBlkStart, iBlkEnd):
        gapBases = 0
        for iBlk in xrange(iBlkStart + 1, iBlkEnd):
            gapBases += gp.exons[iBlk].start - gp.exons[iBlk - 1].end
        return AnnotExon(self, gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end, gapBases)

    def __makeIntron(self, gp, iBlkNext):
        startBases = self.genomeReader.get(gp.tName, gp.exons[iBlkNext - 1].end,
                                           gp.exons[iBlkNext - 1].end + 2)
        endBases = self.genomeReader.get(gp.tName, gp.exons[iBlkNext].start - 2,
                                         gp.exons[iBlkNext].start)
        spliceSites = spliceSitesClassifyStrand(gp.getQStrand(), startBases, endBases)
        return AnnotIntron(self, gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start,
                          startBases, endBases, spliceSites)

    def fromGenePred(self, gp):
        "convert a genePred to an AnnotTranscript"
        return AnnotTranscript(gp.name, gp.chom, gp.strand, gp.txStart, gp.txEnd,
                               self.__buildFeatures(gp))


class AnnotFeatureMap(list):
    "map by exon coordinates of AnnotFeatures objects"

    def __init__(self):
        self.transcripts = []
        self.exonRangeMap = RangeFinder()

    def overlapping(self, chrom, start, end, strand=None):
        "generator over ExonFeatures overlaping the specified range"
        return self.exonRangeMap.overlapping(chrom, start, end, strand)

    def addTranscript(self, annotTrans):
        self.transcripts.append(annotTrans)
        for annotFeat in annotTrans:
            if isinstance(annotFeat, AnnotExon):
                self.exonRangeMap.add(annotTrans.chrom, annotFeat.start, annotFeat.end, annotTrans, annotTrans.strand)

    @staticmethod
    def dbFactory(conn, table, chrom, start, end, genomeReader):
        "constructor from a sqlite3 databases"
        pslDbTable = PslDbTable(conn, table)
        annotFactory = AnnotTranscriptPslFactory(genomeReader)
        annotFeatureMap = AnnotFeatureMap()
        for psl in pslDbTable.getTRangeOverlap(chrom, start, end):
            annotFeatureMap.addTranscript(annotFactory.fromPsl(psl))
        return annotFeatureMap
