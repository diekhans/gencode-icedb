from __future__ import print_function
from pycbio.hgdata.hgLite import GenePredDbTable
from gencode_icedb.genome import spliceSitesClassifyStrand
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.tsl.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures


class AnnotationGenePredFactory(object):
    """
    factory to create annotation features from genePreds
    """

    def __init__(self, genomeReader):
        """genomeReader maybe None if splice sites are not desired """
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
        return ExonFeature(self, gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end, 0, gapBases)

    def __getSpliceSites(self, gp, iBlkNext):
        startBases = self.genomeReader.get(gp.chrom, gp.exons[iBlkNext - 1].end,
                                           gp.exons[iBlkNext - 1].end + 2)
        endBases = self.genomeReader.get(gp.chrom, gp.exons[iBlkNext].start - 2,
                                         gp.exons[iBlkNext].start)
        spliceSites = spliceSitesClassifyStrand(gp.strand, startBases, endBases)
        return startBases, endBases, spliceSites

    def __makeIntron(self, gp, iBlkNext):
        if self.genomeReader is None:
            startBases = endBases = spliceSites = None
        else:
            startBases, endBases, spliceSites = self.__getSpliceSites(gp, iBlkNext)
        return IntronFeature(self, gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start,
                             0, startBases, endBases, spliceSites)

    def fromGenePred(self, gp):
        "convert a genePred to an AnnotTranscript"
        qSize = gp.getLenExons()
        return TranscriptFeatures(gp.chrom, gp.strand, gp.txStart, gp.txEnd, None,
                                  gp.name, 0, qSize, qSize,
                                  self.__buildFeatures(gp))


class AnnotationFeatures(list):
    "table of AnnotTranscript objects"

    def __init__(self):
        self.transcriptByName = []

    def addTranscript(self, annotTrans):
        self.transcriptByName[annotTrans.name].append(annotTrans)
        for annotFeat in annotTrans:
            if isinstance(annotFeat, ExonFeature):
                self.exonRangeMap.add(annotTrans.chrom, annotFeat.start, annotFeat.end, annotTrans, annotTrans.strand)

    @staticmethod
    def dbFactory(conn, table, chrom, start, end, genomeReader):
        "constructor from a sqlite3 databases"
        gpDbTable = GenePredDbTable(conn, table)
        annotFactory = AnnotationGenePredFactory(genomeReader)
        annotFeatureMap = AnnotationFeatures()
        for gp in gpDbTable.getRangeOverlap(chrom, start, end):
            annotFeatureMap.addTranscript(annotFactory.fromGenePred(gp))
        return annotFeatureMap
