from __future__ import print_function
from pycbio.hgdata.hgLite import GenePredDbTable
from gencode_icedb.genome import spliceSitesClassifyStrand
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.tsl.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures


class AnnotationGenePredFactory(object):
    """
    factory to create annotation features from genePreds
    """

    def __init__(self, genomeReader=None, chromSizeFunc=None):
        """genomeReader is used to obtain splice sites and chrom size, maybe
        None if splice sites are reverse complement will not be done. To get chrom sizes
        without splice sites, provide the chromSizeFunc(chrom) function."""
        self.genomeReader = genomeReader
        self.chromSizeFunc = chromSizeFunc

    def __buildFeatures(self, gp):
        iBlkStart = 0
        qStart = 0
        while iBlkStart < len(gp.exons):
            iBlkEnd, qCount = self.__findExonEnd(gp, iBlkStart)
            qEnd = qStart + qCount
            yield self.__makeExon(gp, iBlkStart, iBlkEnd, qStart, qEnd)
            if iBlkEnd < len(gp.exons):
                yield self.__makeIntron(gp, iBlkEnd)
            qStart = qEnd
            iBlkStart = iBlkEnd

    def __tGapSize(self, gp, iBlk):
        "size of gap before the block"
        return gp.exons[iBlk].start - gp.exons[iBlk - 1].end

    def __findExonEnd(self, gp, iBlkStart):
        "finds half-open end of blocks covering current exon and total base, less closed gaps"
        iBlkEnd = iBlkStart + 1
        qCount = gp.exons[iBlkStart].size()
        while (iBlkEnd < len(gp.exons)) and (self.__tGapSize(gp, iBlkEnd) < minIntronSize):
            qCount = gp.exons[iBlkEnd].size()
            iBlkEnd += 1
        return iBlkEnd, qCount

    def __countGap(self, gp, iBlkStart, iBlkEnd):
        gapBases = 0
        for iBlk in xrange(iBlkStart + 1, iBlkEnd):
            gapBases += gp.exons[iBlk].start - gp.exons[iBlk - 1].end
        return gapBases

    def __makeExon(self, gp, iBlkStart, iBlkEnd, qStart, qEnd):
        return ExonFeature(self, gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end,
                           qStart, qEnd, 0, self.__countGap(gp, iBlkStart, iBlkEnd))

    def __getSpliceSites(self, gp, iBlkNext):
        donorSeq = self.genomeReader.get(gp.chrom, gp.exons[iBlkNext - 1].end,
                                         gp.exons[iBlkNext - 1].end + 2)
        acceptorSeq = self.genomeReader.get(gp.chrom, gp.exons[iBlkNext].start - 2,
                                            gp.exons[iBlkNext].start)
        spliceSites = spliceSitesClassifyStrand(gp.strand, donorSeq, acceptorSeq)
        return donorSeq, acceptorSeq, spliceSites

    def __makeIntron(self, gp, iBlkNext):
        if self.genomeReader is None:
            donorSeq = acceptorSeq = spliceSites = None
        else:
            donorSeq, acceptorSeq, spliceSites = self.__getSpliceSites(gp, iBlkNext)
        return IntronFeature(self, gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start,
                             0, donorSeq, acceptorSeq, spliceSites)

    def __getOptionalChromSize(self, chrom):
        if self.genomeReader is not None:
            return self.genomeReader.getChromSize(chrom)
        else:
            return None
    
    def fromGenePred(self, gp):
        "convert a genePred to an AnnotTranscript"
        rnaSize = gp.getLenExons()
        if gp.cdsStart < gp.cdsEnd:
            cdsChromStart, cdsChromEnd = gp.cdsStart, gp.cdsEnd
        else:
            cdsChromStart = cdsChromEnd = None

        return TranscriptFeatures(gp.chrom, '+', gp.txStart, gp.txEnd, self.__getOptionalChromSize(gp.chrom),
                                  gp.name, gp.strand, 0, rnaSize, rnaSize, self.__buildFeatures(gp),
                                  cdsChromStart, cdsChromEnd)


class AnnotationFeatures(list):
    "table of AnnotTranscript objects"

    def __init__(self):
        self.transcriptByName = {}

    def addTranscript(self, annotTrans):
        self.transcriptByName[annotTrans.name].append(annotTrans)
        for annotFeat in annotTrans:
            if isinstance(annotFeat, ExonFeature):
                self.exonRangeMap.add(annotTrans.chrom, annotFeat.start, annotFeat.end, annotTrans, annotTrans.strand)

    @staticmethod
    def dbFactory(conn, table, chrom, chromStart, chromEnd, genomeReader):
        "constructor from a sqlite3 databases"
        gpDbTable = GenePredDbTable(conn, table)
        annotFactory = AnnotationGenePredFactory(genomeReader)
        annotFeatureMap = AnnotationFeatures()
        for gp in gpDbTable.getRangeOverlap(chrom, chromStart, chromEnd):
            annotFeatureMap.addTranscript(annotFactory.fromGenePred(gp))
        return annotFeatureMap
