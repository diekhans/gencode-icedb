from __future__ import print_function
from pycbio.hgdata.hgLite import GenePredDbTable
from pycbio.hgdata.frame import Frame
from gencode_icedb.genome import spliceSitesClassifyStrand
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.tsl.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, Utr5RegionFeature, CdsRegionFeature, Utr3RegionFeature


class AnnotationGenePredFactory(object):
    """
    factory to create annotation features from genePreds
    """

    def __init__(self, genomeReader=None, chromSizeFunc=None):
        """genomeReader is used to obtain splice sites and chrom size, maybe
        None if splice sites are reverse complement will not be done. To get chrom sizes
        without splice sites, provide the chromSizeFunc(chrom) function."""
        self.genomeReader = genomeReader
        # set chromSizeFunc to return size or None if not available
        if chromSizeFunc is not None:
            self.chromSizeFunc = chromSizeFunc
        elif self.genomeReader is not None:
            self.chromSizeFunc = genomeReader.getChromSize
        else:
            self.chromSizeFunc = lambda chrom: None

    def __buildFeatures(self, gp, trans):
        iBlkStart = 0
        rnaStart = 0
        while iBlkStart < len(gp.exons):
            iBlkEnd, qCount = self.__findExonEnd(gp, iBlkStart)
            rnaEnd = rnaStart + qCount
            yield self.__makeExon(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, trans)
            if iBlkEnd < len(gp.exons):
                yield self.__makeIntron(gp, iBlkEnd, rnaEnd, trans)
            rnaStart = rnaEnd
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

    def __getUtr5Annot(self, annot, rnaNext, exon, codingFeatures):
        utr5 = Utr5RegionFeature(exon, annot.start, annot.end, rnaNext, rnaNext + annot.size())
        codingFeatures.append(utr5)
        return utr5.rnaEnd
    
    def __getCdsAnnot(self, annot, rnaNext, exon, frame, codingFeatures):
        cds = CdsRegionFeature(exon, annot.start, annot.end, rnaNext, rnaNext + annot.size(), frame)
        codingFeatures.append(cds)
        return cds.rnaEnd
    
    def __getUtr3Annot(self, annot, rnaNext, exon, codingFeatures):
        utr3 = Utr5RegionFeature(exon, annot.start, annot.end, rnaNext, rnaNext + annot.size())
        codingFeatures.append(utr3)
        return utr3.rnaEnd
    
    def __getCodingFeatures(self, blk, rnaNext, exon, codingFeatures):
        annot = blk.featureSplit()
        if blk.gene.strand == '+':
            if annot.utr5 is not None:
                rnaNext = self.__getUtr5Annot(annot.utr5, rnaNext, exon, codingFeatures)
            if annot.cds is not None:
                rnaNext = self.__getCdsAnnot(annot.cds, rnaNext, exon, Frame(blk.frame), codingFeatures)
            if annot.utr3 is not None:
                rnaNext = self.__getUtr3Annot(annot.utr3, rnaNext, exon, codingFeatures)
        else:
            if annot.utr3 is not None:
                rnaNext = self.__getUtr3Annot(annot.utr3, rnaNext, exon, codingFeatures)
            if annot.cds is not None:
                rnaNext = self.__getCdsAnnot(annot.cds, rnaNext, exon, Frame(blk.frame) - blk.size(), codingFeatures)
            if annot.utr5 is not None:
                rnaNext = self.__getUtr5Annot(annot.utr5, rnaNext, exon, codingFeatures)
        return rnaNext
        
    def __addCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        codingFeatures = []
        rnaNext = rnaStart
        for iBlk in xrange(iBlkStart, iBlkEnd):
            rnaNext = self.__getCodingFeatures(gp.exons[iBlk], rnaNext, exon, codingFeatures)
        assert rnaNext == rnaEnd
        exon.codingFeatures = tuple(codingFeatures)
    
    def __makeExon(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, trans):
        exon = ExonFeature(trans, gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end,
                           rnaStart, rnaEnd)
        if trans.cdsChromStart < trans.cdsChromEnd:
            self.__addCodingFeatures(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon)
        return exon

    def __getSpliceSites(self, gp, iBlkNext):
        donorSeq = self.genomeReader.get(gp.chrom, gp.exons[iBlkNext - 1].end,
                                         gp.exons[iBlkNext - 1].end + 2)
        acceptorSeq = self.genomeReader.get(gp.chrom, gp.exons[iBlkNext].start - 2,
                                            gp.exons[iBlkNext].start)
        spliceSites = spliceSitesClassifyStrand(gp.strand, donorSeq, acceptorSeq)
        return donorSeq, acceptorSeq, spliceSites

    def __makeIntron(self, gp, iBlkNext, rnaEnd, trans):
        if self.genomeReader is None:
            donorSeq = acceptorSeq = spliceSites = None
        else:
            donorSeq, acceptorSeq, spliceSites = self.__getSpliceSites(gp, iBlkNext)
        return IntronFeature(trans, gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start,
                             rnaEnd, rnaEnd, donorSeq, acceptorSeq, spliceSites)

    def fromGenePred(self, gp):
        "convert a genePred to an AnnotTranscript"
        rnaSize = gp.getLenExons()
        if gp.cdsStart < gp.cdsEnd:
            cdsChromStart, cdsChromEnd = gp.cdsStart, gp.cdsEnd
        else:
            cdsChromStart = cdsChromEnd = None

        trans = TranscriptFeatures(gp.chrom, '+', gp.txStart, gp.txEnd, self.chromSizeFunc(gp.chrom),
                                   gp.name, gp.strand, 0, rnaSize, rnaSize,
                                   cdsChromStart, cdsChromEnd)
        trans.features = tuple(self.__buildFeatures(gp, trans))
        return trans


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
