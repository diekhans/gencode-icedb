from __future__ import print_function
from collections import defaultdict
from pycbio.hgdata.hgLite import GenePredDbTable
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from pycbio.hgdata.rangeFinder import RangeFinder
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, Utr5RegionFeature, CdsRegionFeature, Utr3RegionFeature, NonCodingRegionFeature


class AnnotationGenePredFactory(object):
    """
    factory to create annotation features from genePreds
    """

    def __init__(self, genomeReader=None, chromSizeFunc=None):
        """genomeReader is used to obtain splice sites and chrom size, maybe
        None if splice sites and reverse complement will not be done. To get chrom sizes
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
            qCount += gp.exons[iBlkEnd].size()
            iBlkEnd += 1
        return iBlkEnd, qCount

    def __getUtr5Annot(self, annot, rnaNext, exon, rnaFeatures):
        utr5 = Utr5RegionFeature(exon,
                                 exon.chrom.subrange(annot.start, annot.end),
                                 exon.rna.subrange(rnaNext, rnaNext + annot.size()))
        rnaFeatures.append(utr5)
        return utr5.rna.end

    def __getCdsAnnot(self, annot, rnaNext, exon, frame, rnaFeatures):
        cds = CdsRegionFeature(exon,
                               exon.chrom.subrange(annot.start, annot.end),
                               exon.rna.subrange(rnaNext, rnaNext + annot.size()),
                               frame)
        rnaFeatures.append(cds)
        return cds.rna.end

    def __getUtr3Annot(self, annot, rnaNext, exon, rnaFeatures):
        utr3 = Utr3RegionFeature(exon,
                                 exon.chrom.subrange(annot.start, annot.end),
                                 exon.rna.subrange(rnaNext, rnaNext + annot.size()))
        rnaFeatures.append(utr3)
        return utr3.rna.end

    def __getCodingFeatures(self, blk, rnaNext, exon, rnaFeatures):
        annot = blk.featureSplit()
        if blk.gene.strand == '+':
            if annot.utr5 is not None:
                rnaNext = self.__getUtr5Annot(annot.utr5, rnaNext, exon, rnaFeatures)
            if annot.cds is not None:
                rnaNext = self.__getCdsAnnot(annot.cds, rnaNext, exon, Frame(blk.frame), rnaFeatures)
            if annot.utr3 is not None:
                rnaNext = self.__getUtr3Annot(annot.utr3, rnaNext, exon, rnaFeatures)
        else:
            if annot.utr3 is not None:
                rnaNext = self.__getUtr3Annot(annot.utr3, rnaNext, exon, rnaFeatures)
            if annot.cds is not None:
                rnaNext = self.__getCdsAnnot(annot.cds, rnaNext, exon, Frame(blk.frame) - blk.size(), rnaFeatures)
            if annot.utr5 is not None:
                rnaNext = self.__getUtr5Annot(annot.utr5, rnaNext, exon, rnaFeatures)
        return rnaNext

    def __addCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        rnaFeatures = []
        rnaNext = rnaStart
        for iBlk in range(iBlkStart, iBlkEnd):
            rnaNext = self.__getCodingFeatures(gp.exons[iBlk], rnaNext, exon, rnaFeatures)
        assert rnaNext == rnaEnd
        exon.rnaFeatures = tuple(rnaFeatures)

    def __addNonCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        nonCodingFeatures = []
        rnaNext = rnaStart
        for iBlk in range(iBlkStart, iBlkEnd):
            gpExon = gp.exons[iBlk]
            nonCodingFeatures.append(NonCodingRegionFeature(exon,
                                                            exon.chrom.subrange(gpExon.start, gpExon.end),
                                                            exon.rna.subrange(rnaNext, rnaNext + gpExon.size())))
            rnaNext += gpExon.size()
        assert rnaNext == rnaEnd, "rnaNext={}, rnaEnd={}".format(rnaNext, rnaEnd)
        exon.rnaFeatures = tuple(nonCodingFeatures)

    def __makeExon(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, trans):
        exon = ExonFeature(trans,
                           trans.chrom.subrange(gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end),
                           trans.rna.subrange(rnaStart, rnaEnd))
        if trans.cdsChromStart is not None:
            self.__addCodingFeatures(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon)
        else:
            self.__addNonCodingFeatures(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon)
        return exon

    def __getSpliceSites(self, gp, iBlkNext):
        if self.genomeReader is None:
            return (None, None)
        else:
            return spliceJuncsGetSeqs(self.genomeReader, gp.chrom,
                                      gp.exons[iBlkNext - 1].end,
                                      gp.exons[iBlkNext].start, gp.strand)

    def __makeIntron(self, gp, iBlkNext, rnaEnd, trans):
        donorSeq, acceptorSeq = self.__getSpliceSites(gp, iBlkNext)
        return IntronFeature(trans,
                             trans.chrom.subrange(gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start),
                             trans.rna.subrange(rnaEnd, rnaEnd), donorSeq, acceptorSeq)

    def fromGenePred(self, gp):
        "convert a genePred to an AnnotTranscript"
        rnaSize = gp.getLenExons()
        if gp.cdsStart < gp.cdsEnd:
            cdsChromStart, cdsChromEnd = gp.cdsStart, gp.cdsEnd
        else:
            cdsChromStart = cdsChromEnd = None

        chrom = Coords(gp.chrom, gp.txStart, gp.txEnd, '+', self.chromSizeFunc(gp.chrom))
        rna = Coords(gp.name, 0, rnaSize, gp.strand, rnaSize)
        trans = TranscriptFeatures(chrom, rna, cdsChromStart, cdsChromEnd)
        trans.features = tuple(self.__buildFeatures(gp, trans))
        return trans


class AnnotationFeatures(list):
    "table of AnnotTranscript objects"

    def __init__(self):
        self.transcriptsByName = defaultdict(list)
        self.transcriptsByRange = RangeFinder()

    def addTranscript(self, annotTrans):
        self.transcriptsByName[annotTrans.name].append(annotTrans)
        self.transcriptsByRange.add(annotTrans.chrom.name, annotTrans.chrom.start, annotTrans.chrom.end, annotTrans, annotTrans.rna.strand)
        self.append(annotTrans)

    @staticmethod
    def dbFactory(conn, table, chrom, chromStart, chromEnd, genomeReader):
        "constructor from a sqlite3 databases"
        gpDbTable = GenePredDbTable(conn, table)
        annotFactory = AnnotationGenePredFactory(genomeReader)
        annotFeatureMap = AnnotationFeatures()
        for gp in gpDbTable.getRangeOverlap(chrom, chromStart, chromEnd):
            annotFeatureMap.addTranscript(annotFactory.fromGenePred(gp))
        return annotFeatureMap
