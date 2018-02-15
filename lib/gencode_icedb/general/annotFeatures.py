from __future__ import print_function
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, Utr5RegionFeature, CdsRegionFeature, Utr3RegionFeature, NonCodingRegionFeature

# FIXME: consistent variable naming: annot, annotTrans transAnnot ..


class AnnotationGenePredFactory(object):
    """
    Factory to create annotation features from genePreds.
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

    def _buildFeatures(self, gp, trans):
        iBlkStart = 0
        rnaStart = 0
        features = []
        while iBlkStart < len(gp.exons):
            iBlkEnd, qCount = self._findExonEnd(gp, iBlkStart)
            rnaEnd = rnaStart + qCount
            features.append(self._makeExon(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, trans, len(features)))
            if iBlkEnd < len(gp.exons):
                features.append(self._makeIntron(gp, iBlkEnd, rnaEnd, trans, len(features)))
            rnaStart = rnaEnd
            iBlkStart = iBlkEnd
        return features

    def _tGapSize(self, gp, iBlk):
        "size of gap before the block"
        return gp.exons[iBlk].start - gp.exons[iBlk - 1].end

    def _findExonEnd(self, gp, iBlkStart):
        "finds half-open end of blocks covering current exon and total base, less closed gaps"
        iBlkEnd = iBlkStart + 1
        qCount = gp.exons[iBlkStart].size()
        while (iBlkEnd < len(gp.exons)) and (self._tGapSize(gp, iBlkEnd) < minIntronSize):
            qCount += gp.exons[iBlkEnd].size()
            iBlkEnd += 1
        return iBlkEnd, qCount

    def _getUtr5Annot(self, annot, rnaNext, exon, annotFeatures):
        utr5 = Utr5RegionFeature(exon, len(annotFeatures),
                                 exon.chrom.subrange(annot.start, annot.end),
                                 exon.rna.subrange(rnaNext, rnaNext + annot.size()))
        annotFeatures.append(utr5)
        return utr5.rna.end

    def _getCdsAnnot(self, annot, rnaNext, exon, frame, annotFeatures):
        cds = CdsRegionFeature(exon, len(annotFeatures),
                               exon.chrom.subrange(annot.start, annot.end),
                               exon.rna.subrange(rnaNext, rnaNext + annot.size()),
                               frame)
        annotFeatures.append(cds)
        return cds.rna.end

    def _getUtr3Annot(self, annot, rnaNext, exon, annotFeatures):
        utr3 = Utr3RegionFeature(exon, len(annotFeatures),
                                 exon.chrom.subrange(annot.start, annot.end),
                                 exon.rna.subrange(rnaNext, rnaNext + annot.size()))
        annotFeatures.append(utr3)
        return utr3.rna.end

    def _getCodingFeatures(self, blk, rnaNext, exon, annotFeatures):
        annot = blk.featureSplit()
        if blk.gene.strand == '+':
            if annot.utr5 is not None:
                rnaNext = self._getUtr5Annot(annot.utr5, rnaNext, exon, annotFeatures)
            if annot.cds is not None:
                rnaNext = self._getCdsAnnot(annot.cds, rnaNext, exon, Frame(blk.frame), annotFeatures)
            if annot.utr3 is not None:
                rnaNext = self._getUtr3Annot(annot.utr3, rnaNext, exon, annotFeatures)
        else:
            if annot.utr3 is not None:
                rnaNext = self._getUtr3Annot(annot.utr3, rnaNext, exon, annotFeatures)
            if annot.cds is not None:
                rnaNext = self._getCdsAnnot(annot.cds, rnaNext, exon, Frame(blk.frame) - blk.size(), annotFeatures)
            if annot.utr5 is not None:
                rnaNext = self._getUtr5Annot(annot.utr5, rnaNext, exon, annotFeatures)
        return rnaNext

    def _addCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        annotFeatures = []
        rnaNext = rnaStart
        for iBlk in range(iBlkStart, iBlkEnd):
            rnaNext = self._getCodingFeatures(gp.exons[iBlk], rnaNext, exon, annotFeatures)
        assert rnaNext == rnaEnd
        exon.annotFeatures = tuple(annotFeatures)

    def _addNonCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        nonCodingFeatures = []
        rnaNext = rnaStart
        for iBlk in range(iBlkStart, iBlkEnd):
            gpExon = gp.exons[iBlk]
            nonCodingFeatures.append(NonCodingRegionFeature(exon, len(nonCodingFeatures),
                                                            exon.chrom.subrange(gpExon.start, gpExon.end),
                                                            exon.rna.subrange(rnaNext, rnaNext + gpExon.size())))
            rnaNext += gpExon.size()
        assert rnaNext == rnaEnd, "rnaNext={}, rnaEnd={}".format(rnaNext, rnaEnd)
        exon.annotFeatures = tuple(nonCodingFeatures)

    def _makeExon(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, trans, iFeat):
        exon = ExonFeature(trans, iFeat,
                           trans.chrom.subrange(gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end),
                           trans.rna.subrange(rnaStart, rnaEnd))
        if trans.cdsChrom is not None:
            self._addCodingFeatures(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon)
        else:
            self._addNonCodingFeatures(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon)
        return exon

    def _getSpliceSites(self, gp, iBlkNext):
        if self.genomeReader is None:
            return (None, None)
        else:
            return spliceJuncsGetSeqs(self.genomeReader, gp.chrom,
                                      gp.exons[iBlkNext - 1].end,
                                      gp.exons[iBlkNext].start, gp.strand)

    def _makeIntron(self, gp, iBlkNext, rnaEnd, trans, iFeat):
        donorSeq, acceptorSeq = self._getSpliceSites(gp, iBlkNext)
        return IntronFeature(trans, iFeat,
                             trans.chrom.subrange(gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start),
                             trans.rna.subrange(rnaEnd, rnaEnd), donorSeq, acceptorSeq)

    def fromGenePred(self, gp, attrs=None):
        "convert a genePred to an AnnotTranscript"
        chromSize = self.chromSizeFunc(gp.chrom)
        rnaSize = gp.getLenExons()
        if gp.cdsStart < gp.cdsEnd:
            cdsChrom = Coords(gp.chrom, gp.cdsStart, gp.cdsEnd, strand='+', size=chromSize)
        else:
            cdsChrom = None

        chrom = Coords(gp.chrom, gp.txStart, gp.txEnd, '+', chromSize)
        rna = Coords(gp.name, 0, rnaSize, gp.strand, rnaSize)
        trans = TranscriptFeatures(chrom, rna, cdsChrom, attrs)
        trans.features = tuple(self._buildFeatures(gp, trans))
        return trans
