"""
create AnnotationFeature objects from UCSC genePred records.
"""
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.genePred import GenePred
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs
from gencode_icedb.tsl import minIntronSize  # FIXME: this should be in general
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, Utr5RegionFeature, CdsRegionFeature, Utr3RegionFeature, GapAnnotFeature, NonCodingRegionFeature

# FIXME: consistent variable naming: annot, annotTrans transAnnot ..


class GenePredAnnotationFactory(object):
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

    def _buildFeatures(self, gp, transAnnot):
        iBlkStart = 0
        rnaStart = 0
        features = []
        while iBlkStart < len(gp.exons):
            iBlkEnd, qCount = self._findExonEnd(gp, iBlkStart)
            rnaEnd = rnaStart + qCount
            features.append(self._makeExon(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, transAnnot, len(features)))
            if iBlkEnd < len(gp.exons):
                features.append(self._makeIntron(gp, iBlkEnd, rnaEnd, transAnnot, len(features)))
            rnaStart = rnaEnd
            iBlkStart = iBlkEnd
        return features

    def _tGapSize(self, gp, iBlk):
        "size of gap before the block"
        return gp.exons[iBlk].start - gp.exons[iBlk - 1].end

    def _findExonEnd(self, gp, iBlkStart):
        """finds half-open end of blocks covering current exon and total base,
        including closed gaps"""
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(gp.exons)) and (self._tGapSize(gp, iBlkEnd) < minIntronSize):
            iBlkEnd += 1
        return iBlkEnd, gp.exons[iBlkEnd - 1].end - gp.exons[iBlkStart].start

    def _findRnaSize(self, gp):
        iBlkStart = 0
        rnaSize = 0
        while iBlkStart < len(gp.exons):
            iBlkEnd, qCount = self._findExonEnd(gp, iBlkStart)
            rnaSize += qCount
            iBlkStart = iBlkEnd
        return rnaSize

    def _addGapAnnot(self, chromStart, chromEnd, rnaNext, exon, features):
        gap = GapAnnotFeature(exon, len(features),
                              exon.chrom.subrange(chromStart, chromEnd),
                              exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)))
        features.append(gap)
        return gap.rna.end

    def _addUtr5Annot(self, chromStart, chromEnd, rnaNext, exon, features):
        utr5 = Utr5RegionFeature(exon, len(features),
                                 exon.chrom.subrange(chromStart, chromEnd),
                                 exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)))
        features.append(utr5)
        return utr5.rna.end

    def _addCdsAnnot(self, chromStart, chromEnd, rnaNext, exon, frame, features):
        cds = CdsRegionFeature(exon, len(features),
                               exon.chrom.subrange(chromStart, chromEnd),
                               exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)),
                               frame)
        features.append(cds)
        return cds.rna.end

    def _addUtr3Annot(self, chromStart, chromEnd, rnaNext, exon, features):
        utr3 = Utr3RegionFeature(exon, len(features),
                                 exon.chrom.subrange(chromStart, chromEnd),
                                 exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)))
        features.append(utr3)
        return utr3.rna.end

    def _addBlkCodingFeatures(self, blk, rnaNext, exon, features):
        # note that blk may not cover whole exon if there are annotation gaps
        cds = exon.transcript.cdsChrom
        chromNext = blk.start
        if chromNext < cds.start:
            chromEnd = min(blk.end, cds.start)
            if blk.gene.strand == '+':
                rnaNext = self._addUtr5Annot(chromNext, chromEnd, rnaNext, exon, features)
            else:
                rnaNext = self._addUtr3Annot(chromNext, chromEnd, rnaNext, exon, features)
            chromNext = chromEnd
        if (chromNext < blk.end) and (chromNext < cds.end):
            chromEnd = min(blk.end, cds.end)
            rnaNext = self._addCdsAnnot(chromNext, chromEnd, rnaNext, exon, Frame(blk.frame), features)
            chromNext = chromEnd
        if (chromNext < blk.end) and (chromNext >= cds.end):
            if blk.gene.strand == '+':
                rnaNext = self._addUtr3Annot(chromNext, blk.end, rnaNext, exon, features)
            else:
                rnaNext = self._addUtr5Annot(chromNext, blk.end, rnaNext, exon, features)
        return rnaNext

    def _addCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        features = []
        rnaNext = rnaStart
        chromNext = gp.exons[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            if chromNext != gp.exons[iBlk].start:
                rnaNext = self._addGapAnnot(chromNext, gp.exons[iBlk].start, rnaNext, exon, features)
            rnaNext = self._addBlkCodingFeatures(gp.exons[iBlk], rnaNext, exon, features)
            chromNext = features[-1].chrom.end
        exon.annotFeatures = tuple(features)
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _addNonCodingFeature(self, blk, rnaNext, exon, features):
        feat = NonCodingRegionFeature(exon, len(features),
                                      exon.chrom.subrange(blk.start, blk.end),
                                      exon.rna.subrange(rnaNext, rnaNext + blk.size()))
        features.append(feat)
        return feat.rna.end

    def _addNonCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        features = []
        rnaNext = rnaStart
        chromNext = gp.exons[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            if chromNext != gp.exons[iBlk].start:
                rnaNext = self._addGapAnnot(chromNext, gp.exons[iBlk].start, rnaNext, exon, features)
            rnaNext = self._addNonCodingFeature(gp.exons[iBlk], rnaNext, exon, features)
            chromNext = features[-1].chrom.end
        exon.annotFeatures = tuple(features)
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _makeExon(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, transAnnot, iFeat):
        exon = ExonFeature(transAnnot, iFeat,
                           transAnnot.chrom.subrange(gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end),
                           transAnnot.rna.subrange(rnaStart, rnaEnd))
        if transAnnot.cdsChrom is not None:
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

    def _makeIntron(self, gp, iBlkNext, rnaEnd, transAnnot, iFeat):
        donorSeq, acceptorSeq = self._getSpliceSites(gp, iBlkNext)
        return IntronFeature(transAnnot, iFeat,
                             transAnnot.chrom.subrange(gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start),
                             transAnnot.rna.subrange(rnaEnd, rnaEnd), donorSeq, acceptorSeq)

    def fromGenePred(self, gp: GenePred, attrs: ObjDict=None):
        "convert a genePred to an TranscriptFeatures"
        chromSize = self.chromSizeFunc(gp.chrom)
        rnaSize = self._findRnaSize(gp)
        if gp.cdsStart < gp.cdsEnd:
            cdsChrom = Coords(gp.chrom, gp.cdsStart, gp.cdsEnd, strand='+', size=chromSize)
        else:
            cdsChrom = None

        chrom = Coords(gp.chrom, gp.txStart, gp.txEnd, '+', chromSize)
        rna = Coords(gp.name, 0, rnaSize, gp.strand, rnaSize)
        transAnnot = TranscriptFeatures(chrom, rna, transcriptionStrand=gp.strand, cdsChrom=cdsChrom, attrs=attrs)
        transAnnot.features = tuple(self._buildFeatures(gp, transAnnot))
        assert sum([len(e.rna) for e in transAnnot.getFeaturesOfType(ExonFeature)]) == rna.size
        return transAnnot
