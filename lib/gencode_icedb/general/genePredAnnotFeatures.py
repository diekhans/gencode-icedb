"""
create AnnotationFeature objects from UCSC genePred records.
"""
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.genePred import GenePred
from gencode_icedb.tsl import minIntronSize  # FIXME: this should be in general
from gencode_icedb.general.annotFeatureBuilder import AnnotFeatureBuilder

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

    def _buildFeatures(self, gp, builder):
        iBlkStart = 0
        rnaStart = 0
        while iBlkStart < len(gp.exons):
            iBlkEnd, qCount = self._findExonEnd(gp, iBlkStart)
            rnaEnd = rnaStart + qCount
            self._makeExon(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
            if iBlkEnd < len(gp.exons):
                self._makeIntron(gp, iBlkEnd, rnaEnd, builder)
            rnaStart = rnaEnd
            iBlkStart = iBlkEnd

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

    def _addCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        rnaNext = rnaStart
        chromNext = gp.exons[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            blk = gp.exons[iBlk]
            if chromNext != blk.start:
                rnaNext = builder.addGap(chromNext, blk.start, rnaNext)
                chromNext = blk.start
            frame = Frame(blk.frame) if blk.frame >= 0 else None
            rnaNext = builder.addCodingFeatures(chromNext, blk.end, rnaNext, frame)
            chromNext = blk.end
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _addNonCodingFeatures(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        rnaNext = rnaStart
        chromNext = gp.exons[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            blk = gp.exons[iBlk]
            if chromNext != blk.start:
                rnaNext = builder.addGap(chromNext, blk.start, rnaNext)
                chromNext = blk.start
            rnaNext = builder.addNonCodingFeature(chromNext, blk.end, rnaNext)
            chromNext = blk.end
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _makeExon(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        builder.beginExon(gp.exons[iBlkStart].start, gp.exons[iBlkEnd - 1].end,
                          rnaStart, rnaEnd)
        if builder.hasCds:
            self._addCodingFeatures(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
        else:
            self._addNonCodingFeatures(gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
        builder.finishExon()

    def _makeIntron(self, gp, iBlkNext, rnaEnd, builder):
        builder.addIntron(gp.exons[iBlkNext - 1].end, gp.exons[iBlkNext].start, rnaEnd)

    def fromGenePred(self, gp: GenePred, attrs: ObjDict = None):
        "convert a genePred to an TranscriptFeatures"
        chromSize = self.chromSizeFunc(gp.chrom)
        rnaSize = self._findRnaSize(gp)
        if gp.cdsStart < gp.cdsEnd:
            cdsChrom = Coords(gp.chrom, gp.cdsStart, gp.cdsEnd, strand='+', size=chromSize)
        else:
            cdsChrom = None

        chrom = Coords(gp.chrom, gp.txStart, gp.txEnd, '+', chromSize)
        rna = Coords(gp.name, 0, rnaSize, gp.strand, rnaSize)
        builder = AnnotFeatureBuilder(chrom, rna, cdsChrom, gp.strand, attrs, self.genomeReader)
        self._buildFeatures(gp, builder)
        builder.finish()
        return builder.transAnnot
