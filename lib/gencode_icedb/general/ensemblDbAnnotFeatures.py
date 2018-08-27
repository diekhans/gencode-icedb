"""
create AnnotationFeature objects from Ensembl database records (see ensemblDb.py)
"""
from typing import List
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from gencode_icedb.general.ensemblDb import EnsemblGeneTrans, EnsemblTransExon, EnsemblTransAttr
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.annotFeatureBuilder import AnnotFeatureBuilder

# FIXME: consistent variable naming: annot, annotTrans transAnnot ..


class EnsemblDbAnnotationFactory(object):
    """
    Factory to create annotation features from the Ensembl database..
    """
    def __init__(self, genomeReader=None):
        """genomeReader is used to obtain splice sites and chrom size, maybe
        None if splice sites and reverse complement will not be done. To get chrom sizes
        without splice sites, provide the chromSizeFunc(chrom) function."""
        self.genomeReader = genomeReader

    def _buildFeatures(self, geneTransRec, exonRecs, builder):
        iBlkStart = 0
        rnaStart = 0
        while iBlkStart < len(exonRecs):
            iBlkEnd, qCount = self._findExonEnd(exonRecs, iBlkStart)
            rnaEnd = rnaStart + qCount
            self._makeExon(geneTransRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
            if iBlkEnd < len(exonRecs):
                self._makeIntron(geneTransRec, exonRecs, iBlkEnd, rnaEnd, builder)
            rnaStart = rnaEnd
            iBlkStart = iBlkEnd

    def _tGapSize(self, exonRecs, iBlk):
        "size of gap before the block"
        return exonRecs[iBlk].start - exonRecs[iBlk - 1].end

    def _findExonEnd(self, exonRecs, iBlkStart):
        """finds half-open end of blocks covering current exon and total base,
        including closed gaps"""
        iBlkEnd = iBlkStart + 1
        while (iBlkEnd < len(exonRecs)) and (self._tGapSize(exonRecs, iBlkEnd) < minIntronSize):
            iBlkEnd += 1
        return iBlkEnd, exonRecs[iBlkEnd - 1].end - exonRecs[iBlkStart].start

    def _findRnaSize(self, exonRecs):
        iBlkStart = 0
        rnaSize = 0
        while iBlkStart < len(exonRecs):
            iBlkEnd, qCount = self._findExonEnd(exonRecs, iBlkStart)
            rnaSize += qCount
            iBlkStart = iBlkEnd
        return rnaSize

    def _addCodingFeatures(self, geneTransRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        rnaNext = rnaStart
        chromNext = exonRecs[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            exonRec = exonRecs[iBlk]
            if chromNext != exonRecs[iBlk].start:
                rnaNext = builder.addGap(chromNext, exonRec.start, rnaNext)
                chromNext = exonRec.start
            frame = Frame(exonRec.phase) if exonRec.phase >= 0 else None
            rnaNext = builder.addBlkCodingFeatures(chromNext, exonRec.end, rnaNext, frame)
            chromNext = exonRecs[iBlk].end
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _addNonCodingFeatures(self, geneTransRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        rnaNext = rnaStart
        chromNext = exonRecs[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            exonRec = exonRecs[iBlk]
            if chromNext != exonRec.start:
                rnaNext = builder.addGap(chromNext, exonRec.start, rnaNext)
                chromNext = exonRec.start
            rnaNext = builder.addNonCodingFeature(chromNext, exonRec.end, rnaNext)
            chromNext = exonRec.end
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _makeExon(self, geneTransRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        builder.beginExon(exonRecs[iBlkStart].start, exonRecs[iBlkEnd - 1].end,
                          rnaStart, rnaEnd)
        if builder.hasCds:
            self._addCodingFeatures(geneTransRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
        else:
            self._addNonCodingFeatures(geneTransRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
        builder.finishExon()

    def _makeIntron(self, geneTransRec, exonRecs, iBlkNext, rnaEnd, builder):
        builder.addIntron(exonRecs[iBlkNext - 1].end, exonRecs[iBlkNext].start, rnaEnd)

    def _getCdsCoords(self, geneTransRec, exonRecs):
        startExon = endExon = None
        for exonRec in exonRecs:
            if exonRec.exonDbId == geneTransRec.cdsStartExonDbId:
                startExon = exonRec
            if exonRec.exonDbId == geneTransRec.cdsEndExonDbId:
                endExon = exonRec
        cdsStart = startExon.start + geneTransRec.cdsStartOffset
        cdsEnd = endExon.start + geneTransRec.cdsEndOffset
        if geneTransRec.strand == '-':
            cdsStart, cdsEnd = cdsEnd, cdsStart
        return Coords(geneTransRec.chrom, cdsStart, cdsEnd, strand='+', size=geneTransRec.chromSize)

    def fromEnsemblDb(self, geneTransRec: EnsemblGeneTrans,
                      exonRecs: List[EnsemblTransExon],
                      attrsRecs: List[EnsemblTransAttr]):
        """convert records from Ensembl database to TranscriptFeatures object."""
        exonRecs = tuple(sorted(exonRecs, key=lambda e: e.start))
        rnaSize = self._findRnaSize(exonRecs)
        if geneTransRec.cdsStartExonDbId is not None:
            cdsChrom = self._getCdsCoords(geneTransRec, exonRecs)
        else:
            cdsChrom = None

        attrs = None#FIXME
        chrom = Coords(geneTransRec.chrom, geneTransRec.start, geneTransRec.end, '+', geneTransRec.chromSize)
        rna = Coords(geneTransRec.transcriptionId, 0, rnaSize, geneTransRec.strand, rnaSize)
        builder = AnnotFeatureBuilder(chrom, rna, cdsChrom, geneTransRec.strand, attrs, self.genomeReader)
        self._buildFeatures(geneTransRec, exonRecs, builder)
        builder.finish()
        return builder.transAnnot
