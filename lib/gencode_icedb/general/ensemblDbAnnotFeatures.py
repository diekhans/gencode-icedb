"""
create AnnotationFeature objects from Ensembl database records (see ensemblDb.py)
"""
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from gencode_icedb.general.ensemblDb import EnsemblTranscript
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

    def _buildFeatures(self, transRec, exonRecs, builder):
        iBlkStart = 0
        rnaStart = 0
        while iBlkStart < len(exonRecs):
            iBlkEnd, qCount = self._findExonEnd(exonRecs, iBlkStart)
            rnaEnd = rnaStart + qCount
            self._makeExon(transRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
            if iBlkEnd < len(exonRecs):
                self._makeIntron(transRec, exonRecs, iBlkEnd, rnaEnd, builder)
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

    def _getExonFrame(self, transRec, exonRec):
        if (exonRec.startPhase < 0) and (exonRec.endPhase < 0):
            return None
        elif exonRec.startPhase >= 0:
            return Frame(exonRec.startPhase)
        else:
            return Frame(exonRec.endPhase) - transRec.cdsEndOffset

    def _addCodingFeatures(self, transRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        rnaNext = rnaStart
        chromNext = exonRecs[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            exonRec = exonRecs[iBlk]
            if chromNext != exonRecs[iBlk].start:
                rnaNext = builder.addGap(chromNext, exonRec.start, rnaNext)
                chromNext = exonRec.start
            rnaNext = builder.addCodingFeatures(chromNext, exonRec.end, rnaNext,
                                                self._getExonFrame(transRec, exonRec))
            chromNext = exonRecs[iBlk].end
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _addNonCodingFeatures(self, transRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
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

    def _makeExon(self, transRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder):
        builder.beginExon(exonRecs[iBlkStart].start, exonRecs[iBlkEnd - 1].end,
                          rnaStart, rnaEnd)
        if builder.hasCds:
            self._addCodingFeatures(transRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
        else:
            self._addNonCodingFeatures(transRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, builder)
        builder.finishExon()

    def _makeIntron(self, transRec, exonRecs, iBlkNext, rnaEnd, builder):
        builder.addIntron(exonRecs[iBlkNext - 1].end, exonRecs[iBlkNext].start, rnaEnd)

    def _getCdsCoords(self, transRec, exonRecs):
        startExon = endExon = None
        for exonRec in exonRecs:
            if exonRec.exonDbId == transRec.cdsStartExonDbId:
                startExon = exonRec
            if exonRec.exonDbId == transRec.cdsEndExonDbId:
                endExon = exonRec
        cdsStart = startExon.start + transRec.cdsStartOffset
        cdsEnd = endExon.start + transRec.cdsEndOffset
        if transRec.strand == '-':
            cdsStart, cdsEnd = cdsEnd, cdsStart
        return Coords(transRec.chrom, cdsStart, cdsEnd, strand='+', size=transRec.chromSize)

    def fromEnsemblDb(self, ensTrans: EnsemblTranscript):
        """convert records from Ensembl database to TranscriptFeatures object."""
        transRec = ensTrans.transcript
        ensTrans.dump()
        rnaSize = self._findRnaSize(ensTrans.exons)
        if transRec.cdsStartExonDbId is not None:
            cdsChrom = self._getCdsCoords(transRec, ensTrans.exons)
        else:
            cdsChrom = None

        attrs = None  #FIXME
        chrom = Coords(transRec.chrom, transRec.start, transRec.end, '+', transRec.chromSize)
        rna = Coords(transRec.transcriptId, 0, rnaSize, transRec.strand, rnaSize)
        builder = AnnotFeatureBuilder(chrom, rna, cdsChrom, transRec.strand, attrs, self.genomeReader)
        self._buildFeatures(transRec, ensTrans.exons, builder)
        builder.finish()
        return builder.transAnnot
