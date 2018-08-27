"""
create AnnotationFeature objects from Ensembl database records (see ensemblDb.py)
"""
from typing import List
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from gencode_icedb.general.ensemblDb import EnsemblGeneTrans, EnsemblTransExon, EnsemblTransAttr
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, Utr5RegionFeature, CdsRegionFeature, Utr3RegionFeature, GapAnnotFeature, NonCodingRegionFeature

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
        features = []
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
        features = []
        rnaNext = rnaStart
        chromNext = exonRecs[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            if chromNext != exonRecs[iBlk].start:
                rnaNext = self._addGapAnnot(chromNext, exonRecs[iBlk].start, rnaNext, exon, features)
            rnaNext = self._addBlkCodingFeatures(exonRecs[iBlk], rnaNext, exon, features)
            chromNext = features[-1].chrom.end
        assert rnaNext == rnaEnd
        exon.annotFeatures = tuple(features)
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _addNonCodingFeature(self, blk, rnaNext, exon, features):
        feat = NonCodingRegionFeature(exon, len(features),
                                      exon.chrom.subrange(blk.start, blk.end),
                                      exon.rna.subrange(rnaNext, rnaNext + blk.size()))
        features.append(feat)
        return feat.rna.end

    def _addNonCodingFeatures(self, geneTransRec, exonRecs, iBlkStart, iBlkEnd, rnaStart, rnaEnd, exon):
        features = []
        rnaNext = rnaStart
        chromNext = exonRecs[iBlkStart].start
        for iBlk in range(iBlkStart, iBlkEnd):
            if chromNext != exonRecs[iBlk].start:
                rnaNext = self._addGapAnnot(chromNext, exonRecs[iBlk].start, rnaNext, exon, features)
            rnaNext = self._addNonCodingFeature(exonRecs[iBlk], rnaNext, exon, features)
            chromNext = features[-1].chrom.end
        exon.annotFeatures = tuple(features)
        assert rnaNext == rnaEnd, "rnaNext={} != rnaEnd={}".format(rnaNext, rnaEnd)

    def _makeExon(self, gp, iBlkStart, iBlkEnd, rnaStart, rnaEnd, transAnnot, iFeat):
        exon = ExonFeature(transAnnot, iFeat,
                           transAnnot.chrom.subrange(exonRecs[iBlkStart].start, exonRecs[iBlkEnd - 1].end),
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
                                      exonRecs[iBlkNext - 1].end,
                                      exonRecs[iBlkNext].start, gp.strand)

    def _makeIntron(self, geneTransRec, exonRecs, iBlkNext, rnaEnd, transAnnot, iFeat):
        donorSeq, acceptorSeq = self._getSpliceSites(gp, iBlkNext)
        return IntronFeature(trans, iFeat,
                             transAnnot.chrom.subrange(exonRecs[iBlkNext - 1].end, exonRecs[iBlkNext].start),
                             transAnnot.rna.subrange(rnaEnd, rnaEnd), donorSeq, acceptorSeq)

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

        chrom = Coords(geneTransRec.chrom, geneTransRec.start, geneTransRec.end, '+', geneTransRec.chromSize)
        rna = Coords(geneTransRec.transcriptionId, 0, rnaSize, geneTransRec.strand, rnaSize)
        builder = AnnotFeatureBuilder(chrom, rna, cdsChrom, transcriptionStrand, attrs, self.genomeReader)
        self._buildFeatures(geneTransRec, exonRecs, builder)
        builder.finish()
        return build.transAnnot
