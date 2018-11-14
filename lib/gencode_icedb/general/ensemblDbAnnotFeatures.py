"""
create AnnotationFeature objects from Ensembl database records (see ensemblDbQuery.py)
"""
from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from gencode_icedb.general.ensemblDbQuery import EnsemblTranscript
from gencode_icedb.tsl import minIntronSize
from gencode_icedb.general.annotFeatureBuilder import AnnotFeatureBuilder

# FIXME: consistent variable naming: annot, annotTrans transAnnot ..


###
# Representation of frame/phase in Ensembl database.
#
# Referred to asphase, this is not the same as GFF/GTF phgase. Is what is sometimes
# called the frame and has the value of 0, 1, 2.
#
# From Ensembl Bio::EnsEMBL::Exon:
#   * phase - The Ensembl phase convention can be thought of as"the number of
#     bases of the first codon which are on the previous exon".  It is
#     therefore 0, 1 or 2 (or -1 if the exon is non-coding).
#   * end_phase - The number of bases from the last incomplete codon of this
#     exon. Usually, end_phase = (phase + exon_length)%3 but end_phase could
#     be -1 if the exon is half-coding and its 3 prime end is UTR.
#
# Description of how phase is stored in database:
#
#    * translation.start_exon_id references to 5' CDS exon.
#    * translation.seq_start is the offset of 5' codon in the 5' CDS exon:
#      - positive strand: this is the offset from the start of the exon
#      - negative strand: this is the offset from the end of the exon
#      first coding base in the 5' CDS exon.
#    * translation.end_exon_id references the 3' CDS exon.
#    * translation.seq_end is the offset from 3' codon in the 3'CDS exon:
#      - Positive strand this is the offset from the start of the exon
#      - Negative strand this is the offset from the end of the exon
#    * exon.phase has the starting phase in the direction of transcription of a
#      CDS exon if the first base of the exon is in a codon.  It has -1 if the
#      first base is in UTR.  If this is the 5' CDS exon, then assumed to be zero.
#      assumed to be 0, there is no way to represent a partial start codon
#      that is not at the start of an exon.
#    - exon.end_phase has the end phase in the direction of a CDS exon if the
#      last base of the exon is in a codon.  It has -1 if the last base is
#      in UTR.
#
#   Note:
#    * Exons not referenced by start_exon_id or end_exon_id and with
#      exon.phase and exon.end_phase -1 are UTR or non-coding exons.
#    * The above rules require the 5' codon to either start at the being
#      of the exon or have a phase of 0, it is not possible to represent a partial
#      initial codon that is not a the start.  It is possible to represent a partial
#      terminal codon that isn't at the edge of an exon by computing the frame from the
#      start. See polymorphic psueudogene ENST00000377712.3.
#
#
##

class EnsemblDbAnnotationFactory(object):
    """
    Factory to create annotation features from the Ensembl database..
    """
    def __init__(self, genomeReader=None):
        """genomeReader is used to obtain splice sites and chrom size, maybe None if
        splice sites and reverse complement will not be done. To get chrom
        sizes without splice sites, provide the chromSizeFunc(chrom)
        function.
        """
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

    def _findCdsExonRange(self, transRec, exonRecs):
        startExon = endExon = None
        for exonRec in exonRecs:
            if exonRec.exonDbId == transRec.cdsStartExonDbId:
                startExon = exonRec
            if exonRec.exonDbId == transRec.cdsEndExonDbId:
                endExon = exonRec
        assert (startExon is not None) and (endExon is not None)
        return (startExon, endExon)

    def _getCdsCoords(self, transRec, exonRecs):
        startExon, endExon = self._findCdsExonRange(transRec, exonRecs)
        if transRec.strand == '+':
            cdsStart = startExon.start + transRec.cdsStartOffset
            cdsEnd = endExon.start + transRec.cdsEndOffset
        else:
            cdsStart = endExon.end - transRec.cdsEndOffset
            cdsEnd = startExon.end - transRec.cdsStartOffset
        return Coords(transRec.chrom, cdsStart, cdsEnd, strand='+', size=transRec.chromSize)

    def _getInitialCdsExonFrame(self, transRec, exonRec):
        if exonRec.startPhase >= 0:
            return Frame(exonRec.startPhase)
        else:
            return Frame(0)

    def _getTerminalCdsExonFrame(self, transRec, exonRec):
        if exonRec.startPhase >= 0:
            return Frame(exonRec.startPhase)
        elif exonRec.endPhase >= 0:
            if transRec.strand == '-':
                return Frame(exonRec.endPhase) - transRec.cdsEndOffset
            return None
        else:
            return None

    def _getInternalCdsExonFrame(self, transRec, exonRec):
        if exonRec.startPhase >= 0:
            return Frame(exonRec.startPhase)
        else:
            return None

    def _getExonFrame(self, transRec, exonRec):
        if exonRec.exonDbId == transRec.cdsStartExonDbId:
            return self._getInitialCdsExonFrame(transRec, exonRec)
        elif exonRec.exonDbId == transRec.cdsEndExonDbId:
            return self._getTerminalCdsExonFrame(transRec, exonRec)
        elif (exonRec.startPhase >= 0) or (exonRec.endPhase >= 0):
            return self._getInternalCdsExonFrame(transRec, exonRec)
        else:
            return None  # not a coding exon

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

    def _getAttrOrNone(self, code, ensTrans):
        for attr in ensTrans.attrs:
            if attr.code == code:
                return attr.value
        return None

    def _makeAttrs(self, ensTrans):
        "UCSC compatible attribute names (camel-case)"
        # see print_gencode_annotation.pl
        # FIXME are missing ones needed?
        attrs = ObjDict()
        attrs.geneId = ensTrans.transcript.geneId
        attrs.geneName = ensTrans.transcript.geneName
        attrs.geneType = ensTrans.transcript.geneType
        attrs.transcriptId = ensTrans.transcript.transcriptId
        # FIXME attrs.transcriptName
        attrs.transcriptType = ensTrans.transcript.transcriptType
        # FIXME attrs.havanaGeneId
        # FIXME attrs.havanaTranscriptId
        attrs.ccdsId = self._getAttrOrNone("ccds_transcript", ensTrans)
        # FIXME attrs.level
        # FIXME attrs.transcriptClass
        # FIXME attrs.proteinId  (in translation)
        attrs.tsl = self._getAttrOrNone("TSL", ensTrans)
        return attrs

    def _appris_to_tag(self, apprisVal):
        if apprisVal == "principal1":
            return "appris_principal_1"
        elif apprisVal == "principal2":
            return "appris_principal_2"
        elif apprisVal == "principal3":
            return "appris_principal_3"
        elif apprisVal == "principal4":
            return "appris_principal_4"
        elif apprisVal == "principal5":
            return "appris_principal_5"
        elif apprisVal == "alternative1":
            return "appris_alternative_1"
        elif apprisVal == "alternative2":
            return "appris_alternative_2"
        else:
            raise Exception("unknown appris attribute value: {}".format(apprisVal))

    def _makeTags(self, ensTrans):
        # see print_gencode_annotation.pl
        tags = set()
        if self._getAttrOrNone('gencode_basic', ensTrans) is not None:
            tags.add('basic')
        if self._getAttrOrNone('upstream_ATG', ensTrans) is not None:
            tags.add('upstream_ATG')
        apprisVal = self._getAttrOrNone("appris", ensTrans)
        if apprisVal is not None:
            tags.add(self._appris_to_tag(apprisVal))
        return frozenset(tags)

    def _getAttrs(self, ensTrans):
        attrs = self._makeAttrs(ensTrans)
        attrs.tags = self._makeTags(ensTrans)
        return attrs

    def fromEnsemblDb(self, ensTrans: EnsemblTranscript):
        """convert records from Ensembl database to TranscriptFeatures object."""
        transRec = ensTrans.transcript
        rnaSize = self._findRnaSize(ensTrans.exons)
        if transRec.cdsStartExonDbId is not None:
            cdsChrom = self._getCdsCoords(transRec, ensTrans.exons)
        else:
            cdsChrom = None

        attrs = self._getAttrs(ensTrans)
        chrom = Coords(transRec.chrom, transRec.start, transRec.end, '+', transRec.chromSize)
        rna = Coords(transRec.transcriptId, 0, rnaSize, transRec.strand, rnaSize)
        builder = AnnotFeatureBuilder(chrom, rna, cdsChrom, transRec.strand, attrs, self.genomeReader)
        self._buildFeatures(transRec, ensTrans.exons, builder)
        builder.finish()
        return builder.transAnnot
