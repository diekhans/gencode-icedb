"""
Common functions to support analysis.
"""
from __future__ import print_function
import sys
from numpy import median
from gencode_icedb.general import gencodeDb
from gencode_icedb.rsl.rslModel import GencodeSupport, GencodeNovel
from gencode_icedb.rsl.gencodeIntronEvid import IntronSupportLevel, intronEvidSupportLevel
import sqlite3
from gencode_icedb.rsl.rslModel import sqliteConnect, sqliteClose
from pycbio.sys.color import Color
from pycbio.hgdata.hgLite import GencodeTranscriptSourceDbTable, GencodeTranscriptionSupportLevelDbTable, GencodeTagDbTable
from gencode.data.gencodeGenes import sourceToExtendedMethod
from gencode_icedb.tsl.supportDefs import TrascriptionSupportLevel


def isPseudo(suppRec):
    # FIXME: switch to using functions now in CCDS tree
    return (suppRec.transcriptType.find("pseudogene") >= 0) and (suppRec.transcriptType != "polymorphic_pseudogene")


def _intronSupportReaderFilter(isNovelDb, rec):
    """for novel, don't keep ones generated from GENCODE, for genes, don't
    use pseudo"""
    if isNovelDb:
        return rec.numExprs > 0
    else:
        return not isPseudo(rec)


def intronSupportReader(isNovelDb, chrom=None):
    """read with filtering introns support from database.  Should have already
    connected and bound ORM.  If chrom is supplied, only query
    that chrom for testing purposes."""
    ormCls = GencodeNovel if isNovelDb else GencodeSupport
    query = ormCls.select()
    if chrom is not None:
        query = query.where(ormCls.chrom == chrom)
    for rec in query:
        if _intronSupportReaderFilter(isNovelDb, rec):
            yield rec


def intronSupportAlreadyProcessed(rec, processed):
    """check in an intron has already been processed; handles duplication by gene"""
    key = (rec.chrom, rec.intronStart, rec.intronEnd, rec.strand)
    if key in processed:
        return True
    else:
        processed.add(key)
        return False


class SupportTrackColors(object):
    "color codes"
    # FIXME: not just for tracks
    strong = Color.fromRgb8(0, 128, 0)  # dark green
    medium = Color.fromRgb8(0, 0, 255)  # light blue
    weak = Color.fromRgb8(238, 118, 0)  # dark orange
    none = Color.fromRgb8(200, 0, 0)  # red

    @staticmethod
    def __printColorHtml(fh, desc, color):
        print('<li><p style="color:{}">{}</p>'.format(color.toHtmlColor(), desc), file=fh)

    @staticmethod
    def printHtml(fh=sys.stdout):
        "for inclusion in documentation"
        SupportTrackColors.__printColorHtml(fh, "strong support", SupportTrackColors.strong)
        SupportTrackColors.__printColorHtml(fh, "medium support", SupportTrackColors.medium)
        SupportTrackColors.__printColorHtml(fh, "weak support", SupportTrackColors.weak)
        SupportTrackColors.__printColorHtml(fh, "no support", SupportTrackColors.none)

    @staticmethod
    def supportLevelColor(level):
        if level == IntronSupportLevel.STRONG:
            return SupportTrackColors.strong
        elif level == IntronSupportLevel.MEDIUM:
            return SupportTrackColors.medium
        elif level == IntronSupportLevel.WEAK:
            return SupportTrackColors.weak
        elif level == IntronSupportLevel.NONE:
            return SupportTrackColors.none


class TransEvid(object):
    def __init__(self, geneId, transcriptId, bioType):
        self.geneId = geneId
        self.transcriptId = transcriptId
        self.bioType = bioType
        self.method = None
        self.tsl = None
        self.tags = set()
        self.introns = []
        self.intronLevels = []

    def addIntron(self, intron):
        self.introns.append(intron)
        self.intronLevels.append(intronEvidSupportLevel(intron.numUniqueMapReads))

    @property
    def chrom(self):
        return self.introns[0].chrom

    @property
    def strand(self):
        return self.introns[0].strand

    def countLevels(self):
        levelCnts = len(IntronSupportLevel) * [0]
        for intronLevel in self.intronLevels:
            levelCnts[intronLevel.value] += 1
        return levelCnts

    def fullSupportLevel(self):
        "finds lowest level any intron, which is the full transcript support level"
        levelCnts = self.countLevels()
        for level in IntronSupportLevel:
            if levelCnts[level.value] > 0:
                return level

    def medianSupportLevel(self):
        "finds lowest level any intron, which is the full transcript support level"
        return IntronSupportLevel(int(median([il.value for il in self.intronLevels])))


class GencodeIntronEvid(object):
    "data linked to GENCODE transcripts"

    def __init__(self):
        self.byTransIds = {}
        self.chroms = set()
        self.bioTypes = set()

    def __obtainTransInfo(self, rec):
        trans = self.byTransIds.get(rec.transcriptId)
        if trans is None:
            trans = self.byTransIds[rec.transcriptId] = TransEvid(rec.geneId, rec.transcriptId, rec.transcriptType)
        return trans

    def __loadRec(self, rec):
        trans = self.__obtainTransInfo(rec)
        trans.addIntron(rec)
        self.chroms.add(rec.chrom)
        self.bioTypes.add(trans.bioType)

    def loadSupport(self, chrom=None):
        for rec in intronSupportReader(False, chrom):
            self.__loadRec(rec)

    def loadSupportDb(self, intronEvidDb, chrom=None):
        conn = sqliteConnect(intronEvidDb, readonly=True)
        self.loadSupport(chrom)
        sqliteClose(conn)

    def loadGencodeSource(self, gencodeDbConn):
        dbTable = GencodeTranscriptSourceDbTable(gencodeDbConn, gencodeDb.gencode_transcript_source_table)
        for rec in dbTable.queryAll():
            trans = self.byTransIds.get(rec.transcriptId)
            if trans is not None:
                trans.method = sourceToExtendedMethod(rec.source)

    def loadGencodeSupportLevel(self, gencodeDbConn):
        dbTable = GencodeTranscriptionSupportLevelDbTable(gencodeDbConn, gencodeDb.gencode_transcription_support_level_table)
        for rec in dbTable.queryAll():
            trans = self.byTransIds.get(rec.transcriptId)
            if trans is not None:
                trans.tsl = TrascriptionSupportLevel(rec.level)

    def loadGencodeTags(self, gencodeDbConn):
        dbTable = GencodeTagDbTable(gencodeDbConn, gencodeDb.gencode_tag_table)
        for rec in dbTable.queryAll():
            trans = self.byTransIds.get(rec.transcriptId)
            if trans is not None:
                trans.tags.add(rec.tag)

    def loadGencodeDb(self, gencodeDb):
        gencodeDbConn = sqlite3.connect(gencodeDb)
        self.loadGencodeSource(gencodeDbConn)
        self.loadGencodeSupportLevel(gencodeDbConn)
        self.loadGencodeTags(gencodeDbConn)
        gencodeDbConn.close()
