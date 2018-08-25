"""
Read GENCODE annotations from a database.
"""
import six
from pycbio.sys.objDict import ObjDict
from pycbio.db.sqliteOps import sqliteConnect
from pycbio.hgdata.genePredSqlite import GenePredSqliteTable
from pycbio.hgdata.gencodeSqlite import GencodeAttrsSqliteTable, GencodeTagSqliteTable
from gencode_icedb.general.genePredAnnotFeatures import AnnotationGenePredFactory
from gencode_icedb.general.geneAnnot import geneAnnotGroup

# FIXME: rename module to gencodeReader.py
# FIXME: transcriptType filter should be done in sql

# tables in sqlite databases
GENCODE_ANN_TABLE = "gencode_ann"
GENCODE_ATTRS_TABLE = "gencode_attrs"
GENCODE_TRANSCRIPT_SOURCE_TABLE = "gencode_transcript_source"
GENCODE_TRANSCRIPTION_SUPPORT_LEVEL_TABLE = "gencode_transcription_support_level"
GENCODE_TAG_TABLE = "gencode_tag"


def _isChrYPar(annotTrans):
    # ucsc data puts PAR tag on both chrX and chrY, since it tracks only by transcript id.
    return ("PAR" in annotTrans.attrs.tags) and (annotTrans.chrom.name in ("chrY", "Y"))


class UcscGencodeReader(object):
    """Object for accessing a GENCODE sqlite database with UCSC tables.
    """
    def __init__(self, gencodeDbFile, genomeReader=None, filterChrYPar=True,
                 transcriptTypes=None):
        self.conn = None
        self.conn = sqliteConnect(gencodeDbFile)
        self.filterChrYPar = filterChrYPar
        self.transcriptTypes = frozenset(transcriptTypes) if transcriptTypes is not None else None
        self.genePredDbTable = GenePredSqliteTable(self.conn, GENCODE_ANN_TABLE)
        self.attrDbTable = GencodeAttrsSqliteTable(self.conn, GENCODE_ATTRS_TABLE)
        self.tagDbTable = GencodeTagSqliteTable(self.conn, GENCODE_TAG_TABLE)
        self.annotFactory = AnnotationGenePredFactory(genomeReader)

    def close(self):
        self.conn.close()
        self.conn = None

    def _getAttrs(self, gp):
        attrs = ObjDict()
        attrs.update(self.attrDbTable.getrByTranscriptId(gp.name)._asdict())
        attrs.tags = frozenset([t.tag for t in self.tagDbTable.getByTranscriptId(gp.name)])
        return attrs

    def _makeTransAnnot(self, gp):
        "will return None if chrY PAR trans and these are being filtered"
        transAnnot = self.annotFactory.fromGenePred(gp, self._getAttrs(gp))
        if self.filterChrYPar and _isChrYPar(transAnnot):
            return None
        elif (self.transcriptTypes is not None) and (transAnnot.attrs.transcriptType not in self.transcriptTypes):
            return None
        else:
            return transAnnot

    def getGeneIds(self):
        return self.attrDbTable.getGeneIds()

    def getByGeneId(self, geneId):
        transIds = self.attrDbTable.getGeneTranscriptIds(geneId)
        return self.getByTranscriptIds(transIds)

    def getByTranscriptIds(self, transIds):
        """Get an annotation as a TranscriptFeatures by transcript id (or ids).
        """
        if isinstance(transIds, six.string_types):
            transIds = [transIds]
        transAnnots = []
        for transId in transIds:
            for gp in self.genePredDbTable.getByName(transId):
                transAnnot = self._makeTransAnnot(gp)
                if transAnnot is not None:
                    transAnnots.append(transAnnot)
        return transAnnots

    def _getByGencodeId(self, gencodeId):
        """get annotations as TranscriptFeatures by a gene or transcript id"""
        attrs = self.attrDbTable.getByGeneId(gencodeId)
        if len(attrs) > 0:
            return self.getByGeneId(gencodeId)
        else:
            attrs = self.attrDbTable.getByTranscriptId(gencodeId)
            if attrs is None:
                raise Exception("Not a valid GENCODE gene or transcript id: {}".format(gencodeId))
            return self.getByTranscriptIds(gencodeId)

    def _getByGencodeIds(self, gencodeIds):
        """get annotations as TranscriptFeatures by gene or transcript id (or ids)"""
        if isinstance(gencodeIds, six.string_types):
            gencodeIds = [gencodeIds]
        transAnnots = []
        for gencodeId in gencodeIds:
            transAnnots.extend(self._getByGencodeId(gencodeId))
        return transAnnots

    def getGenesByGencodeIds(self, gencodeIds):
        """get annotations as GeneAnnotation object, containing the gene's
        TranscriptFeatures.  This can be by gene or transcript id (or list of
        ids)"""
        return geneAnnotGroup(self._getByGencodeIds(gencodeIds))

    def getTranscriptsOverlapping(self, chrom, start, end, strand=None):
        # FIXME: switch to chords
        """generator get overlapping annotations as TranscriptFeatures"""
        transAnnots = []
        for gp in self.genePredDbTable.getRangeOverlap(chrom, start, end, strand=strand):
            transAnnot = self._makeTransAnnot(gp)
            if transAnnot is not None:
                transAnnots.append(transAnnot)
        return transAnnots
