"""
Read GENCODE annotations from a sqllite3 database with UCSC browser GENCODE
tables.
"""
from pycbio.sys.objDict import ObjDict
from pycbio.db import sqliteOps
from pycbio.hgdata.genePredSqlite import GenePredSqliteTable
from pycbio.hgdata.gencodeSqlite import GencodeAttrsSqliteTable, GencodeTagSqliteTable
from gencode_icedb.general.dataOps import ensureList, isChrYPar, ensemblIdSplit
from gencode_icedb.general.genePredAnnotFeatures import GenePredAnnotationFactory
from gencode_icedb.general.geneAnnot import geneAnnotGroup


# tables in sqlite databases
GENCODE_ANN_TABLE = "gencode_ann"
GENCODE_ATTRS_TABLE = "gencode_attrs"
GENCODE_TRANSCRIPT_SOURCE_TABLE = "gencode_transcript_source"
GENCODE_TRANSCRIPTION_SUPPORT_LEVEL_TABLE = "gencode_transcription_support_level"
GENCODE_TAG_TABLE = "gencode_tag"
GENCODE_GENE_TABLE = "gencode_gene"


class UcscGencodeReader(object):
    """Object for accessing a GENCODE sqlite database with UCSC tables.
    """
    def __init__(self, gencodeDbFile, genomeReader=None, filterChrYPar=True,
                 transcriptTypes=None):
        self.conn = sqliteOps.connect(gencodeDbFile)
        self.filterChrYPar = filterChrYPar
        self.transcriptTypes = frozenset(transcriptTypes) if transcriptTypes is not None else None
        self.genePredDbTable = GenePredSqliteTable(self.conn, GENCODE_ANN_TABLE)
        self.attrDbTable = GencodeAttrsSqliteTable(self.conn, GENCODE_ATTRS_TABLE)
        self.tagDbTable = GencodeTagSqliteTable(self.conn, GENCODE_TAG_TABLE)
        self.annotFactory = GenePredAnnotationFactory(genomeReader)

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
        if self.filterChrYPar and isChrYPar(transAnnot):
            return None
        elif (self.transcriptTypes is not None) and (transAnnot.attrs.transcriptType not in self.transcriptTypes):
            return None
        else:
            return transAnnot

    def getGeneIds(self):
        return sorted(self.attrDbTable.getGeneIds())

    def _getByTranscriptId(self, transId):
        gps = tuple(self.genePredDbTable.getByName(transId))
        if len(gps) == 0:
            raise Exception("no transcripts with id {}".format(transId))
        for gp in gps:
            transAnnot = self._makeTransAnnot(gp)
            if transAnnot is not None:
                yield transAnnot

    def getByTranscriptIds(self, transIds):
        """Get an annotation as a TranscriptFeatures by transcript id or ids.  For PAR regions,
        multiple transcript can be returned.
        """
        transIds = ensureList(transIds)
        transAnnots = []
        for transId in transIds:
            transAnnots.extend(self._getByTranscriptId(transId))
        return transAnnots

    def getByGeneId(self, geneId):
        "get all transcripts associated with a geneId"
        transIds = tuple(self.attrDbTable.getGeneTranscriptIds(geneId))
        if len(transIds) == 0:
            raise Exception("no genes with id {}".format(geneId))
        transAnnots = []
        for transId in transIds:
            transAnnots.extend(self._getByTranscriptId(transId))
        return transAnnots

    def _getByGencodeIds(self, gencodeIds):
        """get annotations as TranscriptFeatures by gene or transcript id (or ids)"""
        gencodeIds = ensureList(gencodeIds)
        geneIds, transIds = ensemblIdSplit(gencodeIds)
        for geneId in geneIds:
            yield from self.getByGeneId(geneId)
        for transId in transIds:
            yield from self._getByTranscriptId(transId)

    def getGenesByGencodeIds(self, gencodeIds):
        """get annotations as GeneAnnotation object, containing the gene's
        TranscriptFeatures.  This can be by gene or transcript id (or list of
        ids)"""
        return geneAnnotGroup(self._getByGencodeIds(gencodeIds))

    def getTranscriptsOverlapping(self, chrom, start, end, strand=None):
        """generator get overlapping annotations as TranscriptFeatures"""
        # FIXME: switch to chords
        transAnnots = []
        for gp in self.genePredDbTable.getRangeOverlap(chrom, start, end, strand=strand):
            transAnnot = self._makeTransAnnot(gp)
            if transAnnot is not None:
                transAnnots.append(transAnnot)
        return transAnnots
