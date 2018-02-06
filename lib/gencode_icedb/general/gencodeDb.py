"""
Read GENCODE annotations from a database.
"""

from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.hgLite import sqliteConnect, GenePredDbTable, GencodeAttrsDbTable, GencodeTagDbTable
from gencode_icedb.general.annotFeatures import AnnotationGenePredFactory

# tables in sqlite databases
GENCODE_ANN_TABLE = "gencode_ann"
GENCODE_ATTRS_TABLE = "gencode_attrs"
GENCODE_TRANSCRIPT_SOURCE_TABLE = "gencode_transcript_source"
GENCODE_TRANSCRIPTION_SUPPORT_LEVEL_TABLE = "gencode_transcription_support_level"
GENCODE_TAG_TABLE = "gencode_tag"


class UcscGencodeReader(object):
    """Object for accessing a GENCODE sqlite database with UCSC tables """
    def __init__(self, gencodeDbFile, genomeReader):
        self.conn = sqliteConnect(gencodeDbFile)
        self.genePredDbTable = GenePredDbTable(self.conn, GENCODE_ANN_TABLE)
        self.attrDbTable = GencodeAttrsDbTable(self.conn, GENCODE_ATTRS_TABLE)
        self.tagDbTable = GencodeTagDbTable(self.conn, GENCODE_TAG_TABLE)
        self.annotFactory = AnnotationGenePredFactory(genomeReader)

    def close(self):
        self.conn.close()
        self.conn = None

    def __del__(self):
        if self.conn is not None:
            self.close()

    def _getMetaData(self, gp):
        metaData = ObjDict()
        metaData.attrs = self.attrDbTable.getrByTranscriptId(gp.name)
        metaData.tags = frozenset([t.tag for t in self.tagDbTable.getByTranscriptId(gp.name)])
        return metaData

    def _makeTransAnnon(self, gp):
        return self.annotFactory.fromGenePred(gp, self._getMetaData(gp))

    def getGeneIds(self):
        return self.attrDbTable.getGeneIds()

    def getByGeneId(self, geneId):
        transIds = list(self.attrDbTable.getGeneTranscriptIds(geneId))
        return list(self.byNameGen(transIds))

    def byNameGen(self, names):
        """generator get annotations as TranscriptFeatures by name"""
        # FIXME: optimize to one query, need gene name query instead
        for name in names:
            for gp in self.genePredDbTable.getByName(name):
                yield self._makeTransAnnon(gp)

    def overlappingGen(self, chrom, start, end, strand=None):
        """generator get overlapping annotations as TranscriptFeatures"""
        for gp in self.genePredDbTable.getRangeOverlap(chrom, start, end, strand=strand):
            yield self._makeTransAnnon(gp)

    def allGen(self):
        """generator get overlapping annotations as TranscriptFeatures"""
        for gp in self.genePredDbTable.getAll():
            yield self._makeTransAnnon(gp)
