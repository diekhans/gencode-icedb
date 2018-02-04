from pycbio.hgdata.hgLite import sqliteConnect, GenePredDbTable, GencodeAttrsDbTable, GencodeTranscriptSourceDbTable, GencodeTranscriptionSupportLevelDbTable, GencodeTagDbTable
from gencode_icedb.general.annotFeatures import AnnotationGenePredDbFactory


# tables in sqlite databases
gencode_ann_table = "gencode_ann"
gencode_attrs_table = "gencode_attrs"
gencode_transcript_source_table = "gencode_transcript_source"
gencode_transcription_support_level_table = "gencode_transcription_support_level"
gencode_tag_table = "gencode_tag"


class GencodeDb(object):
    """collection of access objects for a GENCODE sqlite database"""
    def __init__(self, gencodeDb, genomeReader):
        self.conn = sqliteConnect(gencodeDb)
        self.annotFactory = AnnotationGenePredDbFactory(self.conn, gencode_ann_table, genomeReader)
        self.attrsTbl = GencodeAttrsDbTable(self.conn, gencode_attrs_table)
        self.transSourceTbl = GencodeTranscriptSourceDbTable(self.conn, gencode_transcript_source_table)
        self.tslTbl = GencodeTranscriptionSupportLevelDbTable(self.conn, gencode_transcription_support_level_table)
        self.tagTbl = GencodeTagDbTable(self.conn, gencode_tag_table)

        #FIXME: tslTbl not normally neded
