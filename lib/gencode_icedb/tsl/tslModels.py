"""
PeeWee data models for RNA-Seq metadata and splice junctions.
"""
from peewee import Proxy, Model, Field, PrimaryKeyField, CharField, IntegerField, UUIDField
from gencode_icedb.general.peeweeOps import peeweeConnect, peeweeClose, peeweeClassToTableName, PeeweeModelMixins
from gencode_icedb.tsl.supportDefs import TrascriptionSupportLevel, EvidenceSupport

_database_proxy = Proxy()


def setDatabaseConn(dbconn):
    "bind the proxy to a database"
    _database_proxy.initialize(dbconn)


def tslConnect(sqliteDb, create=False, readonly=True, timeout=None, synchronous=None):
    "connect to sqlite3 database and bind to model"
    return peeweeConnect(sqliteDb, setDatabaseConn, create=create, readonly=readonly, timeout=timeout, synchronous=synchronous)


def tslClose(conn):
    "close database"
    # FIXME should reset proxy
    peeweeClose(conn)


class TrascriptionSupportLevelField(Field):
    """A field storing a TrascriptionSupportLevel python enumeration as a string in the database."""
    db_field = 'text'

    def db_value(self, value):
        return str(value)

    def python_value(self, value):
        return TrascriptionSupportLevel(value)


class EvidenceSupportField(Field):
    """A field storing an EvidenceSupport python enumeration as a string in the database."""
    db_field = 'text'

    def db_value(self, value):
        return str(value)

    def python_value(self, value):
        return EvidenceSupport(value)


class BaseModel(Model, PeeweeModelMixins):
    "base for peewee models, used to bind proxy"
    class Meta:
        database = _database_proxy
        table_function = peeweeClassToTableName


class GencodeSupportEval(BaseModel):
    """Support for GENCODE transcripts from various evidence set.
    """
    id = PrimaryKeyField()
    transcriptId = CharField(index=True,
                             help_text="""GENCODE trancript id""")
    evidSetUuid = UUIDField(index=True,
                            help_text="""UUID of the evidence set""")
    support = EvidenceSupportField(help_text="""Best support by this evidence set""")
    evidCount = IntegerField(help_text="""Number of evidence alignments at this support level.""")
    offset5 = IntegerField(help_text="""Offset of the start of the first evidence exon to the start of """
                           """annotation.  Positive indicates the evidence is longer than annotation, """
                           """negative is shorter. This is in the direction of transcription.  This is """
                           """the longest support of any supporting evidence.""")
    offset3 = IntegerField(help_text="""Offset of the end of the first evidence exon to the end of """
                           """annotation.  Interpretation is the same as offset3.""")
    extend5Exons = IntegerField(help_text="""Maximum extending 5' number of exons.""")
    extend3Exons = IntegerField(help_text="""Maximum extending 3' number of exons.""")


class GencodeTranscriptSupport(BaseModel):
    """Transcript-level support"""
    id = PrimaryKeyField()
    transcriptId = CharField(index=True,
                             help_text="""GENCODE trancript id""")
    level = TrascriptionSupportLevelField(index=True,
                                          help_text="""GENCODE TSL""")
    intLevel = IntegerField(index=True,
                            help_text="""GENCODE TSL as an integer""")
