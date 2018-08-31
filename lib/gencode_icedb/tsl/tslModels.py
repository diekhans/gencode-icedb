"""
PeeWee data models for RNA-Seq metadata and splice junctions.
"""
from peewee import Proxy, Model, Field, PrimaryKeyField, CharField, IntegerField
from gencode_icedb.general.peeweeOps import peeweeConnect, peeweeClose, peeweeClassToTableName, PeeweeModelMixins
from gencode_icedb.tsl.evidenceDb import EvidenceSource
from gencode_icedb.tsl.supportDefs import TrascriptionSupportLevel, EvidenceSupport, GenbankProblemReason

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


class EvidenceSourceField(Field):
    """A field storing an EvidenceSource python enumeration as a string in the database."""
    db_field = 'text'

    def db_value(self, value):
        return str(value)

    def python_value(self, value):
        return EvidenceSource(value)


class EvidenceSupportField(Field):
    """A field storing an EvidenceSupport python enumeration as a string in the database."""
    db_field = 'text'

    def db_value(self, value):
        return str(value)

    def python_value(self, value):
        return EvidenceSupport(value)


class GenbankProblemReasonField(Field):
    """A field storing an GenbankProblemReason python enumeration as a string in the database."""
    db_field = 'text'

    def db_value(self, value):
        return None if value is None else str(value)

    def python_value(self, value):
        if (value is None) or (value is ''):
            return None
        else:
            return GenbankProblemReason(value)


class BaseModel(Model, PeeweeModelMixins):
    "base for peewee models, used to bind proxy"
    class Meta:
        database = _database_proxy
        table_function = peeweeClassToTableName


class GencodeTranscriptSupport(BaseModel):
    """Transcript-level support"""
    id = PrimaryKeyField()
    transcriptId = CharField(index=True,
                             help_text="""GENCODE trancript id""")
    level = TrascriptionSupportLevelField(index=True,
                                          help_text="""GENCODE TSL""")
    intLevel = IntegerField(index=True,
                            help_text="""GENCODE TSL as an integer""")


class GencodeTranscriptSupportDetails(BaseModel):
    """Transcript-level support"""
    id = PrimaryKeyField()
    transcriptId = CharField(index=True,
                             help_text="""GENCODE trancript id""")
    evidSrc = EvidenceSourceField(index=True,
                                  help_text="""Evidence source""")
    evidId = CharField(index=True,
                       help_text="""Evidence identifier""")
    evidSupport = EvidenceSupportField(index=False,
                                       help_text="""Support provided by this piece of evidence""")
    suspect = GenbankProblemReasonField(index=False, null=True,
                                        help_text="""If not NULL evidence is suspect with this reason""")
