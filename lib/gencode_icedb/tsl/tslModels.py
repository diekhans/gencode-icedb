"""
PeeWee data models for RNA-Seq metadata and splice junctions.
"""
from __future__ import print_function
from peewee import Proxy, Model, Field, PrimaryKeyField, CharField
from gencode_icedb.general.peeweeOps import peeweeConnect, peeweeClose
from gencode_icedb.tsl.supportDefs import TrascriptionSupportLevel, EvidenceSupport
from gencode_icedb.general.evidenceDb import EvidenceSource

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


class BaseModel(Model):
    "base for peewee models, used to bind proxy"

    class Meta:
        database = _database_proxy


class TrascriptionSupportLevelField(Field):
    """A field storing a TrascriptionSupportLevel python enumeration as a string in the database."""
    db_field = 'text'

    def db_value(self, value):
        return str(value)

    def python_value(self, value):
        return TrascriptionSupportLevel(value)


class GencodeTranscriptSupport(BaseModel):
    """Transcript-level support"""
    id = PrimaryKeyField()
    transcriptId = CharField(index=True,
                             help_text="""GENCODE trancript id""")
    level = TrascriptionSupportLevelField(index=True,
                                          help_text="""GENCODE TSL""")


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


class GencodeTranscriptSupportDetails(BaseModel):
    """Transcript-level support"""
    id = PrimaryKeyField()
    transcriptId = CharField(index=True,
                             help_text="""GENCODE trancript id""")
    evidSrc = EvidenceSource(index=True,
                             help_text="""Evidence source""")
    evidId = CharField(index=True,
                       help_text="""Evidence identifier""")
    evidSupport = EvidenceSupportField(index=True,
                                       help_text="""Support provided by this piece of evidence""")
