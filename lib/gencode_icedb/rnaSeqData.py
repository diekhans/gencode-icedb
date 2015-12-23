"""
Data associate with RNA seq database, store in an sqlite database
"""
from peewee import *

database_proxy = Proxy()

def setDatabase(database):
    "bind the proxy to a database"
    database_proxy.initialize(database)

class RnaSeqData(Model):
    """specification of a single set of RNA-Seq reads.
    name is a symbolic, unique and is used to create output directories.  Reads can be
    SAM, BAM, or fastq and can be compressed.  readsFile path is relative
    to data root."""
    id = IntegerField(primary_key=True)
    setName = CharField()
    runName = CharField(unique=True)
    created = DateField()
    organism = CharField()
    readsfileurl = CharField(unique=True)
    readsfile = CharField()
    readsfile2url = CharField(unique=True, null=True)
    readsfile2 = CharField(null=True)
    readlength = IntegerField()
    description = CharField()
    tissue = CharField()

    class Meta:
        database = database_proxy

