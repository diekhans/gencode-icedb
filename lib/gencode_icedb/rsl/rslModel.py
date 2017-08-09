"""
PeeWee data models for RNA-Seq metadata and splice junctions.
"""
from __future__ import print_function
import os
from peewee import Proxy, Model, PrimaryKeyField, ForeignKeyField, CharField, TextField
from playhouse.apsw_ext import APSWDatabase
import apsw
from collections import namedtuple
import tabix

_database_proxy = Proxy()


def setDatabaseConn(dbconn):
    "bind the proxy to a database"
    _database_proxy.initialize(dbconn)


def sqliteConnect(rsldb, readonly=False, timeout=None, synchronous=None):
    "connect to sqlite3 database and bind to model"
    kwargs = {}
    if timeout is not None:
        kwargs["timeout"] = timeout
    if readonly:
        kwargs["flags"] = apsw.SQLITE_OPEN_READONLY
    dbconn = APSWDatabase(rsldb, **kwargs)
    setDatabaseConn(dbconn)
    if synchronous is not None:
        sqliteSetSynchronous(dbconn, synchronous)
    return dbconn


def sqliteSetSynchronous(dbconn, mode):
    if mode is False:
        mode = "OFF"
    elif mode is True:
        mode = "NORMAL"
    dbconn.execute_sql("PRAGMA synchronous={}".format(mode))


class RunMetadata(Model):
    """Metadata associated with an RNA-Seq sequencing run.  This
    is obtained from source database (e.g. SRA)
    """
    id = PrimaryKeyField()
    src_symid = CharField(index=True,
                          help_text="""name of source (e.g. SRA)""")
    run_acc = CharField(unique=True,
                        help_text="""sequencing run accession""")
    org_code = CharField(index=True, max_length=2,
                         help_text="""two-character organism code: 'hs' or 'mm'""")
    tumor = CharField(null=True, default=None,
                      help_text="""yes for tumor, no for normal, null if unknown """)
    tissue = CharField(null=True, default=None,
                       help_text="""Tissue or body site, if available.  This is not normalized and maybe hard to interpret.""")

    class Meta:
        database = _database_proxy


class MappingParameters(Model):
    """Mapping parameters, stored into a separate table due to use of same parameters for many runs"""
    id = PrimaryKeyField()
    mapping_param_symid = CharField(unique=True,
                                    help_text="""symbolic name of the parameter set""")
    assembly = CharField(help_text="""genome assemble""")
    gene_set = CharField(help_text="""gene annotations used in mapping""")
    commands = TextField(help_text="""commands used to do mapping""")
    comments = TextField(help_text="""comments""")

    class Meta:
        database = _database_proxy


class MappingMetadata(Model):
    """Metadata associated with of an RNA-Seq mapping and splice junction STAR run.
    """
    id = PrimaryKeyField()
    run_metadata_id = ForeignKeyField(RunMetadata,
                                      help_text="""sequencing run metadata""")
    mapping_symid = CharField(unique=True,
                              help_text="""symbolic name of the mapping, this is defined locally""")
    mapping_acc = CharField(unique=True, null=True, default=None,
                            help_text="""mapping analysis accession, only available if results have been submitted to an archive""")
    mapping_parameters_id = ForeignKeyField(MappingParameters,
                                            help_text="""parameters used in mapping""")

    class Meta:
        database = _database_proxy


class SjSupport(namedtuple("SjSupport", ("chrom", "chromStart", "chromEnd",
                                         "strand", "intronMotif", "annotated",
                                         "numUniqueMapReads", "numMultiMapReads",
                                         "maxOverhang", "mapping_symid"))):
    """SjSupport record created from STAR sjout files.  These are loaded from
    a tabix indexed tab file."""
    @staticmethod
    def factory(row):
        return SjSupport(row[0], int(row[1]), int(row[2]),
                         row[3], row[4], bool(row[5]),
                         int(row[6]), int(row[7]), int(row[8]), row[9])


class SjSupportReader(object):
    """reader for SjSupport from a tabix-indexed file"""
    @staticmethod
    def sjTabFromSjDb(sjDbPath):
        "deduce tab file name from sjDbPath"
        return os.path.splitext(sjDbPath)[0] + ".sjsup.gz"

    def __init__(self, tabFile=None, sjDbConn=None, sjDbPath=None):
        """can open by tab file name, a APSWDatabase object, or path to that
        the database files"""
        self.tb = None
        if tabFile is None:
            if sjDbPath is None:
                sjDbPath = sjDbConn.database
            tabFile = self.sjTabFromSjDb(sjDbPath)
        self.tb = tabix.open(tabFile)

    def close(self):
        if self.tb is not None:
            del self.tb  # doesn't have a close method
            self.tb = None

    def __del__(self):
        if self.tb is not None:
            self.close()

    def query(self, chrom, start, end):
        "Query returning SjSupport.  Use zero-based, half open coordinates"
        for row in self.tb.query(chrom, start, end):
            yield SjSupport.factory(row)
