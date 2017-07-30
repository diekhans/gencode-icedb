"""
PeeWee data models for RNA-Seq metadata and splice junctions.
"""
from peewee import Proxy, Model, ForeignKeyField, IntegerField, CharField, TextField
from peewee import SqliteDatabase

_database_proxy = Proxy()


def setDatabaseConn(dbconn):
    "bind the proxy to a database"
    _database_proxy.initialize(dbconn)


def sqliteConnect(rsldb):
    "connect to sqlite3 database and bind to model"
    dbconn = SqliteDatabase(rsldb)
    setDatabaseConn(dbconn)
    return dbconn


class RunMetadata(Model):
    """Metadata associated with an RNA-Seq sequencing run.  This
    is obtained from source database (e.g. SRA)
    """
    id = IntegerField(primary_key=True)
    srcname = CharField(index=True,
                        help_text="""name of source""")
    acc = CharField(unique=True,
                    help_text="""sequencing run accession""")
    orgcode = CharField(index=True, max_length=2,
                        help_text="""two-character organism code: 'hs' or 'mm'""")
    tissue = CharField(null=True,
                       help_text="""Tissue or body site, if available.  This is not normalized and maybe hard to interpret.""")

    class Meta:
        database = _database_proxy


class MappingParameters(Model):
    """Mapping parameters, stored into a separate table due to use of same parameters for many runs"""
    id = IntegerField(primary_key=True)
    name = CharField(unique=True,
                     help_text="""name of the parameter set""")
    commands = TextField(help_text="""commands used to do mapping""")
    annotations = CharField(help_text="""gene annotations used in mapping""")

    class Meta:
        database = _database_proxy


class MappingMetadata(Model):
    """Metadata associated with of an RNA-Seq mapping and splice junction STAR run.
    """
    id = IntegerField(primary_key=True)
    run_metadata_id = ForeignKeyField(RunMetadata,
                                      help_text="""sequencing run metadata""")
    name = CharField(unique=True,
                     help_text="""name of the mapping, this is defined locally""")
    acc = CharField(unique=True, null=True,
                    help_text="""mapping analysis accession, only available if results have been submitted to an archive""")
    mapping_parameters_id = ForeignKeyField(MappingParameters,
                                            help_text="""parameters used in mapping""")


class PutativeIntron(Model):
    """Location of a putative intron"""
    id = IntegerField(primary_key=True)
    bin = IntegerField(index=True,
                       help_text="""spacial indexing bin""")
    chrom = CharField(index=True,
                      help_text="""name of the mapping, this is defined locally""")
    start = IntegerField(index=True,
                         help_text="""zero-based, half-open start of intron""")
    end = IntegerField(help_text="""zero-based, half-open end of intron""")
    strand = CharField(help_text="""apparent strand (`. if undefined)""")
    motif = CharField(help_text="""motif for intron (??/?? for non-canonical)""")

    class Meta:
        database = _database_proxy


class SpliceJuncSupport(Model):
    """support for a splice junction from a STAR run"""
    id = IntegerField(primary_key=True)
    mapping_metadata_id = ForeignKeyField(MappingMetadata,
                                          help_text="""metadata describing mapping""")
    putative_intron_id = ForeignKeyField(PutativeIntron,
                                         help_text="""putative intron these counts are associated with""")
    annotated = IntegerField(help_text="""0: unannotated, 1: annotated""")
    num_uniq_reads = IntegerField(help_text="""number of uniquely mapping reads crossing the junction""")
    num_multi_reads = IntegerField(help_text="""number of multi-mapping reads crossing the junction""")
    max_overhang = IntegerField(help_text="""maximum spliced alignment overhang""")

    class Meta:
        database = _database_proxy
