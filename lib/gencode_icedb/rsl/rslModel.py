"""
PeeWee data models for RNA-Seq metadata and splice junctions.
"""
from __future__ import print_function
import os
from peewee import Proxy, Model, PrimaryKeyField, ForeignKeyField, CharField, TextField, IntegerField
from playhouse.apsw_ext import APSWDatabase
import apsw
from collections import namedtuple
import pysam

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


def sqliteClose(rslDb):
    "close database if not-null and open"
    # not sure why it might be in closed state even after open, maybe lazy open?
    if not rslDb.is_closed():
        rslDb.close()


def sqliteSetSynchronous(dbconn, mode):
    if mode is False:
        mode = "OFF"
    elif mode is True:
        mode = "NORMAL"
    dbconn.execute_sql("PRAGMA synchronous={}".format(mode))


class BaseModel(Model):
    "base for peewee models, used to bind proxy"

    class Meta:
        database = _database_proxy


class RunMetadata(BaseModel):
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


class MappingParameters(BaseModel):
    """Mapping parameters, stored into a separate table due to use of same parameters for many runs"""
    id = PrimaryKeyField()
    mapping_param_symid = CharField(unique=True,
                                    help_text="""symbolic name of the parameter set""")
    assembly = CharField(help_text="""genome assemble""")
    gene_set = CharField(help_text="""gene annotations used in mapping""")
    commands = TextField(help_text="""commands used to do mapping""")
    comments = TextField(help_text="""comments""")


class MappingMetadata(BaseModel):
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


class SjSupport(namedtuple("SjSupport", ("chrom", "chromStart", "chromEnd",
                                         "strand", "intronMotif", "annotated",
                                         "numUniqueMapReads", "numMultiMapReads",
                                         "maxOverhang", "mapping_symid"))):
    """SjSupport record created from STAR sjout files.  These are loaded from
    a tabix indexed tab file."""
    __slots__ = ()

    def __str__(self):
        return "{}:{}-{} ({}) {}: annot={} uniq={} multi={} over={} expr={}".format(self.chrom, self.chromStart, self.chromEnd, self.strand,
                                                                                    self.intronMotif, self.annotated, self.numUniqueMapReads,
                                                                                    self.numMultiMapReads, self.maxOverhang, self.mapping_symid)

    @property
    def start(self):
        # FIXME make up our mind on names
        return self.chromStart

    def end(self):
        # FIXME make up our mind on names
        return self.chromEnd

    @staticmethod
    def factory(row):
        return SjSupport(row[0], int(row[1]), int(row[2]),
                         row[3], row[4], bool(int(row[5])),
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
        self.tb = pysam.TabixFile(tabFile)
        self.chroms = frozenset(self.tb.contigs)

    def close(self):
        if self.tb is not None:
            try:
                self.tb.close()
            finally:
                self.tb = None

    def __del__(self):
        self.close()

    def fetch(self, chrom, start, end):
        """Query returning SjSupport.  Use zero-based, half open coordinates.
        in chrom is not in index, nothing is returned"""
        if chrom in self.chroms:
            for line in self.tb.fetch(chrom, start, end):
                yield SjSupport.factory(line.split("\t"))


class GencodeSupport(BaseModel):
    """Results from comparing support to GENCODE.  De-normalized and not
    linked for now"""
    # FIXME: not so sure de-normalized is best, need to keep checking for
    # duplication in code and have record that are sometimes GencodeSupport
    # and sometimes GencodeNovel
    id = PrimaryKeyField()
    geneId = CharField(index=True,
                       help_text="""GENCODE gene id""")
    geneName = CharField(index=True,
                         help_text="""gene name/symbol""")
    transcriptId = CharField(index=True,
                             help_text="""GENCODE trancript id""")
    transcriptType = CharField(index=True,
                               help_text="""GENCODE transcript type""")
    chrom = CharField(help_text="""Chromosome""")
    intronStart = IntegerField(help_text="""zero-based start of intron""")
    intronEnd = IntegerField(help_text="""end of intron""")
    strand = CharField(help_text="""strand of gene""")
    intronMotif = CharField(index=True,
                            help_text="""Intron splice junction motif in the forms AT/GC.  If splice site is not a known """
                            """motif, the motif is in lower case. """)
    numExprs = IntegerField(help_text="""total number of experiments having this intron""")
    numUniqueMapReads = IntegerField(index=True,
                                     help_text="""total number of uniquely mapping reads""")
    numMultiMapReads = IntegerField(index=True,
                                    help_text="""total number of multi-mapping reads""")


class GencodeNovel(BaseModel):
    """Support for novel intron from comparing support to GENCODE."""
    id = PrimaryKeyField()
    chrom = CharField(help_text="""Chromosome""")
    intronStart = IntegerField(help_text="""zero-based start of intron""")
    intronEnd = IntegerField(help_text="""end of intron""")
    strand = CharField(help_text="""strand of gene""")
    intronMotif = CharField(index=True,
                            help_text="""Intron splice junction motif in the forms AT/GC.  If splice site is not a known """
                            """motif, the motif is in lower case. """)
    numExprs = IntegerField(help_text="""total number of experiments having this intron""")
    numUniqueMapReads = IntegerField(index=True,
                                     help_text="""total number of uniquely mapping reads""")
    numMultiMapReads = IntegerField(index=True,
                                    help_text="""total number of multi-mapping reads""")
    geneIds = CharField(index=True,
                        help_text="""GENCODE gene ids if intron overlaps a gene""")
