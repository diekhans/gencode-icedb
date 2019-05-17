"""
PeeWee data models for RNA-Seq metadata and splice junctions.
"""
import os
from peewee import Proxy, Model, PrimaryKeyField, ForeignKeyField, CharField, TextField, IntegerField
from gencode_icedb.general.peeweeOps import peeweeConnect, peeweeClose, peeweeClassToTableName, PeeweeModelMixins
from gencode_icedb.general.dbModels import EvidenceSource, AnalysisStatusField
from collections import namedtuple
import pysam

_database_proxy = Proxy()


def setDatabaseConn(dbconn):
    "bind the proxy to a database"
    _database_proxy.initialize(dbconn)


def rslConnect(dburl, create=False, readonly=True, timeout=None, synchronous=None):
    "connect to sqlite3 database and bind to model"
    return peeweeConnect(dbUrl, setDatabaseConn, create=create, readonly=readonly, timeout=timeout, synchronous=synchronous)


def rslClose(conn):
    "close database"
    peeweeClose(conn)
    _database_proxy = Proxy()

class BaseModel(Model, PeeweeModelMixins):
    "base for peewee models, used to bind proxy"

    class Meta:
        database = _database_proxy
        table_function = peeweeClassToTableName


class RslEvidenceSource(EvidenceSource, BaseModel):
    """Sources of evidence associated with RNA-Seq introns."""
    # derived class creates a unique table
    pass


class RslAnalysis(BaseModel):
    """Result of running analysis on RslEvidenceSource
    """
    id = PrimaryKeyField()
    create_time = DateTimeField(help_text="""Date/time registered in database""")
    update_time = DateTimeField(help_text="""Date/time of last update in database""")
    rsl_evidence_source__id = ForeignKeyField(RslEvidenceSource,
                                              help_text="""Associated evidence source""")
    assembly = CharField(help_text="""genome assemble used""")
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
    # FIXME better name the mapping_symid
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
        self.tabFile = tabFile
        if self.tabFile is None:
            if sjDbPath is None:
                sjDbPath = sjDbConn.database
            self.tabFile = self.sjTabFromSjDb(sjDbPath)
        self.tb = pysam.TabixFile(self.tabFile)
        self.chroms = frozenset(self.tb.contigs)

    def close(self):
        if self.tb is not None:
            try:
                self.tb.close()
            finally:
                self.tb = None

    def fetch(self, chrom, start, end):
        """Query returning SjSupport.  Use zero-based, half open coordinates.
        in chrom is not in index, nothing is returned"""
        if chrom in self.chroms:
            for line in self.tb.fetch(chrom, start, end):
                yield SjSupport.factory(line.split("\t"))


class GencodeIntronSupport(BaseModel):
    """Results from comparing support to GENCODE.  De-normalized and not
    linked for now"""
    # FIXME: not so sure de-normalized is best, need to keep checking for
    # duplication in code and have record that are sometimes GencodeIntronSupport
    # and sometimes GencodeIntronNovel
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


class GencodeIntronNovel(BaseModel):
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
