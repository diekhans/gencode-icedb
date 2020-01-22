"""
Common peewee base model classes.  These define the common types from which the specific
table types are derived.
"""
import os
import uuid
import datetime
from peewee import Model, AutoField, ForeignKeyField, CharField, TextField, IntegerField, UUIDField
# see https://github.com/coleifer/peewee/issues/771
from playhouse.apsw_ext import DateTimeField
from gencode_icedb.general.peeweeOps import SymEnumField
from pycbio.sys.symEnum import SymEnum, auto

# Notes:
#  - UUID fields
#    GUID key fields, implemented as UUIDs are added to make it easy to move data
#    between local an InnoDb database and Sqlite3 databases for cluster access.
#    However, some sites indicate using a UUID as a primary or unique key has serious
#    performance implications:
#       http://kccoder.com/mysql/uuid-vs-int-insert-performance/
#    Thus we keep both table-local auto-increment primary keys and UUIDs for each record.
#    Although perhaps some solutions worth looking at:
#       http://mysql.rjweb.org/doc.php/uuid

def dbGuidGen():
    """create a GUID"""
    return uuid.uuid4()


def dbCurrentTime():
    """Current GMT for setting database fields"""
    return datetime.datetime.utcnow()


class EvidenceRepo(SymEnum):
    """Repository source of evidence"""
    __slots__ = ()
    GENBANK = auto()
    ARRAY_EXPRESS = auto()
    SRA = auto()


class EvidenceType(SymEnum):
    """Type of evidence"""
    __slots__ = ()
    GENBANK_RNA = 0    # GenBank RNAs
    GENBANK_EST = 1    # GenBank ESTs
    SHORT_RNASEQ = 2   # Short-read RNA-Seq
    LONG_RNASEQ = 3   # Short-read RNA-Seq


class RnaSeqFileFormat(SymEnum):
    """Format of RNA sequence file"""
    __slots__ = ()
    BAM = auto()
    CRAM = auto()

    @classmethod
    def fromFileName(cls, fname):
        ext = os.path.splitext(fname)
        if ext == ".bam":
            return cls.BAM
        elif ext == ".cram":
            return cls.CRAM
        else:
            raise ValueError("don't know how to convert to a RnaSeqFileFormat: {}".format(fname))


class AnalysisStatus(SymEnum):
    """Status for running an analysis"""
    __slots__ = ()
    NONE = 0             # not assigned
    ACCESS_FAILED = 1    # failed to obtain input data
    ANALYSIS_FAILED = 2  # analysis process failure
    SUCCESS = 3          # analysis successful


class EvidenceSource(Model):
    """Metadata associated with an source of evidence.  This is obtained from
    source database (e.g. SRA).  It maybe unaligned or aligned reads.  Entry
    in this table means data is available, although may not have been
    processed.
    """
    id = AutoField(help_text="""Database table unique id of source""")
    guid = UUIDField(index=True,
                     help_text="""Global unique id assigned to source, allows moving data between databases""")
    create_time = DateTimeField(help_text="""Date and time row created""")
    update_time = DateTimeField(help_text="""Date and time row was last updated""")
    repo = SymEnumField(symEnumCls=EvidenceRepo,
                        help_text="""repository source (e.g. SRA)""")
    sample_id = CharField(index=True, null=True,
                          help_text="""sample id""")
    run_id = CharField(unique=True, index=True,
                       help_text="""run id""")
    src_id = CharField(index=True,
                       help_text="""id used by source database""")
    src_time = DateTimeField(help_text="""Date/time available at source""")
    organism = CharField(index=True,
                         help_text="""organism name in Ensembl format (e.g.  mus_musculus)""")
    assembly = CharField(index=True, null=True,
                         help_text="""Assembly name, if aligned to an assembly, otherwise None""")
    tissue = CharField(null=True, default=None,
                       help_text="""Tissue or body site, if available.  This is not normalized and maybe hard to interpret.""")
    disease = IntegerField(null=True, default=None,
                           help_text="""1 for diseased, usually tumor, 0 for normal, null if unknown""")
    evidence_type = SymEnumField(symEnumCls=EvidenceType,
                                 help_text="""Evidence type category""")
    platform = CharField(help_text="""Sequencing platform""")
    data_format = SymEnumField(symEnumCls=RnaSeqFileFormat,
                               help_text="""Sequence or alignment file format""")
    url1 = CharField(help_text="""URL of first file""")
    url2 = CharField(null=True, default=None,
                     help_text="""URL of second file, for pair-end FASTQ""")
    comments = TextField(help_text="""comments""", null=True)


class EvidenceAnalysis(Model):
    """Analysis of an evidence file. """
    id = AutoField(help_text="""Database table unique id of source""")
    guid = UUIDField(index=True,
                     help_text="""Global unique id assigned to source, allows moving data between databases""")
    create_time = DateTimeField(help_text="""Date and time row created""")
    update_time = DateTimeField(help_text="""Date and time row was last updated""")
    source_id = ForeignKeyField(EvidenceSource, on_delete="CASCADE",
                                help_text="""evidence being analyzed""")
    source_guid = UUIDField(help_text="""GUID of evidence being analyzed""")
    assembly = CharField(help_text="""genome assemble used""")
    status = SymEnumField(AnalysisStatus,
                          help_text="""status of the analysis""")
