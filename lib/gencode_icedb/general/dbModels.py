"""
Common peewee base model classes.  These define the common types from which the specific
table types are derived.
"""
from peewee import Proxy, Model, AutoField, ForeignKeyField, CharField, TextField, IntegerField, DateTimeField, UUIDField
from gencode_icedb.general.peeweeOps import peeweeConnect, peeweeClose, peeweeClassToTableName, PeeweeModelMixins, SymEnumField
from pycbio.sys.symEnum import SymEnum

class EvidenceType(SymEnum):
    """Type of evidence"""
    __slots__ = ()
    GENBANK_RNA = 0    # GenBank RNAs
    GENBANK_EST = 1    # GenBank ESTs


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
    id = AutoField()
    create_time = DateTimeField(help_text="""Date and time row created""")
    update_time = DateTimeField(help_text="""Date and time row was last updated""")
    src_name = CharField(index=True,
                         help_text="""symbolic name of source (e.g. SRA)""")
    acc = CharField(unique=True,
                    help_text="""accession used by database""")
    src_time = DateTimeField(help_text="""Date/time available at source""")
    org_code = CharField(index=True, max_length=2,
                         help_text="""two-character organism code: 'hs' or 'mm'""")
    assembly = CharField(index=True, null=True,
                         help_text="""Assembly name, if aligned to an assembly, otherwise None""")
    tissue = CharField(null=True, default=None,
                       help_text="""Tissue or body site, if available.  This is not normalized and maybe hard to interpret.""")
    tumor = IntegerField(null=True, default=None,
                         help_text="""1 for tumor, 0 for normal, null if unknown""")
    seq_platform = CharField(help_text="""Sequencing platform""")
    data_format = CharField(help_text="""Data format: "SRA", FASTQ, BAM, or CRAM""")
    url1 = CharField(help_text="""URL of first file""")
    url2 = CharField(null=True, default=None,
                     help_text="""URL of second file, for pair-end FASTQ""")
    comments = TextField(help_text="""comments""")

class IntronEvidence(Model):
    """Collection of intron evidence from a source.  Both short and long reads maybe used."""
    id = AutoField()
    create_time = DateTimeField(help_text="""Date and time row created""")
    update_time = DateTimeField(help_text="""Date and time row was last updated""")
    source = ForeignKeyField(EvidenceSource, on_delete="CASCADE",
                             help_text="""evidence being analyzed""")
    mapping_algo = CharField(null=True, default=None,
                             help_text="""mapping technique used (method and parameters), NULL if already aligned.""")
    calling_algo = CharField(help_text="""intron calling technique (method and parameters)""")

    status = SymEnumField(AnalysisStatus,
                          help_text="""status of the analysis""")
