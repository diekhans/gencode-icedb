"""
Types related to support levels.
"""
from pycbio.sys.symEnum import SymEnum

# FIXME: make organization consistent between TSL and RSL stuff

UCSC_RNA_ALN_TBL = "ucsc_rna_aln"
UCSC_EST_ALN_TBL = "ucsc_est_aln"
ENSEMBL_RNA_ALN_TBL = "ensembl_rna_aln"

GENCODE_ANN_TBL = "gencode_ann"
GENCODE_ATTRS_TBL = "gencode_attrs"
GENCODE_TAG_TBL = "gencode_tag                        "
GENCODE_TRANSCRIPT_SOURCE_TBL = "gencode_transcript_source          "
GENCODE_TRANSCRIPTION_SUPPORT_LEVEL_TBL = "gencode_transcription_support_level"


class TrascriptionSupportLevel(SymEnum):
    """Transcription support levels"""
    tslNA = -1
    tsl1 = 1
    tsl2 = 2
    tsl3 = 3
    tsl4 = 4
    tsl5 = 5


class EvidenceType(SymEnum):
    """Type of evidence used in support"""
    __slots__ = ()

    RNA = 1
    EST = 2


class EvidenceComparison(SymEnum):
    """One or more of these attributes describe the support provided to an annotation
    by a given piece of evidence"""
    __slots__ = ()
    good = 1
    polymorphic = 2

    poor = 10  # value greater than this a detail of poor
    exon_boundry_mismatch = 11
    large_indel_size = 12  # give indel exceeds a threshold
    large_indel_content = 14  # indel context if a given exon exceeded
    internal_unaligned = 15  # intron contains unaligned
    feat_mismatch = 16  # different number of features

    # endsMedium = 2
    # endsWeak = 3
    # suspectMRna = 4
    # estN = 5
    # est1 = 6
    # suspectEst = 7
    # nonUniqueMapping = 8
    # poor = 9
