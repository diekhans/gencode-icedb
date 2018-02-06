"""
Types related to support levels.
"""
from pycbio.sys.symEnum import SymEnum

# FIXME: make organization consistent between TSL and RSL stuff


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


class EvidenceSupport(SymEnum):
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

    feat_mismatch = 50  # different number of features
    no_support = 100

    # endsMedium = 2
    # endsWeak = 3
    # suspectMRna = 4
    # estN = 5
    # est1 = 6
    # suspectEst = 7
    # nonUniqueMapping = 8
    # poor = 9
