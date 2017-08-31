"""
Types related to support levels.
"""
from pycbio.sys.symEnum import SymEnum


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
    endsMedium = 2
    endsWeak = 3
    suspectMRna = 4
    estN = 5
    est1 = 6
    suspectEst = 7
    nonUniqueMapping = 8
    poor = 9
