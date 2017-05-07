"""
Types related to support levels.
"""
from pycbio.sys.symEnum import SymEnum
from collections import namedtuple

class EvidenceType(SymEnum):
    """Type of evidence used in support"""
    __slots__ = ()

    RNA = 1
    EST = 2

class EvidenceComparison(SymEnum):
    """One or more of these attributes describe the support provided to an annotation
    by a given piece of evidence"""
    __slots__ = ()
        good
        endsMedium
        endsWeak
        suspectMRna
        estN
        est1
        suspectEst
        nonUniqueMapping
        poor


