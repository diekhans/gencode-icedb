"""
Structure for loading spliceJunctionCollectEvidence along with
gencode annotations for analysis.
"""
from pycbio.sys.symEnum import SymEnum

# FIXME: name confusion with supportAnalysis.GencodeIntronEvid


class IntronSupportLevel(SymEnum):
    "support level on an single intron"
    NONE = 0
    WEAK = 1
    MEDIUM = 2
    STRONG = 3


def intronEvidSupportLevel(numReads):
    """compute the level of support for a given intron"""
    # easy first pass implementation
    if numReads >= 1000:
        return IntronSupportLevel.STRONG
    elif numReads >= 100:
        return IntronSupportLevel.MEDIUM
    elif numReads >= 1:
        return IntronSupportLevel.WEAK
    else:
        return IntronSupportLevel.NONE


class SpliceJuncCat(SymEnum):
    consensus = 1
    known = 2
    unknown = 3

    @staticmethod
    def fromJunc(junc):
        "categorize in the form"
        if junc == "GT/AG":
            return SpliceJuncCat.consensus
        elif junc in ("CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT"):
            return SpliceJuncCat.known
        else:
            return SpliceJuncCat.unknown
