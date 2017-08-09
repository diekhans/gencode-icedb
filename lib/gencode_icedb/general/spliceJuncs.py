from pycbio.sys.symEnum import SymEnum
from pycbio.hgdata import dnaOps


class SpliceJunc(SymEnum):
    "symbolic names for known splice junction patterns"
    unknown = 0
    GT_AG = 1
    CT_AC = 2
    GC_AG = 3
    CT_GC = 4
    AT_AC = 5
    GT_AT = 6


spliceJuncMap = {
    ("gt", "ag"): SpliceJunc.GT_AG,
    ("ct", "ac"): SpliceJunc.CT_AC,
    ("gc", "ag"): SpliceJunc.GC_AG,
    ("ct", "gc"): SpliceJunc.CT_GC,
    ("at", "ac"): SpliceJunc.AT_AC,
    ("gt", "at"): SpliceJunc.GT_AT
}

# mapping of star splice code, 0 is unknown
starSpliceJuncMap = {
    1: SpliceJunc.GT_AG,
    2: SpliceJunc.CT_AC,
    3: SpliceJunc.GC_AG,
    4: SpliceJunc.CT_GC,
    5: SpliceJunc.AT_AC,
    6: SpliceJunc.GT_AT
}


class Spliceosome(SymEnum):
    unknown = 0
    major = 1
    minor = 2


spliceosomeMap = {
    SpliceJunc.GT_AG: Spliceosome.major,
    SpliceJunc.CT_AC: Spliceosome.minor,
    SpliceJunc.GC_AG: Spliceosome.minor,
    SpliceJunc.CT_GC: Spliceosome.minor,
    SpliceJunc.AT_AC: Spliceosome.minor,
    SpliceJunc.GT_AT: Spliceosome.minor
}


def spliceJuncClassify(donor, acceptor):
    return spliceJuncMap.get((donor.lower(), acceptor.lower()), SpliceJunc.unknown)


def spliceJuncClassifyStrand(strand, donorSeq, acceptorSeq):
    if strand == '+':
        return spliceJuncClassify(donorSeq, acceptorSeq)
    else:
        return spliceJuncClassify(dnaOps.reverseComplement(acceptorSeq), dnaOps.reverseComplement(donorSeq))


def spliceosomeClassify(spliceJunc):
    assert isinstance(spliceJunc, SpliceJunc), type(spliceJunc)
    return spliceosomeMap.get(spliceJunc, Spliceosome.unknown)
