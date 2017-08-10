from pycbio.sys.symEnum import SymEnum
from pycbio.hgdata import dnaOps


class SpliceJuncs(SymEnum):
    "symbolic names for known splice junction patterns"
    unknown = 0
    GT_AG = 1
    CT_AC = 2
    GC_AG = 3
    CT_GC = 4
    AT_AC = 5
    GT_AT = 6


spliceJuncsMap = {
    ("gt", "ag"): SpliceJuncs.GT_AG,
    ("ct", "ac"): SpliceJuncs.CT_AC,
    ("gc", "ag"): SpliceJuncs.GC_AG,
    ("ct", "gc"): SpliceJuncs.CT_GC,
    ("at", "ac"): SpliceJuncs.AT_AC,
    ("gt", "at"): SpliceJuncs.GT_AT
}

# mapping of star splice code, 0 is unknown
starSpliceJuncsMap = {
    1: SpliceJuncs.GT_AG,
    2: SpliceJuncs.CT_AC,
    3: SpliceJuncs.GC_AG,
    4: SpliceJuncs.CT_GC,
    5: SpliceJuncs.AT_AC,
    6: SpliceJuncs.GT_AT
}


class Spliceosome(SymEnum):
    unknown = 0
    major = 1
    minor = 2


spliceosomeMap = {
    SpliceJuncs.GT_AG: Spliceosome.major,
    SpliceJuncs.CT_AC: Spliceosome.minor,
    SpliceJuncs.GC_AG: Spliceosome.minor,
    SpliceJuncs.CT_GC: Spliceosome.minor,
    SpliceJuncs.AT_AC: Spliceosome.minor,
    SpliceJuncs.GT_AT: Spliceosome.minor
}


def spliceJuncsClassify(donorSeq, acceptorSeq):
    return spliceJuncsMap.get((donorSeq.lower(), acceptorSeq.lower()), SpliceJuncs.unknown)


def spliceosomeClassify(spliceJunc):
    assert isinstance(spliceJunc, SpliceJuncs), type(spliceJunc)
    return spliceosomeMap.get(spliceJunc, Spliceosome.unknown)


def spliceJuncsGetSeqs(genomeReader, chrom, intronStart, intronEnd, strand):
    """get the donor and acceptor sequences for the intron, considering
    strand.  If the motif is not a know one, make lower case, otherwise
    upper case if known"""
    donorSeq = genomeReader.get(chrom, intronStart, intronStart + 2)
    acceptorSeq = genomeReader.get(chrom, intronEnd - 2, intronEnd)
    if strand == '-':
        donorSeq, acceptorSeq = dnaOps.reverseComplement(acceptorSeq), dnaOps.reverseComplement(donorSeq)
    if spliceJuncsClassify(acceptorSeq, donorSeq) == SpliceJuncs.unknown:
        return (donorSeq.lower(), acceptorSeq.lower())
    else:
        return (donorSeq.upper(), acceptorSeq.upper())
