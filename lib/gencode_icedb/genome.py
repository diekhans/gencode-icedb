"""
Functions around reading and classifying sequence
"""
from __future__ import print_function
from pycbio.hgdata import dnaOps
from pycbio.sys.symEnum import SymEnum


class GenomeReader(object):
    """
    Reads sequences from twoBitReader, which can be substituted with a mock
    version.
    """
    def __init__(self, twoBitReader):
        self.twoBitReader = twoBitReader
        self.chrom = self.size = None

    def __obtainChrom(self, chrom):
        if self.chrom != chrom:
            self.twoBitSeq = self.twoBitReader[chrom]
            self.chrom = chrom
            self.size = len(self.twoBitSeq)

    def get(self, chrom, start, end, strand=None):
        self.__obtainChrom(chrom)
        if strand == '-':
            start, end = dnaOps.reverseCoords(start, end, self.size)
        seq = self.twoBitSeq[start:end]
        if strand == '-':
            seq = dnaOps.reverseComplement(seq)
        return seq

    def getChromSize(self, chrom):
        self.__obtainChrom(chrom)
        return self.size


    # FIXME: move to a different module
class SpliceJunc(SymEnum):
    "symbolic names for known splice junction patterns"
    unknown = 0
    GT_AG   = 1
    CT_AC   = 2
    GC_AG   = 3
    CT_GC   = 4
    AT_AC   = 5
    GT_AT   = 6

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
    major   = 1
    minor   = 2

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
