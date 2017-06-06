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
# FIXME: base classification on splicing mechanism, not base-pairs
SpliceSites = SymEnum("SpliceSite",
                      ("GT_AG", "GC_AG", "AT_AC", "unknown"))
spliceSitesMap = {
    ("gt", "ag"): SpliceSites.GT_AG,
    ("gc", "ag"): SpliceSites.GC_AG,
    ("at", "ac"): SpliceSites.AT_AC,
}

Spliceosome = SymEnum("Spliceosome",
                      ("major", "minor", "unknown"))

spliceosomeMap = {
    SpliceSites.GT_AG: Spliceosome.major,
    SpliceSites.GC_AG: Spliceosome.minor,
    SpliceSites.AT_AC: Spliceosome.minor,
}


def spliceSitesClassify(donor, acceptor):
    return spliceSitesMap.get((donor.lower(), acceptor.lower()), SpliceSites.unknown)


def spliceSitesClassifyStrand(strand, donorSeq, acceptorSeq):
    if strand == '+':
        return spliceSitesClassify(donorSeq, acceptorSeq)
    else:
        return spliceSitesClassify(dnaOps.reverseComplement(acceptorSeq), dnaOps.reverseComplement(donorSeq))


def spliceosomeClassify(spliceSites):
    assert isinstance(spliceSites, SpliceSites)
    return spliceosomeMap.get(spliceSites, Spliceosome.unknown)
