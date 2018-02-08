"""
library functions and class for classify tests
"""
from collections import defaultdict
from pycbio.tsv import TsvReader
from pycbio.hgdata.coords import Coords


class TestGene(list):
    "Has rows from TSV for each transcript in a gene, by transcriptId"
    def __init__(self):
        self.chrom = self.start = self.end = None

    def add(self, trans):
        self.append(trans)
        if self.chrom is None:
            self.chrom = trans.chrom
            self.start = trans.txStart
            self.end = trans.txEnd
        else:
            self.start = min(self.start, trans.txStart)
            self.end = max(self.end, trans.txEnd)


class TestCases(defaultdict):
    """TestGene gene id from TSV
    """
    def __init__(self, testCaseTsv):
        super().__init__(TestGene)
        self.byTranscriptId = dict()
        for trans in TsvReader(testCaseTsv,
                               typeMap={"txStart": int,
                                        "txEnd": int,
                                        "level": int}):
            self[trans.geneId].add(trans)
            self.byTranscriptId[trans.transcriptId] = trans

    def getRanges(self):
        return [Coords(g.chrom, g.start, g.end) for g in self.values()]
