"""
Various functions related to STAR short read mapper.
"""
from collections import namedtuple


class StarSjOutRow(namedtuple("StarSjOutRow",
                              ("chrom", "start", "end", "strand", "motif", "annotated",
                               "num_uniq_reads", "num_multi_reads", "max_overhang"))):
    """Row from a STAR sjout file, with coordinates as-is (1-based, closed)"""
    __slots__ = ()

    @staticmethod
    def factory(row):
        "construct named tuple from row"
        return StarSjOutRow(row[0], int(row[1]), int(row[2]), int(row[3]), int(row[4]), int(row[5]), int(row[6]), int(row[7]), int(row[8]))


_starStrandCodeToChar = {0: '.', 1: '+', 2: '-'}


def starStrandCodeToChar(starStrandCode):
    """convert STAR strand code to a character"""
    return _starStrandCodeToChar[starStrandCode]


_starMotifCodeToStr = {0: "??/??", 1: "GT/AG", 2: "CT/AC", 3: "GC/AG", 4: "CT/GC", 5: "AT/AC", 6: "GT/AT"}


def starMotifCodeToStr(starMotifCode):
    """convert STAR motif code to a string"""
    return _starMotifCodeToStr[starMotifCode]
