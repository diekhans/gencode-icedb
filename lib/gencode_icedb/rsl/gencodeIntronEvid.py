"""
Structure for loading spliceJunctionCollectEvidence along with
gencode annotations for analysis.
"""
from pycbio import tsv
from pycbio.sys.symEnum import SymEnum
from collections import defaultdict


class IntronSupportLevel(SymEnum):
    "support level on an single intron"
    NONE = 0
    WEAK = 1
    MEDIUM = 2
    STRONG = 3


def intronEvidSupportLevel(numReads, isNovel):
    """compute the level of support for a given intron"""
    # easy first pass implementation
    mult = 10 if isNovel else 1
    if numReads >= 1000 * mult:
        return IntronSupportLevel.STRONG
    elif numReads >= 100 * mult:
        return IntronSupportLevel.MEDIUM
    elif numReads >= 1:
        return IntronSupportLevel.WEAK
    else:
        return IntronSupportLevel.NONE


# FIXME: no longer needed

class GencodeIntronEvidReader(tsv.TsvReader):
    "TSV reader for evidence"

    # chrom intronStart intronEnd novel annotStrand rnaSeqStrand intronMotif
    # numUniqueMapReads numMultiMapReads transcripts

    spliceJunctionTsvTypeMap = {
        "intronStart": int,
        "intronEnd": int,
        "novel": lambda v: bool(int(v)),
        "annotStrand": tsv.strOrNoneType,
        "rnaSeqStrand": tsv.strOrNoneType,
        "numUniqueMapReads": int,
        "numMultiMapReads": int,
        "transcripts": (lambda s: s.split(',') if len(s) > 0 else [],
                        lambda l: ",".join(l)),
    }

    def __init__(self, evidTsv):
        super(GencodeIntronEvidReader, self).__init__(evidTsv, typeMap=GencodeIntronEvidReader.spliceJunctionTsvTypeMap)


class GencodeIntronEvidSet(list):
    """
    Set of all intron evidence for a release.
    """
    def __init__(self, evidTsv):
        self.novel = []
        # this mixes transcript loci by PAR, sorted out with get
        self.byTranscriptId = defaultdict(list)

        self.__loadEvid(evidTsv)

    def __loadEvid(self, evidTsv):
        for evid in GencodeIntronEvidReader(evidTsv):
            self.__loadEvidRec(evid)

    def __loadEvidRec(self, evid):
        self.append(evid)
        if len(evid.transcripts) == 0:
            self.novel.append(evid)
        else:
            for ti in evid.transcripts:
                self.byTranscriptId[ti].append(evid)

    def getByTranscribeLocus(self, transLocus):
        "get introns for a transcript locus, sorted in genomic order"
        evids = self.byTranscriptId[transLocus.transcript.id]
        if transLocus.hasMultiLoci():
            # find ones on same chrom
            evids = [evid for evid in evids if evid.chrom == transLocus.chrom]
        return tuple(sorted(evids, key=lambda r: r.intronStart))


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
