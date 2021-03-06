"""
Code for counting STAR SJ support.
"""
from collections import defaultdict, namedtuple
from gencode_icedb.general.spliceJuncs import SpliceJuncs, spliceJuncsGetSeqs

# FIXME: make naming of supp vs sjsupp consistent
# FIXME including intronMotif in coords means they are different if one source
#       has them and one doesn't.  Bad if a hash key.
# FIXME: switch to using gencodeDb module


class IntronCoords(namedtuple("Intron", ("chrom", "chromStart", "chromEnd", "strand", "intronMotif"))):
    """"Coordinates of an intron that can be used as an hash key"""

    @staticmethod
    def fromSjSupport(sjSupp, intronMotif):
        return IntronCoords(sjSupp.chrom, sjSupp.chromStart, sjSupp.chromEnd, sjSupp.strand, intronMotif)

    def __lt__(self, other):
        "less than, ignoring strand"
        if self.chrom < other.chrom:
            return True
        if self.chromStart < other.chromStart:
            return True
        if self.chromEnd < other.chromEnd:
            return True
        return False

    def __str__(self):
        return "{}:{}-{} ({}) {}".format(self.chrom, self.chromStart, self.chromEnd, self.strand, self.intronMotif)

    @staticmethod
    def fromIntronFeat(intronFeat):
        transFeat = intronFeat.transcript
        return IntronCoords(transFeat.chrom.name, intronFeat.chrom.start, intronFeat.chrom.end, transFeat.rna.strand, intronFeat.sjBases)

    @staticmethod
    def fromIntronSupp(supp):
        return IntronCoords(supp.chrom, supp.chromStart, supp.chromEnd, supp.strand, supp.intronMotif)


class IntronSupportCounts(object):
    """sum of counts for an intron"""
    __slots__ = ("numExprs", "numUniqueMapReads", "numMultiMapReads")

    def __init__(self):
        self.numExprs = self.numUniqueMapReads = self.numMultiMapReads = 0

    def __str__(self):
        return "{} {}".format(self.numUniqueMapReads, self.numMultiMapReads)

    def sum(self, cnts):
        "can sum with SjCounts or SjSupport objects"
        self.numExprs += 1
        self.numUniqueMapReads += cnts.numUniqueMapReads
        self.numMultiMapReads += cnts.numMultiMapReads


class IntronSupportCounter(defaultdict):
    """Collect counts for introns, index by IntronCoords """
    def __init__(self, genomeReader):
        super(IntronSupportCounter, self).__init__(IntronSupportCounts)
        self.genomeReader = genomeReader
        self.genoneMotifCache = {}  # indexed by genome coords

    def _getIntronMotifFromGenome(self, chrom, chromStart, chromEnd, strand):
        key = (chrom, chromStart, chromEnd, strand)
        intronMotif = self.genoneMotifCache.get(key)
        if intronMotif is None:
            intronMotif = spliceJuncsGetSeqs(self.genomeReader, chrom, chromStart, chromEnd, strand)
            self.genoneMotifCache[key] = intronMotif
        return intronMotif

    def _getSjSuppIntronMotif(self, sjSupp):
        "Return motif for start support. if it's not one known by STAR, look up in genome."
        if sjSupp.intronMotif == "??/??":
            return self._getIntronMotifFromGenome(sjSupp.chrom, sjSupp.chromStart, sjSupp.chromEnd, sjSupp.strand)
        else:
            return sjSupp.intronMotif

    def sumSjSupp(self, sjSupp):
        intron = IntronCoords(sjSupp.chrom, sjSupp.chromStart, sjSupp.chromEnd, sjSupp.strand,
                              self._getSjSuppIntronMotif(sjSupp))
        self[intron].sum(sjSupp)  # defaultdict will create

    def __getIntronFeatMotif(self, intronFeat):
        """will be lower-case if not a known slice junction"""
        intronMotif = "{}/{}".format(intronFeat.donorSeq, intronFeat.acceptorSeq)
        if intronFeat.spliceJuncs == SpliceJuncs.unknown:
            intronMotif = intronMotif.lower()
        else:
            intronMotif = intronMotif.upper()
        return intronMotif

    def getIntronFeatCounts(self, intronFeat):
        """given an Intron feature, get the IntronCoords and IntronSupportCounts objects, creating zero counts if needed"""
        intron = IntronCoords(intronFeat.transcript.chrom.name, intronFeat.chrom.start, intronFeat.chrom.end,
                              intronFeat.transcript.rna.strand, self.__getIntronFeatMotif(intronFeat))
        return (intron, self[intron])   # defaultdict will create
