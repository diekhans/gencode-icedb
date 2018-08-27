"""
Common code used to build a TranscriptFeature object for an annotation.
"""
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, TranscriptFeatures, Utr5RegionFeature, CdsRegionFeature, Utr3RegionFeature, GapAnnotFeature, NonCodingRegionFeature
from gencode_icedb.general.spliceJuncs import spliceJuncsGetSeqs


class AnnotFeatureBuilder(object):
    def __init__(self, chrom, rna, cdsChrom, transcriptionStrand, attrs, genomeReader):
        self.genomeReader = genomeReader
        self.transAnnot = TranscriptFeatures(chrom, rna, transcriptionStrand=transcriptionStrand, cdsChrom=cdsChrom, attrs=attrs)
        self.transAnnot.features = []
        self.exon = None

    def finish(self):
        "finish off transcript"
        self.transAnnot.features = tuple(self.transAnnot.features)
        assert sum([len(e.rna) for e in self.transAnnot.getFeaturesOfType(ExonFeature)]) == self.transAnnot.rna.size
        return self.transAnnot

    @property
    def hasCds(self):
        return self.transAnnot.cdsChrom is not None

    def beginExon(self, chromStart, chromEnd, rnaStart, rnaEnd):
        """start an new exon"""
        assert self.exon is None
        self.exon = ExonFeature(self.transAnnot, len(self.transAnnot.features),
                                self.transAnnot.chrom.subrange(chromStart, chromEnd),
                                self.transAnnot.rna.subrange(rnaStart, rnaEnd))
        self.exon.annotFeatures = []
        self.transAnnot.features.append(self.exon)

    def finishExon(self):
        """finish up current exon"""
        self.exon.annotFeatures = tuple(self.exon.annotFeatures)
        self.exon = None

    def addGap(self, chromStart, chromEnd, rnaNext):
        """add a gap annotation to exon, returning next RNA position"""
        gap = GapAnnotFeature(self.exon, len(self.exon.annotFeatures),
                              self.exon.chrom.subrange(chromStart, chromEnd),
                              self.exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)))
        self.exon.annotFeatures.append(gap)
        return gap.rna.end

    def _addUtr5Region(self, chromStart, chromEnd, rnaNext):
        utr5 = Utr5RegionFeature(self.exon, len(self.exon.annotFeatures),
                                 self.exon.chrom.subrange(chromStart, chromEnd),
                                 self.exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)))
        self.exon.annotFeatures.append(utr5)
        return utr5.rna.end

    def _addCdsRegion(self, chromStart, chromEnd, rnaNext, frame):
        cds = CdsRegionFeature(self.exon, len(self.exon.annotFeatures),
                               self.exon.chrom.subrange(chromStart, chromEnd),
                               self.exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)),
                               frame)
        self.exon.annotFeatures.append(cds)
        return cds.rna.end

    def _addUtr3Region(self, chromStart, chromEnd, rnaNext):
        utr3 = Utr3RegionFeature(self.exon, len(self.exon.annotFeatures),
                                 self.exon.chrom.subrange(chromStart, chromEnd),
                                 self.exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)))
        self.exon.annotFeatures.append(utr3)
        return utr3.rna.end

    def addCodingFeatures(self, chromStart, chromEnd, rnaNext, frame):
        """Add coding-related feature for a given contiguous block for the current
        exon.  A block may not cover whole exon if there are annotation
        gaps"""
        cds = self.exon.transcript.cdsChrom
        chromNext = chromStart
        if chromNext < cds.start:
            chromStop = min(chromEnd, cds.start)
            if self.exon.rna.strand == '+':
                rnaNext = self._addUtr5Region(chromNext, chromStop, rnaNext)
            else:
                rnaNext = self._addUtr3Region(chromNext, chromStop, rnaNext)
            chromNext = chromStop
        if (chromNext < chromEnd) and (chromNext < cds.end):
            chromStop = min(chromEnd, cds.end)
            rnaNext = self._addCdsRegion(chromNext, chromStop, rnaNext, frame)
            chromNext = chromStop
        if (chromNext < chromEnd) and (chromNext >= cds.end):
            if self.exon.rna.strand == '+':
                rnaNext = self._addUtr3Region(chromNext, chromEnd, rnaNext)
            else:
                rnaNext = self._addUtr5Region(chromNext, chromEnd, rnaNext)
        return rnaNext

    def addNonCodingFeature(self, chromStart, chromEnd, rnaNext):
        """Add a non-coding feature for a given contiguous block for the
        current exon.  A block may not cover whole exon if there are
        annotation gaps"""
        feat = NonCodingRegionFeature(self.exon, len(self.exon.annotFeatures),
                                      self.exon.chrom.subrange(chromStart, chromEnd),
                                      self.exon.rna.subrange(rnaNext, rnaNext + (chromEnd - chromStart)))
        self.exon.annotFeatures.append(feat)
        return feat.rna.end

    def _getSpliceSites(self, chromStart, chromEnd):
        if self.genomeReader is None:
            return (None, None)
        else:
            return spliceJuncsGetSeqs(self.genomeReader,
                                      self.transAnnot.chrom.name,
                                      chromStart, chromEnd,
                                      self.transAnnot.rna.strand)

    def addIntron(self, chromStart, chromEnd, rnaEnd):
        """Add an intron."""
        assert self.exon is None
        donorSeq, acceptorSeq = self._getSpliceSites(chromStart, chromEnd)
        intron = IntronFeature(self.transAnnot, len(self.transAnnot.features),
                               self.transAnnot.chrom.subrange(chromStart, chromEnd),
                               self.transAnnot.rna.subrange(rnaEnd, rnaEnd),
                               donorSeq, acceptorSeq)
        self.transAnnot.features.append(intron)
