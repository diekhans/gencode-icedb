"""
Common code used to build a TranscriptFeature object for an annotation.
"""
from pycbio.hgdata.coords import Coords


class AnnotFeatureBuilder(object):
    def __init__(self, chrom, rna, cdsChrom, transcriptionStrand, attrs):
        self.transAnnot = TranscriptFeatures(chrom, rna, transcriptionStrand=transcriptionStrand, cdsChrom=cdsChrom, attrs=attrs)
        self.transFeatures = []
        self.exon = None
        self.exonFeatures = None

    def finish(self):
        "finish off transcript"
