"""
Annotation of a gene locus.
"""
from pycbio.hgdata.coords import Coords
from pycbio.sys.objDict import ObjDict


class GeneAnnotation(object):
    """Annotation of a gene locus.  This is a collection of TranscriptFeatures,
    along with bounds and attributes.  In the case of the PAR, one of these objects
    is created for each copy of the gene.

    Attributes, if supplied, they are of type ObjDict and thus addressable by
    field name. The attrs property will be None if not supply, so this must be
    checked first.

    Bounds are only for transcripts that have been added to object.
    """
    # FIXME: rename `chrom'
    __slots__ = ("geneId", "chrom", "transcriptionStrand", "attrs", "transcripts")

    def __init__(self, geneId, attrs=None):
        assert (attrs is None) or isinstance(attrs, ObjDict)
        self.geneId = geneId
        self.chrom = None
        self.transcriptionStrand = None
        self.attrs = attrs
        self.transcripts = []

    def _updateBounds(self, transAnnot):
        if transAnnot.chrom.name != self.chrom.name:
            raise Exception("Bug: mix of chromosomes provided: {} ({}) and {} ({})".format(transAnnot.rna.name, transAnnot.chrom.name, self.geneId, self.chrom.name))
        if transAnnot.chrom.strand != '+':
            raise Exception("Bug: assumes positive chromosome strand: {} and {}".format(transAnnot.rna.name, self.geneId))
        if transAnnot.transcriptionStrand != self.transcriptionStrand:
            raise Exception("Bug: mix of transcriptionStrand provided: {} and {}".format(transAnnot.rna.name, self.geneId))
        self.chrom = Coords(self.chrom.name,
                            min(transAnnot.chrom.start, self.chrom.start),
                            max(transAnnot.chrom.end, self.chrom.end),
                            self.chrom.strand, self.chrom.size)

    def addTranscript(self, transAnnot):
        assert transAnnot.geneAnnot is None
        if self.chrom is None:
            self.chrom = transAnnot.chrom
            self.transcriptionStrand = transAnnot.transcriptionStrand
        else:
            self._updateBounds(transAnnot)
        self.transcripts.append(transAnnot)
        transAnnot.geneAnnot = self


def _obtainGene(transAnnot, byGeneIdChrom):
    key = (transAnnot.chrom.name, transAnnot.attrs.geneId)
    geneAnnot = byGeneIdChrom.get(key)
    if geneAnnot is None:
        geneAnnot = byGeneIdChrom[key] = GeneAnnotation(transAnnot.attrs.geneId)
    return geneAnnot


def geneAnnotGroup(transAnnots):
    """Group TranscriptAnnotation objects into a list of GeneAnnotation objects."""
    # including chrom handles PAR
    byGeneIdChrom = {}
    for transAnnot in transAnnots:
        geneAnnot = _obtainGene(transAnnot, byGeneIdChrom)
        geneAnnot.addTranscript(transAnnot)
    return [byGeneIdChrom[key] for key in sorted(byGeneIdChrom.keys())]
