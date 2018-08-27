"""
Features of a transcript annotation or alignment.
"""
import sys
from copy import deepcopy
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.frame import Frame
from pycbio.hgdata.bed import Bed
from pycbio.sys.objDict import ObjDict
from gencode_icedb.general.spliceJuncs import SpliceJuncs, spliceJuncsClassify

# FIXME: this could easily be extended to handle non-transcript features, so change name.
# FIXME: should rna be None, rather than start/end=None,None


def _reverseComplementChildren(rcParent, features):
    "reverse complement a list of child Features, return None if features is None"
    if features is None:
        return None
    rcFeatures = []
    for feat in reversed(features):
        rcFeatures.append(feat.reverseComplement(rcParent, len(rcFeatures)))
    return tuple(rcFeatures)


def _getFeatureType0(featureTypes):
    """If featureTypes is a type, return it, otherwise return the first one"""
    return featureTypes if isinstance(featureTypes, type) else featureTypes[0]


def _getFeaturesOfType(feats, featureTypes):
    """Get of features that are of one of the specified types in the
    list of features..  The featureTypes argument can be a single type or a tuple of
    types"""
    return [f for f in feats if isinstance(f, featureTypes)]


class Feature(object):
    """A feature of annotation or alignment to the genome.  The RNA range is
    in genomic coordinates and it's length may not match the length of the
    feature in gapped feature types, due to gap closing.  For introns, are the
    interbase coordinates of the intron, which might not be equal if there are
    unaligned bases in the intron.  It is possible for either the chrom or RNA
    range to be None.

    Attributes maybe associated with any feature.  It supplied, they are of
    type ObjDict and thus addressable by field name. The attrs property will
    be None if not supply, so this must be checked first.
    """
    # FIXME: this should be named simply Feature
    # FIXME: not sure like the *Loc names
    __slots__ = ("parent", "iParent", "chrom", "rna", "attrs")
    name = None  # overridden by derived case

    def __init__(self, parent, iParent, chrom, rna, attrs):
        assert (chrom is None) or (isinstance(chrom, Coords) and (chrom.start is not None) and (chrom.end is not None))
        assert (rna is None) or (isinstance(rna, Coords) and (rna.start is not None) and (rna.end is not None))
        assert (attrs is None) or isinstance(attrs, ObjDict)
        self.parent = parent
        self.iParent = iParent
        self.chrom = chrom
        self.rna = rna
        self.attrs = attrs

    def __str__(self):
        "default str(), uses name"
        return "{} {}".format(self.name, self.coordsStr())

    def coordsStr(self):
        "get coordinates string"
        cBounds = (self.chrom.start, self.chrom.end) if self.chrom is not None else (None, None)
        rBounds = (self.rna.start, self.rna.end) if self.rna is not None else (None, None)
        return "{}-{} rna={}-{}".format(cBounds[0], cBounds[1], rBounds[0], rBounds[1])

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        return (str(self),)

    def _getChildrenStrTree(self, features):
        "build tuple for a list of children"
        r = []
        for feat in features:
            r.append(feat.toStrTree())
        return tuple(r)

    def dump(self, fh=sys.stderr, indent=0, msg=None):
        """print the tree for debugging purposes., optionally prefixing first line with
        msg and indenting beneath it"""
        if msg is not None:
            fh.write(msg + "\n")
        # use override-able implementation
        self._dumpImpl(fh, indent)

    def _dumpImpl(self, fh, indent):
        """print current node and children, assumes current indent has been done,
        This is overridden by derived classes as needed"""
        fh.write(str(self) + '\n')

    @staticmethod
    def _dumpChildren(fh, indent, features):
        "utility to dump a list of children"
        for feat in features:
            fh.write(indent * " ")
            feat._dumpImpl(fh, indent)

    @property
    def transcript(self):
        "get the transcript feature by walking the parents"
        if isinstance(self, TranscriptFeatures):
            return self
        elif self.parent is None:
            raise TypeError("no TranscriptFeatures ancestor found: {}".format(self))
        else:
            return self.parent.transcript

    def reverseComplement(self, rcParent, iRcParent):
        "default reverseComplement when __init__ takes just coordinates"
        return self.__class__(rcParent, iRcParent, self.chrom.reverse(), self.rna.reverse(), deepcopy(self.attrs))

    def nextFeature(self, featureTypes=None):
        """Return the next feature in the sequence of features of one of the
        given types.  This maybe a child of a different parent. The
        featureTypes argument can be one type or a tuple of types"""
        if featureTypes is None:
            featureTypes = Feature
        feat = self
        while True:
            feat = feat._nextFeatureImpl()
            if (feat is None) or isinstance(feat, featureTypes):
                return feat

    def prevFeature(self, featureTypes=None):
        """Return the previous feature in the sequence of features of one of the
        given types.  This maybe a child of a different parent. The
        featureTypes argument can be one type or a tuple of types"""
        if featureTypes is None:
            featureTypes = Feature
        feat = self
        while True:
            feat = feat._prevFeatureImpl()
            if (feat is None) or (featureTypes is None) or isinstance(feat, featureTypes):
                return feat

    def _nextFeatureImpl(self):
        raise TypeError("nextFeature not valid on features of type {}".format(type(self).__name__))

    def _prevFeatureImpl(self):
        raise TypeError("prevFeature not valid on features of type {}".format(type(self).__name__))
        self._prevNextNotValue()


class AlignmentFeature(Feature):
    "ABC for alignment features"
    def __init__(self, parent, iParent, chrom, rna, attrs=None):
        assert isinstance(parent, StructureFeature)
        super(AlignmentFeature, self).__init__(parent, iParent, chrom, rna, attrs)

    def _nextFeatureImpl(self):
        """Return the next feature in this sequence of AlignmentFeatures."""
        if self.iParent < len(self.parent.alignFeatures) - 1:
            return self.parent.alignFeatures[self.iParent + 1]
        else:
            pnext = self.parent
            while True:
                pnext = pnext._nextFeatureImpl()
                if pnext is None:
                    return None
                elif len(pnext.alignFeatures) > 0:
                    return pnext.alignFeatures[0]

    def _prevFeatureImpl(self):
        """Return the next feature in this sequence of AlignmentFeatures.  This
        maybe a child of a different StructureFeature."""
        if self.iParent > 0:
            return self.parent.alignFeatures[self.iParent - 1]
        else:
            pprev = self.parent
            while True:
                pprev = pprev._prevFeatureImpl()
                if pprev is None:
                    return None
                elif len(pprev.alignFeatures) > 1:
                    return pprev.alignFeatures[-1]


class AlignedFeature(AlignmentFeature):
    """Ungapped alignment block."""
    name = "aln"
    __slots__ = ()


class ChromInsertFeature(AlignmentFeature):
    """Chromosome insert region (or RNA deletion)."""
    name = "cins"
    __slots__ = ()

    def __init__(self, parent, iParent, chrom, attrs=None):
        super(ChromInsertFeature, self).__init__(parent, iParent, chrom, None, attrs)

    def reverseComplement(self, rcParent, iRcParent):
        return ChromInsertFeature(rcParent, iRcParent, self.chrom.reverse(), deepcopy(self.attrs))


class RnaInsertFeature(AlignmentFeature):
    """RNA insert region (or chrom deletion)."""
    name = "rins"
    __slots__ = ()

    def __init__(self, parent, iParent, rna, attrs=None):
        super(RnaInsertFeature, self).__init__(parent, iParent, None, rna, attrs)

    def reverseComplement(self, rcParent, iRcParent):
        return RnaInsertFeature(rcParent, iRcParent, self.rna.reverse(), deepcopy(self.attrs))


class AnnotationFeature(Feature):
    "ABC for annotation features"
    def __init__(self, parent, iParent, chrom, rna, attrs=None):
        super(AnnotationFeature, self).__init__(parent, iParent, chrom, rna, attrs)

    def _nextFeatureImpl(self):
        """Return the next feature in this sequence of AnnotationFeatures."""
        if self.iParent < len(self.parent.annotFeatures) - 1:
            return self.parent.annotFeatures[self.iParent + 1]
        else:
            pnext = self.parent
            while True:
                pnext = pnext._nextFeatureImpl()
                if pnext is None:
                    return None
                elif len(pnext.annotFeatures) > 0:
                    return pnext.annotFeatures[0]

    def _prevFeatureImpl(self):
        """Return the next feature in this sequence of AnnotationFeatures."""
        if self.iParent > 0:
            return self.parent.annotFeatures[self.iParent - 1]
        else:
            pprev = self.parent
            while True:
                pprev = pprev._prevFeatureImpl()
                if pprev is None:
                    return None
                elif len(pprev.annotFeatures) > 0:
                    return pprev.annotFeatures[-1]


class GapAnnotFeature(AnnotationFeature):
    """A gap in CDS annotation for unspecified reasons, often causes by join small gaps in exons"""
    name = "gap"
    __slots__ = ()

    def __init__(self, parent, iParent, chrom, rna, attrs=None):
        super(GapAnnotFeature, self).__init__(parent, iParent, chrom, rna, attrs)


class Utr5RegionFeature(AnnotationFeature):
    "A un-gapped 5'UTR region in an exon"
    name = "5'UTR"
    __slots__ = ()


class CdsRegionFeature(AnnotationFeature):
    "A un-gapped CDS region in an exon"
    name = "CDS"
    __slots__ = ("frame")

    def __init__(self, parent, iParent, chrom, rna, frame, attrs=None):
        assert(isinstance(frame, Frame))
        super(CdsRegionFeature, self).__init__(parent, iParent, chrom, rna, attrs)
        self.frame = frame

    def __str__(self):
        return "{} {} {}".format(self.name, self.coordsStr(), self.frame)

    def reverseComplement(self, rcParent, iRcParent):
        rcRnaLoc = self.rna.reverse()
        return CdsRegionFeature(rcParent, iRcParent, self.chrom.reverse(), rcRnaLoc, self.frame, deepcopy(self.attrs))


class Utr3RegionFeature(AnnotationFeature):
    "A un-gapped 3'UTR region in an exon"
    name = "3'UTR"
    __slots__ = ()


class NonCodingRegionFeature(AnnotationFeature):
    "An un-gapped non-coding region in an exon"
    name = "NC"
    __slots__ = ()


class StructureFeature(Feature):
    "ABC for structural features"
    __slots__ = ("annotFeatures", "alignFeatures")

    def __init__(self, parent, iParent, chrom, rna, attrs=None):
        assert isinstance(parent, TranscriptFeatures)
        super(StructureFeature, self).__init__(parent, iParent, chrom, rna, attrs)
        self.annotFeatures = ()
        self.alignFeatures = ()

    def _nextFeatureImpl(self):
        """Return the next feature in this sequence of StructureFeatures."""
        if self.iParent < len(self.parent.features) - 1:
            return self.parent.features[self.iParent + 1]
        else:
            return None

    def _prevFeatureImpl(self):
        """Return the next feature in this sequence of StructureFeatures."""
        if self.iParent > 0:
            return self.parent.features[self.iParent - 1]
        else:
            return None

    def _getChildList(self, featureTypes):
        featureType0 = _getFeatureType0(featureTypes)
        if issubclass(featureType0, AlignmentFeature):
            return self.alignFeatures
        elif issubclass(featureType0, AnnotationFeature):
            return self.annotFeatures
        else:
            raise TypeError("Type '{}' is a valid child feature type of '{}', must specify an '{}' or an '{}' type or subtype".
                            format(featureType0.__name__, type(self).__name__,
                                   AlignmentFeature.__name__, AnnotationFeature.__name__))

    def getFeaturesOfType(self, featureTypes):
        """Get child of features of the specified type.  These must all be
        either AlignmentFeature or AnnotationFeature types."""
        return _getFeaturesOfType(self._getChildList(featureTypes), featureTypes)

    def firstFeature(self, featureTypes=Feature):
        """Get the first child of feature of the specified type.  These must all be
        either AlignmentFeature or AnnotationFeature types."""
        for feat in self._getChildList(featureTypes):
            if isinstance(feat, featureTypes):
                return feat
        return None

    def lastFeature(self, featureTypes=Feature):
        """Get the last chilkd feature of the specified type.  These must all be
        either AlignmentFeature or AnnotationFeature types."""
        for feat in reversed(self._getChildList(featureTypes)):
            if isinstance(feat, featureTypes):
                return feat
        return None

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.annotFeatures:
            r.append(self._getChildrenStrTree(self.annotFeatures))
        if self.alignFeatures:
            r.append(self._getChildrenStrTree(self.alignFeatures))
        return tuple(r)

    def _dumpImpl(self, fh, indent):
        fh.write(str(self) + '\n')
        self._dumpChildren(fh, indent + 2, self.annotFeatures)
        self._dumpChildren(fh, indent + 2, self.alignFeatures)


class ExonFeature(StructureFeature):
    """exon with target gaps closed."""
    name = "exon"

    def __init__(self, parent, iParent, chrom, rna, attrs=None):
        super(ExonFeature, self).__init__(parent, iParent, chrom, rna, attrs)

    def reverseComplement(self, rcParent, iRcParent):
        rcExon = ExonFeature(rcParent, iRcParent, self.chrom.reverse(), self.rna.reverse(), deepcopy(self.attrs))
        rcExon.annotFeatures = _reverseComplementChildren(rcExon, self.annotFeatures)
        rcExon.alignFeatures = _reverseComplementChildren(rcExon, self.alignFeatures)
        return rcExon

    def rnaOverlaps(self, feat2):
        "does RNA range overlap another exon?"
        return (self.rna.start < feat2.rna.end) and (self.rna.end > feat2.rna.start)

    def _countAlignedBases(self):
        alignedCnt = 0
        for feat in self.alignFeatures:
            if isinstance(feat, AlignedFeature):
                alignedCnt += len(feat.rna)
        return alignedCnt

    @property
    def alignedBases(self):
        """if there are alignment subfeatures, return the actual number of
        aligned bases, otherwise, the RNA for an annotation"""
        if len(self.alignFeatures) > 0:
            return self._countAlignedBases()
        else:
            return len(self.rna)


class IntronFeature(StructureFeature):
    """Intron from annotation or alignment, splice junction information maybe
    None. The alignFeatures field are for query insertions in intron, and maybe None.
    The donorSeq and acceptorSeq values are in the direction of transcription and
    are lower-case if splice junction motif is unknown and upper case if known.
    The spliceJuncs field is a SpliceJuncs object or None.
    """
    name = "intron"
    __slots__ = ("donorSeq", "acceptorSeq", "spliceJuncs")

    def __init__(self, parent, iParent, chrom, rna, donorSeq, acceptorSeq, attrs=None):
        super(IntronFeature, self).__init__(parent, iParent, chrom, rna, attrs)
        self.donorSeq = self.acceptorSeq = self.spliceJuncs = None
        if donorSeq is not None:
            self.spliceJuncs = spliceJuncsClassify(donorSeq, acceptorSeq)
            if self.spliceJuncs == SpliceJuncs.unknown:
                self.donorSeq, self.acceptorSeq = donorSeq.lower(), acceptorSeq.lower()
            else:
                self.donorSeq, self.acceptorSeq = donorSeq.upper(), acceptorSeq.upper()

    def __str__(self):
        if self.donorSeq is None:
            sjDesc = ""
        else:
            sjDesc = " sjBases={}...{} ({})".format(self.donorSeq, self.acceptorSeq, str(self.spliceJuncs))
        return "intron {}{}".format(self.coordsStr(), sjDesc)

    @property
    def sjBases(self):
        """Get splice junction patterns, or None if not available"""
        if self.donorSeq is None:
            return None
        else:
            return "{}/{}".format(self.donorSeq, self.acceptorSeq)

    def reverseComplement(self, rcParent, iRcParent):
        rcIntron = IntronFeature(rcParent, iRcParent, self.chrom.reverse(), self.rna.reverse(),
                                 self.donorSeq, self.acceptorSeq, deepcopy(self.attrs))
        rcIntron.alignFeatures = _reverseComplementChildren(rcIntron, self.alignFeatures)
        return rcIntron

    def rnaIntersect(self, feat2):
        "does RNA interbase range intersect another feature? (overlap allowing zero length)"
        return (self.rna.start <= feat2.rna.end) and (self.rna.end >= feat2.rna.start)


class TranscriptFeatures(Feature):
    """
    Set of features for a transcript derived from an alignment or annotation.
    The transcriptionStrand column is the direction of transcription.  For alignments,
    this may differ from the rna strand for 3' ESTs.
    """
    name = "trans"
    __slots__ = ("chrom", "rna", "transcriptionStrand", "cdsChrom", "features", "geneAnnot")

    def __init__(self, chrom, rna, transcriptionStrand, cdsChrom=None, attrs=None):
        super(TranscriptFeatures, self).__init__(None, None, chrom, rna, attrs)
        if (chrom.size is not None) and (cdsChrom is not None) and (cdsChrom.size is None):
            raise TypeError("if chrom has size, cdsChrom must have size if specified")
        self.chrom = chrom
        self.rna = rna
        self.transcriptionStrand = transcriptionStrand
        self.cdsChrom = cdsChrom
        self.features = None
        self.geneAnnot = None  # will be set if contained in a GeneAnnotation object

    def __str__(self):
        return "t={}/{}, rna={}/{} {} <{}>".format(str(self.chrom), self.chrom.strand,
                                                   str(self.rna), self.rna.strand, self.rna.size,
                                                   self.transcriptionStrand)

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.features:
            r.append(self._getChildrenStrTree(self.features))
        return tuple(r)

    def _dumpImpl(self, fh, indent):
        fh.write(str(self) + '\n')
        self._dumpChildren(fh, indent + 2, self.features)

    def toBed(self, itemRgb=""):
        """convert transcript and CDS to a Bed object"""
        # FIXME: shouldn't me method, but I/O object
        blocks = []
        for feat in self.features:
            if isinstance(feat, ExonFeature):
                blocks.append(Bed.Block(feat.chrom.start, feat.chrom.end))
        if self.cdsChrom is not None:
            cdsStart, cdsEnd = self.cdsChrom.start, self.cdsChrom.end
        else:
            cdsStart = cdsEnd = self.chrom.end
        return Bed(self.chrom.name, self.chrom.start, self.chrom.end,
                   self.rna.name, 0, self.rna.strand, cdsStart, cdsEnd,
                   itemRgb, blocks)

    @property
    def alignedBases(self):
        alignedCnt = 0
        for feature in self.features:
            if isinstance(feature, ExonFeature):
                alignedCnt += feature.alignedBases
        return alignedCnt

    def reverseComplement(self):
        "return a new TranscriptFeatures object that is reverse complemented"
        if self.chrom.size is None:
            raise TypeError("can't reverse-complement transcript without chrom.size: {}".format(self))
        if self.rna.size is None:
            raise TypeError("can't reverse-complement transcript without rna.size: {}".format(self))
        rcCdsChromLoc = self.cdsChrom.reverse() if self.cdsChrom is not None else None

        rcTrans = TranscriptFeatures(self.chrom.reverse(), self.rna.reverse(),
                                     self.transcriptionStrand, rcCdsChromLoc, deepcopy(self.attrs))
        rcTrans.features = _reverseComplementChildren(rcTrans, self.features)
        return rcTrans

    def getFeaturesOfType(self, featureTypes):
        """Get child of features of the specified type. Types should all be
        valid sibling types of each other."""
        featureType0 = _getFeatureType0(featureTypes)
        if issubclass(featureType0, StructureFeature):
            return _getFeaturesOfType(self.features, featureTypes)
        else:
            return [af for sf in self.features
                    for af in sf.getFeaturesOfType(featureTypes)]

    def firstFeature(self, featureTypes=StructureFeature):
        """Get the first child of feature of the specified type. Types should all be
        valid sibling types of each other.  If a StructureFeature is not requests,
        checks the children of the structure features.. """
        featureType0 = _getFeatureType0(featureTypes)
        if issubclass(featureType0, StructureFeature):
            for sf in self.features:
                if isinstance(sf, featureTypes):
                    return sf
        else:
            for sf in self.features:
                feat = sf.firstFeature(featureTypes)
                if feat is not None:
                    return feat
        return None

    def lastFeature(self, featureTypes=Feature):
        """Get the last child of feature of the specified type. Types should all be
        valid sibling types of each other."""
        featureType0 = _getFeatureType0(featureTypes)
        if issubclass(featureType0, StructureFeature):
            for sf in reversed(self.features):
                if isinstance(sf, featureTypes):
                    return sf
        else:
            for sf in reversed(self.features):
                feat = sf.lastFeature(featureTypes)
                if feat is not None:
                    return feat
            return None
