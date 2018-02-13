"""
Features of a transcript annotation or alignment.
"""
from __future__ import print_function
from copy import deepcopy
from pycbio.hgdata.coords import Coords
from pycbio.hgdata import dnaOps
from pycbio.hgdata.frame import Frame
from pycbio.hgdata.bed import Bed
from pycbio.sys.objDict import ObjDict
from gencode_icedb.general.spliceJuncs import SpliceJuncs, spliceJuncsClassify


def _reverseComplementChildren(rcParent, features):
    "reverse complement a list of child TransFeatures, return None if features is None"
    if features is None:
        return None
    rcFeatures = []
    for i in range(len(features) - 1, -1, -1):
        rcFeatures.append(features[i].reverseComplement(rcParent, len(rcFeatures)))
    return tuple(rcFeatures)


def _getFeaturesOfType(feats, types):
    """Get of features that are of one of the specified types in the
    list of features..  The types argument can be one type of a list of
    types"""
    if isinstance(types, type):
        types = (types,)
    return [f for f in feats if isinstance(f, types)]


class TransFeature(object):
    """A feature of annotation or alignment to the genome.  The RNA range is
    in genomic coordinates and it's length may not match the length of the
    feature in gapped feature types, due to gap closing.  For introns, are the
    interbase coordinates of the intron, which might not be equal if there are
    unaligned bases in the intron.  It is possible for either the chrom or RNA
    range to be None.

    Attributes maybe associated with any feature.  It supplied, they are of
    type ObjDict and thus addressable by field name. The attrs property will
    be None if not supply, so this much be checked first.
    """
    __slots__ = ("parent", "iParent", "chromLoc", "rnaLoc", "attrs")
    name = None  # overridden by derived case

    def __init__(self, parent, iParent, chromLoc, rnaLoc, attrs):
        assert (chromLoc is None) or (isinstance(chromLoc, Coords) and (chromLoc.start is not None) and (chromLoc.end is not None))
        assert (rnaLoc is None) or (isinstance(rnaLoc, Coords) and (rnaLoc.start is not None) and (rnaLoc.end is not None))
        assert (attrs is None) or isinstance(attrs, ObjDict)
        self.parent = parent
        self.iParent = iParent
        self.chromLoc = chromLoc
        self.rnaLoc = rnaLoc
        self.attrs = attrs

    def __str__(self):
        "default str(), uses name"
        return "{} {}".format(self.name, self.coordsStr())

    def coordsStr(self):
        "get coordinates string"
        cBounds = (self.chromLoc.start, self.chromLoc.end) if self.chromLoc is not None else (None, None)
        rBounds = (self.rnaLoc.start, self.rnaLoc.end) if self.rnaLoc is not None else (None, None)
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

    def dump(self, fh, indent=0, msg=None):
        """print the tree for debugging purposes., optionally prefixing first line with
        msg and indenting beneath it"""
        if msg is not None:
            pmsg = msg + ": "
            fh.write(pmsg)
            indent += len(pmsg)
        self._dumpDerive(fh, indent)

    def _dumpDerive(self, fh, indent):
        """print current node and children, assumes current indent has been done,
        This is overridden by derived classes as needed"""
        fh.write(str(self) + '\n')

    @staticmethod
    def _dumpChildren(fh, indent, features):
        "utility to dump a list of children"
        for feat in features:
            fh.write(indent * " ")
            feat._dumpDerive(fh, indent)

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
        return self.__class__(rcParent, iRcParent, self.chromLoc.reverse(), self.rnaLoc.reverse(), deepcopy(self.attrs))


class AlignmentFeature(TransFeature):
    "ABC for alignment features"
    def __init__(self, parent, iParent, chromLoc, rnaLoc, attrs=None):
        assert isinstance(parent, StructureFeature)
        super(AlignmentFeature, self).__init__(parent, iParent, chromLoc, rnaLoc, attrs)

    def prevFeat(self):
        """Return the next feature in this sequence of AlignmentFeatures.  This
        maybe a child of a different StructureFeature."""
        if self.iParent > 0:
            return self.parent.alignFeatures[self.iParent - 1]
        else:
            pprev = self.parent.prevFeat()
            if pprev is not None:
                return pprev.alignFeatures[:-1]
            else:
                return None

    def nextFeat(self):
        """Return the next feature in this sequence of AlignmentFeatures.
        This maybe a child of a different StructureFeature."""
        if self.iParent < len(self.parent.alignFeatures) - 1:
            return self.parent.alignFeatures[self.iParent + 1]
        else:
            pnext = self.parent.nextFeat()
            if pnext is not None:
                return pnext.alignFeatures[0]
            else:
                return None


class AlignedFeature(AlignmentFeature):
    """Ungapped alignment block."""
    name = "aln"
    __slots__ = ()


class ChromInsertFeature(AlignmentFeature):
    """Chromosome insert region (or RNA deletion)."""
    name = "cins"
    __slots__ = ()

    def __init__(self, parent, iParent, chromLoc, attrs=None):
        super(ChromInsertFeature, self).__init__(parent, iParent, chromLoc, None, attrs)

    def reverseComplement(self, rcParent, iRcParent):
        return ChromInsertFeature(rcParent, iRcParent, self.chromLoc.reverse(), deepcopy(self.attrs))


class RnaInsertFeature(AlignmentFeature):
    """RNA insert region (or chrom deletion)."""
    name = "rins"
    __slots__ = ()

    def __init__(self, parent, iParent, rnaLoc, attrs=None):
        super(RnaInsertFeature, self).__init__(parent, iParent, None, rnaLoc, attrs)

    def reverseComplement(self, rcParent, iRcParent):
        return RnaInsertFeature(rcParent, iRcParent, self.rnaLoc.reverse(), deepcopy(self.attrs))


class AnnotationFeature(TransFeature):
    "ABC for annotation features"
    def __init__(self, parent, iParent, chromLoc, rnaLoc, attrs=None):
        super(AnnotationFeature, self).__init__(parent, iParent, chromLoc, rnaLoc, attrs)

    def prevFeat(self):
        """Return the next feature in this sequence of AnnotationFeatures.  This maybe a
        child of a different StructureFeature."""
        if self.iParent > 0:
            return self.parent.annotFeatures[self.iParent - 1]
        else:
            pprev = self.parent.prevFeat()
            if pprev is not None:
                return pprev.annotFeatures[:-1]
            else:
                return None

    def nextFeat(self):
        """Return the next feature in this sequence of AnnotationFeatures.
        This maybe a child of a different StructureFeature."""
        if self.iParent < len(self.parent.annotFeatures) - 1:
            return self.parent.annotFeatures[self.iParent + 1]
        else:
            pnext = self.parent.nextFeat()
            if pnext is not None:
                return pnext.annotFeatures[0]
            else:
                return None


class Utr5RegionFeature(AnnotationFeature):
    "A un-gapped 5'UTR region in an exon"
    name = "5'UTR"
    __slots__ = ()


class CdsRegionFeature(AnnotationFeature):
    "A un-gapped CDS region in an exon"
    name = "CDS"
    __slots__ = ("frame")

    def __init__(self, parent, iParent, chromLoc, rnaLoc, frame, attrs=None):
        assert(isinstance(frame, Frame))
        super(CdsRegionFeature, self).__init__(parent, iParent, chromLoc, rnaLoc, attrs)
        self.frame = frame

    def __str__(self):
        return "{} {} {}".format(self.name, self.coordsStr(), self.frame)

    def reverseComplement(self, rcParent, iRcParent):
        rcRnaLoc = self.rnaLoc.reverse()
        rcFrame = self.frame + len(rcRnaLoc)
        return CdsRegionFeature(rcParent, iRcParent, self.chromLoc.reverse(), rcRnaLoc, rcFrame, deepcopy(self.attrs))


class Utr3RegionFeature(AnnotationFeature):
    "A un-gapped 3'UTR region in an exon"
    name = "3'UTR"
    __slots__ = ()


class NonCodingRegionFeature(AnnotationFeature):
    "An un-gapped non-coding region in an exon"
    name = "NC"
    __slots__ = ()


class StructureFeature(TransFeature):
    "ABC for structural features"
    __slots__ = ("annotFeatures", "alignFeatures")

    def __init__(self, parent, iParent, chromLoc, rnaLoc, attrs=None):
        assert isinstance(parent, TranscriptFeatures)
        super(StructureFeature, self).__init__(parent, iParent, chromLoc, rnaLoc, attrs)
        self.annotFeatures = ()
        self.alignFeatures = ()

    def prevFeat(self):
        """Return the next feature in this sequence of StructureFeatures."""
        if self.iParent > 0:
            return self.parent.features[self.iParent - 1]
        else:
            return None

    def nextFeat(self):
        """Return the next feature in this sequence of StructureFeatures."""
        if self.iParent < len(self.parent.features) - 1:
            return self.parent.features[self.iParent + 1]
        else:
            return None

    def getAnnotationFeaturesOfType(self, types):
        """Get of AnnotationFeatures that are of one of the specified types in the
        list of features..  The types argument can be one type of a list of
        types"""
        return _getFeaturesOfType(self.annotFeatures, types)

    def getAlignmentFeaturesOfType(self, types):
        """Get of AlignmentFeatures that are of one of the specified types in the
        list of features..  The types argument can be one type of a list of
        types"""
        return _getFeaturesOfType(self.alignFeatures, types)

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.annotFeatures:
            r.append(self._getChildrenStrTree(self.annotFeatures))
        if self.alignFeatures:
            r.append(self._getChildrenStrTree(self.alignFeatures))
        return tuple(r)

    def _dumpDerived(self, fh, indent):
        fh.write(str(self) + '\n')
        self._dumpChildren(fh, indent + 2, self.annotFeatures)
        self._dumpChildren(fh, indent + 2, self.alignFeatures)


class ExonFeature(StructureFeature):
    """exon with target gaps closed."""
    name = "exon"

    def __init__(self, parent, iParent, chromLoc, rnaLoc, attrs=None):
        super(ExonFeature, self).__init__(parent, iParent, chromLoc, rnaLoc, attrs)

    def reverseComplement(self, rcParent, iRcParent):
        rcExon = ExonFeature(rcParent, iRcParent, self.chromLoc.reverse(), self.rnaLoc.reverse(), deepcopy(self.attrs))
        rcExon.annotFeatures = _reverseComplementChildren(rcExon, self.annotFeatures)
        rcExon.alignFeatures = _reverseComplementChildren(rcExon, self.alignFeatures)
        return rcExon

    def rnaOverlaps(self, feat2):
        "does RNA range overlap another exon?"
        return (self.rnaLoc.start < feat2.rnaLoc.end) and (self.rnaLoc.end > feat2.rnaLoc.start)

    def _countAlignedBases(self):
        alignedCnt = 0
        for feat in self.alignFeatures:
            if isinstance(feat, AlignedFeature):
                alignedCnt += len(feat.rnaLoc)
        return alignedCnt

    @property
    def alignedBases(self):
        """if there are alignment subfeatures, return the actual number of
        aligned bases, otherwise, the RNA for an annotation"""
        if len(self.alignFeatures) > 0:
            return self._countAlignedBases()
        else:
            return len(self.rnaLoc)


class IntronFeature(StructureFeature):
    """Intron from annotation or alignment, splice junction information maybe
    None. The alignFeatures field are for query insertions in intron, and maybe None.
    The donorSeq and acceptorSeq values are in the direction of transcription and
    are lower-case if splice junction motif is unknown and upper case if known.
    The spliceJuncs field is a SpliceJuncs object or None.
    """
    name = "intron"
    __slots__ = ("donorSeq", "acceptorSeq", "spliceJuncs")

    def __init__(self, parent, iParent, chromLoc, rnaLoc, donorSeq, acceptorSeq, attrs=None):
        super(IntronFeature, self).__init__(parent, iParent, chromLoc, rnaLoc, attrs)
        self.donorSeq = self.acceptorSeq = self.spliceJuncs = None
        if donorSeq is not None:
            self.spliceJuncs = spliceJuncsClassify(donorSeq, acceptorSeq)
            if self.spliceJuncs == SpliceJuncs.unknown:
                self.donorSeq, self.acceptorSeq = donorSeq.lower(), acceptorSeq.lower()
            else:
                self.donorSeq, self.acceptorSeq = donorSeq.upper(), acceptorSeq.upper()

    def __str__(self):
        if self.donorSeq is None:
            sjDesc = None
        else:
            sjDesc = "{}...{} ({})".format(self.donorSeq, self.acceptorSeq, str(self.spliceJuncs))
        return "intron {} sjBases={}".format(self.coordsStr(), sjDesc)

    @property
    def sjBases(self):
        """Get splice junction patterns, or None if not available"""
        if self.donorSeq is None:
            return None
        else:
            return "{}/{}".format(self.donorSeq, self.acceptorSeq)

    def reverseComplement(self, rcParent, iRcParent):
        rcIntron = IntronFeature(rcParent, iRcParent, self.chromLoc.reverse(), self.rnaLoc.reverse(),
                                 self.donorSeq, self.acceptorSeq, deepcopy(self.attrs))
        rcIntron.alignFeatures = _reverseComplementChildren(rcIntron, self.alignFeatures)
        return rcIntron

    def rnaIntersect(self, feat2):
        "does RNA interbase range intersect another feature? (overlap allowing zero length)"
        return (self.rnaLoc.start <= feat2.rnaLoc.end) and (self.rnaLoc.end >= feat2.rnaLoc.start)


class TranscriptFeatures(TransFeature):
    """
    Set of features for a transcript derived from an alignment or annotation,
    features are kept in chromosome order (positive strand).
    """
    name = "trans"
    __slots__ = ("chromLoc", "rnaLoc", "cdsChromLoc", "features")

    def __init__(self, chromLoc, rnaLoc, cdsChromLoc=None, attrs=None):
        super(TranscriptFeatures, self).__init__(None, None, chromLoc, rnaLoc, attrs)
        if (chromLoc.size is not None) and (cdsChromLoc is not None) and (cdsChromLoc.size is None):
            raise Exception("if chromLoc has size, cdsChromLoc must have size if specified")
        self.chromLoc = chromLoc
        self.rnaLoc = rnaLoc
        self.cdsChromLoc = cdsChromLoc
        self.features = None

    def __str__(self):
        return "t={}/{}, rna={}/{} {}".format(str(self.chromLoc), self.chromLoc.strand,
                                              str(self.rnaLoc), self.rnaLoc.strand, self.rnaLoc.size)

    def toStrTree(self):
        """recursively convert to a recursive tuple of strings representing
        the feature tree"""
        r = [str(self)]
        if self.features:
            r.append(self._getChildrenStrTree(self.features))
        return tuple(r)

    def _dumpDerived(self, fh, indent):
        fh.write(str(self) + '\n')
        self._dumpChildren(fh, indent + 2, self.features)

    def toBed(self, itemRgb=""):
        """convert transcript and CDS to a Bed object"""
        # FIXME: shouldn't me method, but I/O object
        blocks = []
        for feat in self.features:
            if isinstance(feat, ExonFeature):
                blocks.append(Bed.Block(feat.chromLoc.start, feat.chromLoc.end))
        if self.cdsChromLoc is not None:
            cdsStart, cdsEnd = self.cdsChromLoc.start, self.cdsChromLoc.end
        else:
            cdsStart = cdsEnd = self.chromLoc.end
        return Bed(self.chromLoc.name, self.chromLoc.start, self.chromLoc.end,
                   self.rnaLoc.name, 0, self.rnaLoc.strand, cdsStart, cdsEnd,
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
        if self.chromLoc.size is None:
            raise Exception("can't reverse-complement transcript without chromLoc.size: {}".format(self))
        rcCdsChromLoc = self.cdsChromLoc.reverse() if self.cdsChromLoc is not None else None

        rcTrans = TranscriptFeatures(self.chromLoc.reverse(), self.rnaLoc.reverse(),
                                     rcCdsChromLoc, deepcopy(self.attrs))
        rcTrans.features = _reverseComplementChildren(rcTrans, self.features)
        return rcTrans

    def getStructureFeaturesOfType(self, types):
        """Get of StructureFeatures that are of one of the specified types in the
        list of features..  The types argument can be one type of a list of
        types"""
        return _getFeaturesOfType(self.features, types)

    def getAnnotationFeaturesOfType(self, types):
        """Get of AnnotationFeatures that are of one of the specified types in the
        list of features..  The types argument can be one type of a list of
        types"""
        return [af for sf in self.features
                for af in sf.getAnnotationFeaturesOfType(types)]

    def getAlignmentFeaturesOfType(self, types):
        """Get of AlignmentFeatures that are of one of the specified types in the
        list of features..  The types argument can be one type of a list of
        types"""
        return [af for sf in self.features
                for af in sf.getAlignmentFeaturesOfType(types)]
