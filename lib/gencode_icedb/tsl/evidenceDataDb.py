"""
Read evidence alignments from tabix files.
"""
import pysam
from pycbio.sys.symEnum import SymEnum, auto
from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.psl import Psl
from gencode_icedb.general.evidFeatures import EvidencePslFactory
from gencode_icedb.general.transFeatures import ExonFeature
import pipettor


class EvidenceSource(SymEnum):
    """Source of evidence used in support"""
    __slots__ = ()
    # FIXME: move to metadata or drop
    UCSC_RNA = auto()
    ENSEMBL_RNA = auto()
    MIXED_RNA = auto()
    UCSC_EST = auto()
    NANOPORE_DRNA = auto()
    NANOPORE_CDNA = auto()
    ISOSEQ_CDNA = auto()


def evidenceAlignsIndexPsl(pslFile):
    """Index the pslFile. It must be sorted."""
    pipettor.run(["tabix", "--force", "--sequence=14", "--begin=16", "--end=17", "--zero-based", pslFile])


class EvidenceAlignsReader(object):
    """Object for accessing overlapping alignment evidence data from a tabix file.
    """
    def __init__(self, evidSetUuid, evidPslTabix, genomeReader=None, genbankProblems=None):
        self.evidSetUuid = evidSetUuid
        self.tabix = pysam.TabixFile(evidPslTabix)
        self.contigs = frozenset(self.tabix.contigs)
        self.nameSubset = None  # used for testing and debugging.
        self.genbankProblems = genbankProblems
        self.evidFactory = EvidencePslFactory(genomeReader)

    def setNameSubset(self, nameSubset):
        """Set file on query names.  Can be a string, list, or set, or None to
        clear.  This is use for testing and debugging"""
        if isinstance(nameSubset, str):
            nameSubset = [nameSubset]
        if nameSubset is not None:
            frozenset(nameSubset)
        self.nameSubset = nameSubset

    def close(self):
        if self.tabix is not None:
            self.tabix.close()
            self.tabix = None

    def _makeTrans(self, psl):
        genbankProblem = self.genbankProblems.getProblem(psl.qName) if self.genbankProblems is not None else None
        attrs = ObjDict(genbankProblem=genbankProblem)
        return self.evidFactory.fromPsl(psl, attrs=attrs, orientChrom=True)

    def _usePsl(self, psl, strands):
        return (((self.nameSubset is None) or (psl.qName in self.nameSubset))
                and (psl.strand in strands))

    _posStrands = frozenset(('+', '++'))
    _negStrands = frozenset(('-', '+-', '-+', '--'))
    _allStrands = _posStrands.union(_negStrands)

    def _getSelectStrands(self, transcriptionStrand):
        # strand -- is used when PSLs of 3' ESTs have been reversed.
        if transcriptionStrand is None:
            return self._allStrands
        elif transcriptionStrand == '+':
            return self._posStrands
        else:
            return self._negStrands

    def _genOverlapping(self, coords, strands, minExons):
        for line in self.tabix.fetch(coords.name, coords.start, coords.end):
            psl = Psl.fromRow(line.split('\t'))
            if self._usePsl(psl, strands):
                trans = self._makeTrans(psl)
                if len(trans.getFeaturesOfType(ExonFeature)) >= minExons:
                    yield trans

    def genOverlapping(self, coords, transcriptionStrand=None, minExons=0):
        """Generator of overlapping alignments as TranscriptFeatures, possibly filtered
        by nameSubset.
        """
        if coords.name in self.contigs:
            yield from self._genOverlapping(coords, self._getSelectStrands(transcriptionStrand), minExons)
