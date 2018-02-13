"""
Analyze transcript annotations against evidence.

The current algorithm is:
    - tslNA:  transcripts that are:
              - Olfactory receptors
              - HLA transcripts
              - Immunoglobins
              - T-Cell receptors
              - Single exon genes
  - multi-exon transcripts:
     - tsl1 - all splice junctions of the transcript are supported by at least one non-suspect mRNA.
     - tsl2 - the best supporting mRNA is flagged as suspect or the support is from multiple ESTs
     - tsl3 - the only support is from a single EST
     - tsl4 - the best supporting EST is flagged as suspect
     - tsl5 - no single transcript supports the model structure
"""
import re
from collections import namedtuple, defaultdict
from pycbio.sys import fileOps
from gencode_icedb.general.transFeatures import ExonFeature, ChromInsertFeature, RnaInsertFeature
from gencode_icedb.tsl.supportDefs import EvidenceSupport, TrascriptionSupportLevel
from gencode_icedb.tsl.evidenceDb import EvidenceSource
from gencode_icedb.general.gencodeDb import findAnnotationBounds


# limits on size of as single indel in an exon.
exonPolymorhicSizeLimit = 12
# fraction of allowed total indel size relative exon length
exonPolymorhicFactionLimit = 0.05


# FIXME these are from the ccds2/modules/gencode/src/lib/gencode/data/gencodeGenes.py, migrate to new module
def _transIsSingleExon(annotTrans):
    return len(annotTrans.getStructureFeaturesOfType(ExonFeature)) <= 1


def _geneIsTslIgnored(annotTrans):
    bioType = annotTrans.attrs.transcriptType
    geneName = annotTrans.attrs.geneName
    # FIXME: this was part of ccds gencodeGenes module, we need to make that independent and use it here
    # trans.isPseudo()
    if (bioType != "polymorphic_pseudogene") and (bioType.find("pseudogene") >= 0) and (bioType.find("transcribed") < 0):
        return True
    # trans.isOlfactoryReceptor()
    if re.match("^OR[0-9].+", geneName):
        return True
    # trans.isHLA()
    if geneName.find("HLA-") == 0:
        return True
    # trans.isImmunoglobin()
    if bioType.startswith("IG_"):
        return True
    # trans.isTCellReceptor())
    if bioType.startswith("TR_"):
        return True
    return False


def _checkMegExonIndels(evidExon):
    "check for allowed indel polymorphism"

    def getIndelSize(aln):
        if isinstance(aln, ChromInsertFeature):
            return len(aln.chromLoc)
        elif isinstance(aln, RnaInsertFeature):
            return len(aln.rnaLoc)
        else:
            return 0  # not an indel

    totalIndelSize = 0
    for aln in evidExon.alignFeatures:
        indelSize = getIndelSize(aln)
        if indelSize > exonPolymorhicSizeLimit:
            return EvidenceSupport.large_indel_size
        totalIndelSize += indelSize
    if totalIndelSize > exonPolymorhicFactionLimit * len(evidExon.chromLoc):
        return EvidenceSupport.large_indel_content
    else:
        return EvidenceSupport.polymorphic


def _compareMegExon(annotExon, evidExon):
    # boundaries are check with introns, so just check indels
    if len(evidExon.alignFeatures) > 1:
        return _checkMegExonIndels(evidExon)
    else:
        return EvidenceSupport.good


def _compareIntron(annotIntron, evidIntron):
    if not annotIntron.chromLoc.eqAbsLoc(evidIntron.chromLoc):
        return EvidenceSupport.exon_boundry_mismatch
    elif len(evidIntron.alignFeatures) > 1:
        return EvidenceSupport.internal_unaligned
    else:
        return EvidenceSupport.good


def _compareFeature(annotFeat, evidFeat):
    assert type(annotFeat) == type(evidFeat)
    if isinstance(annotFeat, ExonFeature):
        return _compareMegExon(annotFeat, evidFeat)
    else:
        return _compareIntron(annotFeat, evidFeat)


def compareMegWithEvidence(annotTrans, evidTrans):
    """Compare a multi-exon annotation with a given piece of evidence"""
    if len(evidTrans.features) != len(annotTrans.features):
        return EvidenceSupport.feat_mismatch
    worstCmpr = EvidenceSupport.good
    for iFeat in range(len(annotTrans.features)):
        cmpr = _compareFeature(annotTrans.features[iFeat], evidTrans.features[iFeat])
        if cmpr > worstCmpr:  # > is worse
            worstCmpr = cmpr
    return worstCmpr


class EvidenceCache(object):
    """Cache for a range of evidence, normally for one gene"""
    def __init__(self, evidenceReader, bounds):
        self.evidenceReader = evidenceReader
        self.bounds = bounds
        self.evidBySrc = {evidSrc: None for evidSrc in EvidenceSource}

    def get(self, evidSrc):
        "get a type of evidence"
        if self.evidBySrc[evidSrc] is None:
            overGen = self.evidenceReader.genOverlapping(evidSrc, self.bounds.name, self.bounds.start, self.bounds.end,
                                                         rnaStrand=self.bounds.strand, minExons=2)
            self.evidBySrc[evidSrc] = tuple(sorted(overGen, key=lambda a: (a.rnaLoc.name, a.chromLoc.name, a.chromLoc.start)))
        return self.evidBySrc[evidSrc]


class AnnotationEvidenceEval(namedtuple("AnnotationEvidence",
                                        ("annotId", "evidSrc", "evidId", "support", "suspect"))):
    """Evaluation of an annotation against an evidence alignment."""
    slots = ()


class AnnotationEvidenceCollector(defaultdict):
    """Collection of evidence supporting a transcript annotation,
    indexed by EvidenceSource"""
    def __init__(self, annotTrans):
        self.annotTrans = annotTrans
        super(AnnotationEvidenceCollector, self).__init__(list)

    def add(self, evidSrc, evidEval):
        self[evidSrc].append(evidEval)


def _countFullSupport(evidEvals):
    """count support from normal and suspect evidence"""
    supporting = [ev for ev in evidEvals if ev.support <= EvidenceSupport.polymorphic]
    return (len([ev for ev in supporting if ev.suspect is None]),
            len([ev for ev in supporting if ev.suspect is not None]))


def _calculateRnaTsl(evidCollector):
    """support from RNAs in a given set"""
    goodCnt, suspectCnt = _countFullSupport(evidCollector)
    if goodCnt >= 1:
        return TrascriptionSupportLevel.tsl1
    elif suspectCnt >= 1:
        return TrascriptionSupportLevel.tsl2
    else:
        return TrascriptionSupportLevel.tsl5


def _calculateEstTsl(evidCollector):
    """support from ESTs in a given set"""
    goodCnt, suspectCnt = _countFullSupport(evidCollector)
    if goodCnt >= 2:
        return TrascriptionSupportLevel.tsl2
    elif goodCnt == 1:
        return TrascriptionSupportLevel.tsl3
    elif suspectCnt >= 1:
        return TrascriptionSupportLevel.tsl4
    else:
        return TrascriptionSupportLevel.tsl5


def _calculateTsl(evidCollector):
    """compute TSL from evidence in evidCollector"""
    ucscRnaTsl = _calculateRnaTsl(evidCollector[EvidenceSource.UCSC_RNA])
    ensemblRnaTsl = _calculateRnaTsl(evidCollector[EvidenceSource.ENSEMBL_RNA])
    estRnaTsl = _calculateEstTsl(evidCollector[EvidenceSource.UCSC_EST])
    return min(ucscRnaTsl, ensemblRnaTsl, estRnaTsl)


def writeTsvHeaders(tslTsvFh, detailsTsvFh=None):
    fileOps.prRowv(tslTsvFh, "transcriptId", "level")
    if detailsTsvFh is not None:
        fileOps.prRowv(detailsTsvFh, "transcriptId", "evidSrc", "evidId", "evidSupport", "suspect")


def _writeDetails(detailsTsvFh, annotTrans, evidSrc, evidTrans, evidSupport, suspect):
    fileOps.prRowv(detailsTsvFh, annotTrans.rnaLoc.name, evidSrc, evidTrans.rnaLoc.name, evidSupport,
                   "" if suspect is None else suspect)


def _compareWithEvidence(annotTrans, evidSrc, evidTrans, evidCollector, detailsTsvFh):
    evidSupport = compareMegWithEvidence(annotTrans, evidTrans)
    suspect = evidTrans.attrs.genbankProblem
    if detailsTsvFh is not None:
        _writeDetails(detailsTsvFh, annotTrans, evidSrc, evidTrans, evidSupport, suspect)
    if evidSupport < EvidenceSupport.feat_mismatch:
        evidCollector.add(evidSrc,
                          AnnotationEvidenceEval(annotTrans, evidSrc, evidTrans.rnaLoc.name, evidSupport, suspect))


def _collectTransSupport(annotTrans, evidCache, detailsTsvFh):
    evidCollector = AnnotationEvidenceCollector(annotTrans)
    for evidSrc in EvidenceSource:
        for evidTrans in evidCache.get(evidSrc):
            _compareWithEvidence(annotTrans, evidSrc, evidTrans, evidCollector, detailsTsvFh)
    return evidCollector


def _classifyTrans(annotTrans, evidCache, tslTsvFh, detailsTsvFh):
    if _transIsSingleExon(annotTrans) or _geneIsTslIgnored(annotTrans):
        tsl = TrascriptionSupportLevel.tslNA
    else:
        evidCollector = _collectTransSupport(annotTrans, evidCache, detailsTsvFh)
        tsl = _calculateTsl(evidCollector)
    fileOps.prRowv(tslTsvFh, annotTrans.rnaLoc.name, tsl)


def classifyGeneTranscripts(evidenceReader, geneAnnotTranses, tslTsvFh, detailsTsvFh=None):
    """Classify a list of transcripts, which must be all on the same chromosome
    and should be from the same gene."""
    evidCache = EvidenceCache(evidenceReader, findAnnotationBounds(geneAnnotTranses))
    for annotTrans in geneAnnotTranses:
        _classifyTrans(annotTrans, evidCache, tslTsvFh, detailsTsvFh)
