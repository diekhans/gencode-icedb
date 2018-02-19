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
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, ChromInsertFeature, RnaInsertFeature
from gencode_icedb.tsl.supportDefs import EvidenceSupport, TrascriptionSupportLevel
from gencode_icedb.tsl.evidenceDb import EvidenceSource
from gencode_icedb.general.gencodeDb import findAnnotationBounds


# limits on size of as single indel in an exon.
exonPolymorhicSizeLimit = 12

# fraction of allowed total indel size relative exon length
exonPolymorhicFactionLimit = 0.05


# FIXME these are from the ccds2/modules/gencode/src/lib/gencode/data/gencodeGenes.py, migrate to new module
def _transIsSingleExon(annotTrans):
    return len(annotTrans.getFeaturesOfType(ExonFeature)) <= 1


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


def _checkExonIndels(evidExon):
    "check for allowed indel polymorphism"

    def getIndelSize(aln):
        if isinstance(aln, ChromInsertFeature):
            return len(aln.chrom)
        elif isinstance(aln, RnaInsertFeature):
            return len(aln.rna)
        else:
            return 0  # not an indel

    if len(evidExon.alignFeatures) == 1:
        return EvidenceSupport.good
    totalIndelSize = 0
    for aln in evidExon.alignFeatures:
        indelSize = getIndelSize(aln)
        if indelSize > exonPolymorhicSizeLimit:
            return EvidenceSupport.large_indel_size
        totalIndelSize += indelSize
    if totalIndelSize > exonPolymorhicFactionLimit * len(evidExon.chrom):
        return EvidenceSupport.large_indel_content
    else:
        return EvidenceSupport.polymorphic


def _checkIntronIndels(evidIntron):
    if len(evidIntron.alignFeatures) > 1:
        return EvidenceSupport.internal_unaligned
    else:
        return EvidenceSupport.good


def _checkEvidFeatQuality(evidFeat):
    if isinstance(evidFeat, ExonFeature):
        return _checkExonIndels(evidFeat)
    else:
        return _checkIntronIndels(evidFeat)


def checkEvidQuality(evidTrans):
    "initial validation of the quality of evidence"
    worstSupport = EvidenceSupport.good
    for evidFeat in evidTrans.features:
        worstSupport = max(_checkEvidFeatQuality(evidFeat), worstSupport)
    return worstSupport


def findEvidExonRangeStart(annotTrans, evidTrans):
    firstAnnotExon = annotTrans.features[0]
    firstEvidExon = evidTrans.features[0]
    while (firstEvidExon is not None) and (firstEvidExon.chrom.end < firstAnnotExon.chrom.start):
        firstEvidExon = firstEvidExon.nextFeature(ExonFeature)
    if (firstEvidExon is None) or (firstEvidExon.chrom.start > firstAnnotExon.chrom.end):
        return None  # not found or after firstAnnotExon
    return firstEvidExon  # overlapping first


def findEvidExonRangeEnd(annotTrans, evidTrans):
    lastAnnotExon = annotTrans.features[-1]
    lastEvidExon = evidTrans.features[-1]
    while (lastEvidExon is not None) and (lastEvidExon.chrom.start > lastAnnotExon.chrom.end):
        lastEvidExon = lastEvidExon.prevFeature(ExonFeature)
    if (lastEvidExon is None) or (lastEvidExon.chrom.end < lastAnnotExon.chrom.start):
        return None  # not found or before lastAnnotExon
    return lastEvidExon  # overlapping last


def findEvidExonRange(annotTrans, evidTrans):
    """Find the first and last evidence exons that overlap the first and last
    exons of the transcript.  Return None,None when there isn't an overlap.
    This range is what is compared."""
    # Walk from start to find start and end to find end, as multiple might overlap,
    # however don't go past bounds of the annotation, in case evidence and annotation
    # are interleaved.  This would be simper without the extend mode.
    return (findEvidExonRangeStart(annotTrans, evidTrans),
            findEvidExonRangeEnd(annotTrans, evidTrans))


def _compareFeatureCounts(annotTrans, evidTrans, allowExtension):
    """fast check if number of features warrant detailed checking"""
    if allowExtension:
        if len(evidTrans.features) < len(annotTrans.features):
            return EvidenceSupport.feat_count_mismatch
    else:
        if len(evidTrans.features) != len(annotTrans.features):
            return EvidenceSupport.feat_count_mismatch
    return EvidenceSupport.good


def sameChromBounds(feat1, feat2):
    return feat1.chrom.eqAbsLoc(feat2.chrom)


def _compareIntron(annotIntron, evidIntron):
    if not sameChromBounds(annotIntron, evidIntron):
        return EvidenceSupport.exon_boundry_mismatch
    else:
        return EvidenceSupport.good


def _compareExonStart(annotExon, evidExon, allowExtension):
    if evidExon.chrom.start != annotExon.chrom.start:
        # start does not match, is this an error?
        if annotExon.prevFeature() is not None:
            return EvidenceSupport.feat_mismatch  # not start of annotation
        elif evidExon.prevFeature() is not None:
            return EvidenceSupport.feat_mismatch  # evidence extends beyond
        else:
            return EvidenceSupport.good
    elif (annotExon.prevFeature() is None) and (evidExon.prevFeature() is not None) and (not allowExtension):
        # start matches, but there is unallowed preceding evidence
        return EvidenceSupport.feat_mismatch
    else:
        return EvidenceSupport.good


def _compareExonEnd(annotExon, evidExon, allowExtension):
    if evidExon.chrom.end != annotExon.chrom.end:
        # end does not match, is this an error?
        if annotExon.nextFeature() is not None:
            return EvidenceSupport.feat_mismatch  # not end of annotation
        elif evidExon.nextFeature() is not None:
            return EvidenceSupport.feat_mismatch  # evidence extends beyond
        else:
            return EvidenceSupport.good
    elif (annotExon.nextFeature() is None) and (evidExon.nextFeature() is not None) and (not allowExtension):
        # start matches, but there is unallowed preceding evidence
        return EvidenceSupport.feat_mismatch
    else:
        return EvidenceSupport.good


def _compareExon(annotExon, evidExon, allowExtension):
    # Bounds must match exactly except for outside ends of initial or terminal
    # annotation exon.  In extension mode, if there are evidence features
    # beyond the ends, then both those must match.  In non-extension mode,
    # there will not be features outside of the annotation, so the check is the same
    if not evidExon.chrom.overlaps(annotExon.chrom):
        return EvidenceSupport.feat_mismatch  # don't even overlap
    return max(_compareExonStart(annotExon, evidExon, allowExtension),
               _compareExonEnd(annotExon, evidExon, allowExtension))


def _compareFeature(annotFeat, evidFeat, allowExtension):
    # evidence has already been validated for indels.
    if type(annotFeat) != type(evidFeat):
        return EvidenceSupport.feat_mismatch
    if isinstance(annotFeat, IntronFeature):
        return _compareIntron(annotFeat, evidFeat)
    elif isinstance(annotFeat, ExonFeature):
        return _compareExon(annotFeat, evidFeat, allowExtension)
    else:
        raise Exception("BUG: unexpected structural feature type: {}", type(annotFeat))


def _compareFeatures(annotTrans, firstEvidExon, lastEvidExon, allowExtension):
    worstSupport = EvidenceSupport.good
    annotFeat = annotTrans.firstFeature()
    evidFeat = firstEvidExon
    while worstSupport < EvidenceSupport.poor:
        worstSupport = max(_compareFeature(annotFeat, evidFeat, allowExtension),
                           worstSupport)
        if evidFeat is lastEvidExon:
            break
        evidFeat = evidFeat.nextFeature()
        annotFeat = annotFeat.nextFeature()
        if annotFeat is None:
            worstSupport = max(EvidenceSupport.feat_mismatch, worstSupport)
            break  # mismatch if features are not out of sync
    return worstSupport


def _compareMegWithEvidenceImpl(annotTrans, evidTrans, allowExtension=False):
    """Compare a multi-exon annotation with a given piece of evidence"""
    # check full evidence first; > is worse
    worstSupport = checkEvidQuality(evidTrans)
    if worstSupport >= EvidenceSupport.poor:
        return worstSupport
    # fast path check
    worstSupport = max(_compareFeatureCounts(annotTrans, evidTrans, allowExtension),
                       worstSupport)
    if worstSupport >= EvidenceSupport.poor:
        return worstSupport
    # get range of exons to compare
    firstEvidExon, lastEvidExon = findEvidExonRange(annotTrans, evidTrans)
    if (firstEvidExon is None) or (lastEvidExon is None):
        return EvidenceSupport.feat_mismatch

    worstSupport = max(_compareFeatures(annotTrans, firstEvidExon, lastEvidExon, allowExtension),
                       worstSupport)
    return worstSupport


def compareMegWithEvidence(annotTrans, evidTrans, allowExtension=False):
    """Compare a multi-exon annotation with a given piece of evidence, If
    allowExtension is true, evidence features at allowed beyond the start and
    end of the annotation.  This does not check txStart/txEnd is consistent
    with evidence, so this can use to find extensions of exons."""
    try:
        return _compareMegWithEvidenceImpl(annotTrans, evidTrans, allowExtension=allowExtension)
    except Exception as ex:
        raise Exception("Bug evaluating {} with {}".format(annotTrans.rna.name, evidTrans.rna.name)) from ex


class EvidenceCache(object):
    """Cache for a range of evidence, normally for one gene"""
    def __init__(self, evidenceReader, bounds):
        self.evidenceReader = evidenceReader
        self.bounds = bounds
        self.evidBySrc = {evidSrc: None for evidSrc in evidenceReader.sources}

    @property
    def sources(self):
        return self.evidenceReader.sources

    def _load(self, evidSrc):
        overGen = self.evidenceReader.genOverlapping(evidSrc, self.bounds.name, self.bounds.start, self.bounds.end,
                                                     rnaStrand=self.bounds.strand, minExons=2)
        self.evidBySrc[evidSrc] = tuple(sorted(overGen, key=lambda a: (a.rna.name, a.chrom.name, a.chrom.start)))

    def get(self, evidSrc):
        "get a type of evidence"
        assert evidSrc in self.sources
        if self.evidBySrc[evidSrc] is None:
            self._load(evidSrc)
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
    supporting = [ev for ev in evidEvals if ev.support < EvidenceSupport.poor]
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
    fileOps.prRowv(detailsTsvFh, annotTrans.rna.name, evidSrc, evidTrans.rna.name, evidSupport,
                   "" if suspect is None else suspect)


def _compareWithEvidence(annotTrans, evidSrc, evidTrans, evidCollector, detailsTsvFh):
    evidSupport = compareMegWithEvidence(annotTrans, evidTrans)
    suspect = evidTrans.attrs.genbankProblem
    if detailsTsvFh is not None:
        _writeDetails(detailsTsvFh, annotTrans, evidSrc, evidTrans, evidSupport, suspect)
    if evidSupport < EvidenceSupport.poor:
        evidCollector.add(evidSrc,
                          AnnotationEvidenceEval(annotTrans, evidSrc, evidTrans.rna.name, evidSupport, suspect))


def _collectTransSupport(annotTrans, evidCache, detailsTsvFh):
    evidCollector = AnnotationEvidenceCollector(annotTrans)
    for evidSrc in evidCache.sources:
        for evidTrans in evidCache.get(evidSrc):
            _compareWithEvidence(annotTrans, evidSrc, evidTrans, evidCollector, detailsTsvFh)
    return evidCollector


def _classifyTrans(annotTrans, evidCache, tslTsvFh, detailsTsvFh):
    if _transIsSingleExon(annotTrans) or _geneIsTslIgnored(annotTrans):
        tsl = TrascriptionSupportLevel.tslNA
    else:
        evidCollector = _collectTransSupport(annotTrans, evidCache, detailsTsvFh)
        tsl = _calculateTsl(evidCollector)
    fileOps.prRowv(tslTsvFh, annotTrans.rna.name, tsl)


def classifyGeneTranscripts(evidenceReader, geneAnnotTranses, tslTsvFh, detailsTsvFh=None):
    """Classify a list of transcripts, which must be all on the same chromosome
    and should be from the same gene."""
    evidCache = EvidenceCache(evidenceReader, findAnnotationBounds(geneAnnotTranses))
    for annotTrans in geneAnnotTranses:
        _classifyTrans(annotTrans, evidCache, tslTsvFh, detailsTsvFh)
