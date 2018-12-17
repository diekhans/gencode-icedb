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

# FIXME: need to decide on best limits for various cases.

debug = False
# debug = True

# limits on size of as single indel in an exon.
tightExonPolymorphicSizeLimit = 12

# fraction of allowed total indel size relative exon length
tightExonPolymorphicFactionLimit = 0.1

# loose for comparability with old TSL code, which most ignored indels
looseExonPolymorphicSizeLimit = 5000
looseExonPolymorphicFactionLimit = 1.0


# FIXME these are from the ccds2/modules/gencode/src/lib/gencode/data/gencodeGenes.py, migrate to new module
def _transIsSingleExon(transAnnot):
    return len(transAnnot.getFeaturesOfType(ExonFeature)) <= 1


def _geneIsTslIgnored(transAnnot):
    bioType = transAnnot.attrs.transcriptType
    geneName = transAnnot.attrs.geneName
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


def sameChromBounds(feat1, feat2):
    """is the location on the chromosome identical, regardless of strand"""
    return feat1.chrom.eqAbsLoc(feat2.chrom)


class EvidenceQualityEval(object):
    """Evaluates quality of evidence aliment.

    Parameterized by attributes determining fuzziness of match.
    - exonPolymorphicSizeLimit - maximum size of a single indel in an exon.
    - exonPolymorphicFactionLimit - maximum fraction of allowed total indel
      size relative exon length.
    """
    def __init__(self, exonPolymorphicSizeLimit, exonPolymorphicFactionLimit):
        self.exonPolymorphicSizeLimit = exonPolymorphicSizeLimit
        self.exonPolymorphicFactionLimit = exonPolymorphicFactionLimit

    def _checkExonIndels(self, evidExon):
        "check for allowed indel polymorphism"

        def isInitialTerminalUnaligned(evidExon, aln):
            return (aln.rna.start == 0) or (aln.rna.end == evidExon.transcript.rna.size)

        def getIndelSize(aln):
            if isinstance(aln, ChromInsertFeature):
                return len(aln.chrom)
            elif isinstance(aln, RnaInsertFeature) and not isInitialTerminalUnaligned(evidExon, aln):
                return len(aln.rna)
            else:
                return 0  # not an indel

        if len(evidExon.alignFeatures) == 1:
            return EvidenceSupport.good

        totalIndelSize = 0
        for aln in evidExon.alignFeatures:
            indelSize = getIndelSize(aln)
            if indelSize > self.exonPolymorphicSizeLimit:
                return EvidenceSupport.large_indel_size
            totalIndelSize += indelSize

        if totalIndelSize == 0:
            return EvidenceSupport.good
        elif totalIndelSize > self.exonPolymorphicFactionLimit * len(evidExon.chrom):
            return EvidenceSupport.large_indel_content
        else:
            return EvidenceSupport.polymorphic

    def _checkIntronIndels(self, evidIntron):
        if len(evidIntron.alignFeatures) > 1:
            return EvidenceSupport.internal_unaligned
        else:
            return EvidenceSupport.good

    def _checkEvidFeatQuality(self, evidFeat):
        if isinstance(evidFeat, ExonFeature):
            return self._checkExonIndels(evidFeat)
        else:
            return self._checkIntronIndels(evidFeat)

    def check(self, evidTrans):
        """Initial validation of the quality of evidence, returning the best support it
        could provide."""
        worstSupport = EvidenceSupport.good
        for evidFeat in evidTrans.features:
            worstSupport = max(self._checkEvidFeatQuality(evidFeat), worstSupport)
        return worstSupport


class MegSupportEvaluator(object):
    """Evaluate a multi-exon annotation against an evidence alignment.  If
    allowExtension is specified, allow evidence to extend beyond annotation.
    qualEval is an instance of EvidenceQualityEval that defined the mention
    of the evaluation.
    """
    def __init__(self, qualEval, allowExtension=False):
        self.qualEval = qualEval
        self.allowExtension = allowExtension

    def _findEvidExonRangeStart(self, transAnnot, evidTrans):
        firstAnnotExon = transAnnot.features[0]
        firstEvidExon = evidTrans.features[0]
        while (firstEvidExon is not None) and (firstEvidExon.chrom.end < firstAnnotExon.chrom.start):
            firstEvidExon = firstEvidExon.nextFeature(ExonFeature)
        if (firstEvidExon is None) or (firstEvidExon.chrom.start > firstAnnotExon.chrom.end):
            return None  # not found or after firstAnnotExon
        return firstEvidExon  # overlapping first

    def _findEvidExonRangeEnd(self, transAnnot, evidTrans):
        lastAnnotExon = transAnnot.features[-1]
        lastEvidExon = evidTrans.features[-1]
        while (lastEvidExon is not None) and (lastEvidExon.chrom.start > lastAnnotExon.chrom.end):
            lastEvidExon = lastEvidExon.prevFeature(ExonFeature)
        if (lastEvidExon is None) or (lastEvidExon.chrom.end < lastAnnotExon.chrom.start):
            return None  # not found or before lastAnnotExon
        return lastEvidExon  # overlapping last

    def _findEvidExonRange(self, transAnnot, evidTrans):
        """Find the first and last evidence exons that overlap the first and last
        exons of the transcript.  Return None,None when there isn't an overlap.
        This range is what is compared."""
        # Walk from start to find start and end to find end, as multiple might overlap,
        # however don't go past bounds of the annotation, in case evidence and annotation
        # are interleaved.  This would be simper without the extend mode.
        return (self._findEvidExonRangeStart(transAnnot, evidTrans),
                self._findEvidExonRangeEnd(transAnnot, evidTrans))

    def _compareFeatureCounts(self, transAnnot, evidTrans):
        """fast check if number of features warrant detailed checking"""
        if self.allowExtension:
            if len(evidTrans.features) < len(transAnnot.features):
                return EvidenceSupport.feat_count_mismatch
        else:
            if len(evidTrans.features) != len(transAnnot.features):
                return EvidenceSupport.feat_count_mismatch
        return EvidenceSupport.good

    def _compareIntron(self, annotIntron, evidIntron):
        if not sameChromBounds(annotIntron, evidIntron):
            return EvidenceSupport.exon_boundry_mismatch
        else:
            return EvidenceSupport.good

    def _compareExonStart(self, annotExon, evidExon):
        if evidExon.chrom.start != annotExon.chrom.start:
            # start does not match, is this an error?
            if annotExon.prevFeature() is not None:
                return EvidenceSupport.feat_mismatch  # not start of annotation
            elif evidExon.prevFeature() is not None:
                return EvidenceSupport.feat_mismatch  # evidence extends beyond
            else:
                return EvidenceSupport.good
        elif (annotExon.prevFeature() is None) and (evidExon.prevFeature() is not None) and (not self.allowExtension):
            # start matches, but there is unallowed preceding evidence
            return EvidenceSupport.feat_mismatch
        else:
            return EvidenceSupport.good

    def _compareExonEnd(self, annotExon, evidExon):
        if evidExon.chrom.end != annotExon.chrom.end:
            # end does not match, is this an error?
            if annotExon.nextFeature() is not None:
                return EvidenceSupport.feat_mismatch  # not end of annotation
            elif evidExon.nextFeature() is not None:
                return EvidenceSupport.feat_mismatch  # evidence extends beyond
            else:
                return EvidenceSupport.good
        elif (annotExon.nextFeature() is None) and (evidExon.nextFeature() is not None) and (not self.allowExtension):
            # start matches, but there is unallowed preceding evidence
            return EvidenceSupport.feat_mismatch
        else:
            return EvidenceSupport.good

    def _compareExon(self, annotExon, evidExon):
        # Bounds must match exactly except for outside ends of initial or terminal
        # annotation exon.  In extension mode, if there are evidence features
        # beyond the ends, then both those must match.  In non-extension mode,
        # there will not be features outside of the annotation, so the check is the same
        if not evidExon.chrom.overlaps(annotExon.chrom):
            return EvidenceSupport.feat_mismatch  # don't even overlap
        return max(self._compareExonStart(annotExon, evidExon),
                   self._compareExonEnd(annotExon, evidExon))

    def _compareFeature(self, annotFeat, evidFeat):
        # evidence has already been validated for indels.
        if type(annotFeat) != type(evidFeat):
            return EvidenceSupport.feat_mismatch
        if isinstance(annotFeat, IntronFeature):
            return self._compareIntron(annotFeat, evidFeat)
        elif isinstance(annotFeat, ExonFeature):
            return self._compareExon(annotFeat, evidFeat)
        else:
            raise Exception("BUG: unexpected structural feature type: {}", type(annotFeat))

    def _compareFeatures(self, transAnnot, firstEvidExon, lastEvidExon):
        worstSupport = EvidenceSupport.good
        annotFeat = transAnnot.firstFeature()
        evidFeat = firstEvidExon
        while worstSupport < EvidenceSupport.poor:
            worstSupport = max(self._compareFeature(annotFeat, evidFeat), worstSupport)
            if evidFeat is lastEvidExon:
                break
            evidFeat = evidFeat.nextFeature()
            annotFeat = annotFeat.nextFeature()
            if annotFeat is None:
                worstSupport = max(EvidenceSupport.feat_mismatch, worstSupport)
                break  # mismatch if features are not out of sync
        return worstSupport

    def _compareMegWithEvidenceImpl(self, transAnnot, evidTrans):
        """Compare a multi-exon annotation with a given piece of evidence"""
        # check full evidence first; > is worse
        worstSupport = self.qualEval.check(evidTrans)
        if worstSupport >= EvidenceSupport.poor:
            return worstSupport  # worse than allowed, give up now

        # fast path check
        worstSupport = max(self._compareFeatureCounts(transAnnot, evidTrans), worstSupport)
        if worstSupport >= EvidenceSupport.poor:
            return worstSupport
        # get range of exons to compare
        firstEvidExon, lastEvidExon = self._findEvidExonRange(transAnnot, evidTrans)
        if (firstEvidExon is None) or (lastEvidExon is None):
            return EvidenceSupport.feat_mismatch

        worstSupport = max(self._compareFeatures(transAnnot, firstEvidExon, lastEvidExon), worstSupport)
        return worstSupport

    def compare(self, transAnnot, evidTrans):
        """Compare a multi-exon annotation with a given piece of evidence,"""
        try:
            if debug:
                transAnnot.dump(msg="annotation")
                evidTrans.dump(msg="evidence")
            return self._compareMegWithEvidenceImpl(transAnnot, evidTrans)
        except Exception as ex:
            raise Exception("Bug evaluating {} with {}".format(transAnnot.rna.name, evidTrans.rna.name)) from ex


class EvidenceCache(object):
    """Cache for a range of evidence, normally for one gene"""
    def __init__(self, evidenceReader, bounds, transcriptionStrand):
        self.evidenceReader = evidenceReader
        self.bounds = bounds
        self.transcriptionStrand = transcriptionStrand
        self.evidBySrc = {evidSrc: None for evidSrc in evidenceReader.sources}  # initialize to None

    @property
    def sources(self):
        return self.evidenceReader.sources

    def _load(self, evidSrc):
        overGen = self.evidenceReader.genOverlapping(evidSrc, self.bounds.name, self.bounds.start, self.bounds.end,
                                                     transcriptionStrand=self.transcriptionStrand, minExons=2)
        self.evidBySrc[evidSrc] = tuple(sorted(overGen, key=lambda a: (a.rna.name, a.chrom.name, a.chrom.start)))

    def get(self, evidSrc):
        "get a type of evidence"
        assert evidSrc in self.sources
        if self.evidBySrc[evidSrc] is None:
            self._load(evidSrc)
        return self.evidBySrc[evidSrc]


class AnnotationEvidenceEval(namedtuple("AnnotationEvidenceEval",
                                        ("annotId", "evidSrc", "evidId", "support", "suspect"))):
    """Evaluation of an annotation against an evidence alignment."""
    slots = ()


class AnnotationEvidenceCollector(defaultdict):
    """Collection of evidence supporting a transcript annotation,
    indexed by EvidenceSource"""
    def __init__(self, transAnnot):
        self.transAnnot = transAnnot
        super(AnnotationEvidenceCollector, self).__init__(list)

    def add(self, evidSrc, evidEval):
        self[evidSrc].append(evidEval)


class FullLengthSupportEvaluator(object):
    """
    Full-length support evaluation.
    qualEval is an instance of EvidenceQualityEval that defined the method
    of the evaluation.
    allowExtension indicates if new exons should be allowed
    """
    def __init__(self, evidenceReader, qualEval, allowExtension=False):
        self.evidenceReader = evidenceReader
        self.evaluator = MegSupportEvaluator(qualEval, allowExtension)

    def _countFullSupport(self, evidEvals):
        """count support from normal and suspect evidence"""
        supporting = [ev for ev in evidEvals if ev.support < EvidenceSupport.poor]
        return (len([ev for ev in supporting if ev.suspect is None]),
                len([ev for ev in supporting if ev.suspect is not None]))

    def _calculateRnaTsl(self, evidCollector):
        """support from RNAs in a given set"""
        goodCnt, suspectCnt = self._countFullSupport(evidCollector)
        if goodCnt >= 1:
            return TrascriptionSupportLevel.tsl1
        elif suspectCnt >= 1:
            return TrascriptionSupportLevel.tsl2
        else:
            return TrascriptionSupportLevel.tsl5

    def _calculateEstTsl(self, evidCollector):
        """support from ESTs in a given set"""
        goodCnt, suspectCnt = self._countFullSupport(evidCollector)
        if goodCnt >= 2:
            return TrascriptionSupportLevel.tsl2
        elif goodCnt == 1:
            return TrascriptionSupportLevel.tsl3
        elif suspectCnt >= 1:
            return TrascriptionSupportLevel.tsl4
        else:
            return TrascriptionSupportLevel.tsl5

    def _calculateTsl(self, evidCollector):
        """compute TSL from evidence in evidCollector"""
        ucscRnaTsl = self._calculateRnaTsl(evidCollector[EvidenceSource.UCSC_RNA])
        ensemblRnaTsl = self._calculateRnaTsl(evidCollector[EvidenceSource.ENSEMBL_RNA])
        estRnaTsl = self._calculateEstTsl(evidCollector[EvidenceSource.UCSC_EST])
        return min(ucscRnaTsl, ensemblRnaTsl, estRnaTsl)

    def _writeDetails(self, detailsTsvFh, transAnnot, evidSrc, evidTrans, evidSupport, suspect):
        fileOps.prRowv(detailsTsvFh, transAnnot.rna.name, evidSrc, evidTrans.rna.name, evidSupport,
                       "" if suspect is None else suspect)

    def _compareWithEvidence(self, transAnnot, evidSrc, evidTrans, evidCollector, detailsTsvFh):
        evidSupport = self.evaluator.compare(transAnnot, evidTrans)
        suspect = evidTrans.attrs.genbankProblem
        if detailsTsvFh is not None:
            self._writeDetails(detailsTsvFh, transAnnot, evidSrc, evidTrans, evidSupport, suspect)
        if evidSupport < EvidenceSupport.poor:
            evidCollector.add(evidSrc,
                              AnnotationEvidenceEval(transAnnot, evidSrc, evidTrans.rna.name, evidSupport, suspect))

    def collectTransSupport(self, transAnnot, evidCache, detailsTsvFh=None):
        """lower-level function to collect evidence for one transcript,
        returning an AnnotationEvidenceCollector object"""
        evidCollector = AnnotationEvidenceCollector(transAnnot)
        for evidSrc in evidCache.sources:
            for evidTrans in evidCache.get(evidSrc):
                self._compareWithEvidence(transAnnot, evidSrc, evidTrans, evidCollector, detailsTsvFh)
        return evidCollector

    def _classifyTrans(self, transAnnot, evidCache, tslTsvFh, detailsTsvFh):
        if _transIsSingleExon(transAnnot) or _geneIsTslIgnored(transAnnot):
            tsl = TrascriptionSupportLevel.tslNA
        else:
            evidCollector = self.collectTransSupport(transAnnot, evidCache, detailsTsvFh)
            tsl = self._calculateTsl(evidCollector)
        fileOps.prRowv(tslTsvFh, transAnnot.rna.name, tsl)

    @staticmethod
    def writeTsvHeaders(tslTsvFh, detailsTsvFh=None):
        fileOps.prRowv(tslTsvFh, "transcriptId", "level")
        if detailsTsvFh is not None:
            fileOps.prRowv(detailsTsvFh, "transcriptId", "evidSrc", "evidId", "evidSupport", "suspect")

    def classifyGeneTranscripts(self, geneAnnot, tslTsvFh, detailsTsvFh=None):
        """Classify a list of transcripts, which must be all on the same chromosome.
        They should be from the same gene locus or overlapping loci for caching efficiency."""
        evidCache = EvidenceCache(self.evidenceReader, geneAnnot.chrom, geneAnnot.transcriptionStrand)
        for transAnnot in geneAnnot.transcripts:
            self._classifyTrans(transAnnot, evidCache, tslTsvFh, detailsTsvFh)
