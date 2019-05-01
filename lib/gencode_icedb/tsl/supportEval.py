"""
Evaluation of transcripts against evidence.
"""
from pycbio.sys import fileOps
from collections import defaultdict
from gencode_icedb.general.transFeatures import ExonFeature, ChromInsertFeature, RnaInsertFeature
from gencode_icedb.tsl.supportDefs import EvidenceSupport
from gencode_icedb.tsl.supportDefs import isGeneIgnored, isTransIgnored
from gencode_icedb.tsl.supportEvalDb import SupportEvidEvalResult, SupportEvalResult


# FIXME: inconsistent names:  transAnnot, transExon, evidTrans, evidExon, geneAnnot

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

# Categories we keep
keepEvidenceSupport = (EvidenceSupport.good, EvidenceSupport.polymorphic, EvidenceSupport.extends_exons)


def sameChromBounds(feat1, feat2):
    """is the location on the chromosome identical, regardless of strand"""
    return feat1.chrom.eqAbsLoc(feat2.chrom)


def keepEvidEval(support):
    """good or shows extensions"""
    return ((support < EvidenceSupport.poor)
            or (support == EvidenceSupport.extends_exons))


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
    def __init__(self, evidSetUuid, qualEval, allowExtension=False):
        self.evidSetUuid = evidSetUuid
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
        This range is what is compared,"""
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

    def _countExtraExons(self, dir, evidExon):
        cnt = 0
        while True:
            evidExon = evidExon.prevFeature(ExonFeature) if dir < 0 else evidExon.nextFeature(ExonFeature)
            if evidExon is None:
                break
            cnt += 1
        return cnt

    def _compareFirstExon(self, annotExon, evidExon):
        "return (support, offset, extraExons)"
        if (evidExon.chrom.end != annotExon.chrom.end):
            # must support end of exon
            return (EvidenceSupport.feat_mismatch, 0, 0)
        elif evidExon.prevFeature() is not None:
            # evidence has exons before start, must be allowed and exon start match
            if (not self.allowExtension) or (evidExon.chrom.start != annotExon.chrom.start):
                return (EvidenceSupport.feat_mismatch, 0, 0)
            else:
                return (EvidenceSupport.extends_exons, 0, self._countExtraExons(-1, evidExon))
        else:
            return (EvidenceSupport.good, annotExon.chrom.start - evidExon.chrom.start, 0)

    def _compareLastExon(self, annotExon, evidExon):
        "returns (support, offset, extraExons)"
        if (evidExon.chrom.start != annotExon.chrom.start):
            # must support start of exon
            return (EvidenceSupport.feat_mismatch, 0, 0)
        elif evidExon.nextFeature() is not None:
            # evidence has exons after end, must be allowed and exon end match
            if (not self.allowExtension) or (evidExon.chrom.end != annotExon.chrom.end):
                return (EvidenceSupport.feat_mismatch, 0, 0)
            else:
                return (EvidenceSupport.extends_exons, 0, self._countExtraExons(1, evidExon))
        else:
            return (EvidenceSupport.good, evidExon.chrom.end - annotExon.chrom.end, 0)

    def _compareInternalExon(self, annotExon, evidExon):
        if (evidExon.chrom.start != annotExon.chrom.start) or (evidExon.chrom.end != annotExon.chrom.end):
            return EvidenceSupport.feat_mismatch
        else:
            return EvidenceSupport.good

    def _mkSupportEvidEvalResults(self, transAnnot, evidTrans, support,
                                  firstOffset=0, lastOffset=0, firstExtend=0, lastExtend=0):
        if transAnnot.transcriptionStrand == '+':
            offset5, offset3, extend5Exons, extend3Exons = firstOffset, lastOffset, firstExtend, lastExtend
        else:
            offset5, offset3, extend5Exons, extend3Exons = lastOffset, firstOffset, lastExtend, firstExtend
        return SupportEvidEvalResult(transAnnot.rna.name, self.evidSetUuid, evidTrans.rna.name, support,
                                     offset5, offset3, extend5Exons, extend3Exons)

    def _compareFeatures(self, transAnnot, evidTrans, firstEvidExon, lastEvidExon):
        # This assumes that the number of evidence exons has been checked to be
        # equal or more in cases of allowExtension.  This is not done in transcription
        # order.  This will never work on single-exon transcripts, which need a different
        # algorithm anyway
        # FIXME: is lastEvidExon needed?
        annotExon = transAnnot.firstFeature(ExonFeature)
        annotLastExon = transAnnot.lastFeature(ExonFeature)
        evidExon = firstEvidExon

        # initial exon
        firstSupport, firstOffset, firstExtend = self._compareFirstExon(annotExon, evidExon)
        worstSupport = firstSupport
        if not keepEvidEval(worstSupport):
            return self._mkSupportEvidEvalResults(transAnnot, evidTrans, worstSupport)

        # internal exons
        annotExon = annotExon.nextFeature(ExonFeature)
        evidExon = evidExon.nextFeature(ExonFeature)
        while annotExon != annotLastExon:
            worstSupport = max(self._compareInternalExon(annotExon, evidExon), worstSupport)
            if not keepEvidEval(worstSupport):
                return self._mkSupportEvidEvalResults(transAnnot, evidTrans, worstSupport)
            annotExon = annotExon.nextFeature(ExonFeature)
            evidExon = evidExon.nextFeature(ExonFeature)

        # final exon
        lastSupport, lastOffset, lastExtend = self._compareLastExon(annotExon, evidExon)
        worstSupport = max(lastSupport, worstSupport)
        if not keepEvidEval(worstSupport):
            return self._mkSupportEvidEvalResults(transAnnot, evidTrans, worstSupport)
        return self._mkSupportEvidEvalResults(transAnnot, evidTrans, worstSupport, firstOffset, lastOffset, firstExtend, lastExtend)

    def _compareMegWithEvidenceImpl(self, transAnnot, evidTrans):
        """Compare a multi-exon annotation with a given piece of evidence"""
        # check full evidence first; > is worse
        worstSupport = self.qualEval.check(evidTrans)
        if not keepEvidEval(worstSupport):
            return self._mkSupportEvidEvalResults(transAnnot, evidTrans, worstSupport)

        # fast path check
        worstSupport = max(self._compareFeatureCounts(transAnnot, evidTrans), worstSupport)
        if not keepEvidEval(worstSupport):
            return self._mkSupportEvidEvalResults(transAnnot, evidTrans, worstSupport)

        # get range of support exons to compare
        firstEvidExon, lastEvidExon = self._findEvidExonRange(transAnnot, evidTrans)
        if (firstEvidExon is None) or (lastEvidExon is None):
            return self._mkSupportEvidEvalResults(transAnnot, evidTrans, EvidenceSupport.feat_mismatch)

        return self._compareFeatures(transAnnot, evidTrans, firstEvidExon, lastEvidExon)

    def compare(self, transAnnot, evidTrans):
        """Compare a multi-exon annotation with a given piece of evidence, return a
        SupportEvidEvalResult object."""
        try:
            if debug:
                transAnnot.dump(msg="annotation")
                evidTrans.dump(msg="evidence")
            return self._compareMegWithEvidenceImpl(transAnnot, evidTrans)
        except Exception as ex:
            raise Exception("Bug evaluating {} with {}".format(transAnnot.rna.name, evidTrans.rna.name)) from ex


class FullLengthSupportEvaluator(object):
    """
    Full-length support evaluation.
    :param qualEval: is an instance of EvidenceQualityEval that defined the method
           of the evaluation.
    :param allowExtension: indicates if new exons should be allowed
    """
    def __init__(self, evidenceReader, qualEval, allowExtension=False):
        self.evidenceReader = evidenceReader
        self.evaluator = MegSupportEvaluator(evidenceReader.evidSetUuid,
                                             qualEval, allowExtension)

    def _compareWithEvidence(self, transAnnot, evidTrans, detailsTsvFh):
        evidEvalResult = self.evaluator.compare(transAnnot, evidTrans)
        if detailsTsvFh is not None:
            fileOps.prRow(detailsTsvFh, evidEvalResult.toRow())
        if keepEvidEval(evidEvalResult.support):
            return evidEvalResult
        else:
            return None

    def _collectSupportEvid(self, transAnnot, evidCache, detailsTsvFh):
        for evidTrans in evidCache:
            evr = self._compareWithEvidence(transAnnot, evidTrans, detailsTsvFh)
            if evr is not None:
                yield evr

    def _betterAlign(self, evr, evr0):
        if ((evr0 is None)
            or (evr.support < evr0.support)
            or ((evr.extend5Exons + evr.extend3Exons) > (evr0.extend5Exons + evr0.extend3Exons))
            or ((evr.offset5 + evr.offset3) > (evr0.offset5 + evr0.offset3))):
            return evr
        else:
            return evr0

    def collectSupportEvid(self, transAnnot, evidCache, detailsTsvFh=None):
        """lower-level function to collect all evidence for one transcript,
        returning an list of SupportEvidEvalResult objects. If there are
        multiple occurrences of the same id, pick the best.  This happens when
        RNAs are combined from two sources."""
        evrById = {}
        if not isTransIgnored(transAnnot):
            for result in self._collectSupportEvid(transAnnot, evidCache, detailsTsvFh):
                evrById[result.evidId] = self._betterAlign(result, evrById.get(result.evidId))
        return list(sorted(evrById.values(), key=lambda r: r.evidId))

    @staticmethod
    def writeTsvHeaders(supportEvalTsvFh, detailsTsvFh=None):
        fileOps.prRow(supportEvalTsvFh, SupportEvalResult.tsvHeader())
        if detailsTsvFh is not None:
            fileOps.prRow(detailsTsvFh, SupportEvidEvalResult.tsvHeader())

    def getEvidenceCache(self, annot):
        """get array of all evidence overlapping a gene or transcript annotation"""
        return tuple(self.evidenceReader.genOverlapping(annot.chrom, transcriptionStrand=annot.transcriptionStrand, minExons=2))

    def _combineSupportEvid(self, evidEvalResults):
        """combine SupportEvidEvalResult objects into a SupportEvalResult object."""
        evr = evidEvalResults[0]
        bestSupport, bestOffset5, bestOffset3, bestExtend5Exon, bestExtend3Exon = evr.support.value, evr.offset5, evr.offset3, evr.extend5Exons, evr.extend3Exons
        cnt = 1
        for evr in evidEvalResults[1:]:
            bestSupport = min(bestSupport, evr.support)
            bestOffset5 = max(bestOffset5, evr.offset5)
            bestOffset3 = max(bestOffset3, evr.offset3)
            bestExtend5Exon = max(bestExtend5Exon, evr.extend5Exons)
            bestExtend3Exon = max(bestExtend3Exon, evr.extend3Exons)
            cnt += 1
        return SupportEvalResult(evr.transcriptId, self.evidenceReader.evidSetUuid,
                                 EvidenceSupport(bestSupport), cnt,
                                 bestOffset5, bestOffset3,
                                 bestExtend5Exon, bestExtend3Exon)

    def _splitBySupport(self, evidEvalResults):
        evidEvalBySupport = defaultdict(list)
        for evr in evidEvalResults:
            evidEvalBySupport[evr.support].append(evr)
        return evidEvalBySupport

    def _mergeSupportEvidResults(self, evidEvalResults):
        """Merge results into the various categories we are going to keep.
        Results are kept for both supporting and ones suggesting modification,
        such as extension. """
        evidEvalBySupport = self._splitBySupport(evidEvalResults)
        evalResults = []
        for evidSupport in keepEvidenceSupport:
            if evidSupport in evidEvalBySupport:
                evalResults.append(self._combineSupportEvid(evidEvalBySupport[evidSupport]))
        # sort for reproducible
        return sorted(evalResults)

    def _evaulateTrans(self, transAnnot, evidCache, supportEvalTsvFh, detailsTsvFh):
        """get best supporting evidence and write to file, along with a
        count."""
        evidEvalResults = self.collectSupportEvid(transAnnot, evidCache, detailsTsvFh)
        for evalResult in self._mergeSupportEvidResults(evidEvalResults):
            fileOps.prRow(supportEvalTsvFh, evalResult.toRow())

    def evaluateGeneTranscripts(self, geneAnnot, supportEvalTsvFh, detailsTsvFh=None):
        """Evaluate a list of transcripts, which must be all on the same chromosome.
        They should be from the same gene locus or overlapping loci for caching efficiency."""
        if not isGeneIgnored(geneAnnot):
            evidCache = self.getEvidenceCache(geneAnnot)
            for transAnnot in geneAnnot.transcripts:
                self._evaulateTrans(transAnnot, evidCache, supportEvalTsvFh, detailsTsvFh)
