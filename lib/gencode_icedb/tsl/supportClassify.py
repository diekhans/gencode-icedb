"""
Analyze transcript annotations using EvidenceSupport that has been collected.

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
from collections import defaultdict
from pycbio.sys import fileOps
from gencode_icedb.tsl.supportDefs import EvidenceSupport, TrascriptionSupportLevel
from gencode_icedb.tsl.supportDefs import transIsSingleExon, geneIsTslIgnored
from gencode_icedb.tsl.tslModels import GencodeTranscriptSupport


class AnnotationSupport(defaultdict):
    """Collection of evidence supporting a transcript annotation, index by
    evidence set uuid"""
    def __init__(self, transAnnot):
        self.transAnnot = transAnnot
        super(AnnotationSupportCollector, self).__init__(list)

    def add(self, supportEval):
        self[supportEval.evidSetUuid].append(supportEval)


class SupportCollector(defaultdict):
    """Loads support for"""
    pass


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

    def _countFullSupport(self, supportEvals):
        """count support from normal and suspect evidence"""
        supporting = [ev for ev in supportEvals if ev.support < EvidenceSupport.poor]
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


class SupportClassifier(object):
    """
    Classify transcripts and/or collect statistics based on evidence from a set
    of database.
    """
    # Process by evidence database for all genes to not reopening evidence databases
    pass
