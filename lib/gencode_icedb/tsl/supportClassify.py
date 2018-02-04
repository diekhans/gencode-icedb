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
from collections import namedtuple
from pycbio.sys import fileOps
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.hgLite import sqliteConnect
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, ChromInsertFeature, RnaInsertFeature
from gencode_icedb.tsl.supportDefs import EvidenceSource, EvidenceEval
from gencode_icedb.tsl.supportDefs import UCSC_RNA_ALN_TBL, UCSC_EST_ALN_TBL, ENSEMBL_RNA_ALN_TBL
from gencode_icedb.general.evidFeatures import EvidencePslDbFactory

#FIXME: needed to separate from geneAnnot DB more (as might be comming from ensembl db)

# limits on size of as single indel in an exon.
exonPolymorhicSizeLimit = 12
# fraction of allowed total indel size relative exon length
exonPolymorhicFactionLimit = 0.05

def _checkMegExonIndels(evidExon):
    "check for allowed indel polymorphism"

    def getIndelSize(aln):
        if isinstance(aln, ChromInsertFeature):
            return len(aln.chrom)
        elif isinstance(aln, RnaInsertFeature):
            return len(aln.rna)
        else:
            return 0  # not an indel

    totalIndelSize = 0
    for aln in evidExon.alignFeatures:
        indelSize = getIndelSize(aln)
        if indelSize > exonPolymorhicSizeLimit:
            return EvidenceEval.large_indel_size
        totalIndelSize += indelSize
    if totalIndelSize > exonPolymorhicFactionLimit * len(evidExon.chrom):
        return EvidenceEval.large_indel_content
    else:
        return EvidenceEval.polymorphic


def _compareMegExon(annotExon, evidExon):
    # boundaries are check with introns, so just check indels
    if len(evidExon.alignFeatures) > 1:
        return _checkMegExonIndels(evidExon)
    else:
        return EvidenceEval.good


def _compareIntron(annotIntron, evidIntron):
    if annotIntron.chrom != evidIntron.chrom:
        return EvidenceEval.exon_boundry_mismatch
    elif len(evidIntron.alignFeatures) > 1:
        return EvidenceEval.internal_unaligned
    else:
        return EvidenceEval.good


def _compareFeature(annotFeat, evidFeat):
    assert type(annotFeat) == type(evidFeat)
    if isinstance(annotFeat, ExonFeature):
        return _compareMegExon(annotFeat, evidFeat)
    else:
        return _compareIntron(annotFeat, evidFeat)


def compareMegWithEvidence(annotTrans, evidTrans):
    """Compare a multi-exon annotation with a given piece of evidence"""
    if len(evidTrans.features) != len(annotTrans.features):
        return EvidenceEval.feat_mismatch
    worstCmpr = EvidenceEval.good
    for iFeat in range(len(annotTrans.features)):
        cmpr = _compareFeature(annotTrans.features[iFeat], evidTrans.features[iFeat])
        if cmpr > worstCmpr:  # > is worse
            worstCmpr = cmpr
    return worstCmpr

def _findAnnotBounds(geneTranses):
    name = geneTranses[0].chrom.name
    start = geneTranses[0].chrom.start
    end = geneTranses[0].chrom.end
    strand = geneTranses[0].rna.strand
    for geneTrans in geneTranses:
        if geneTrans.chrom.name != name:
            raise Exception("Bug: mix of chromosomes provided")
        if geneTrans.chrom.strand != '+':
            raise Exception("Bug: assumes positive chromosome strand")
        if geneTrans.rna.strand != strand:
            raise Exception("Bug: mix of RNA strand provided")
        start = min(geneTrans.chrom.start, start)
        end = max(geneTrans.chrom.end, end)
    return Coords(name, start, end, strand)


class EvidenceReader(object):
    """Container for sources of evidence."""

    evidSrcToTbl = {
        EvidenceSource.UCSC_RNA: UCSC_RNA_ALN_TBL,
        EvidenceSource.ENSEMBL_RNA: UCSC_EST_ALN_TBL,
        EvidenceSource.UCSC_EST: ENSEMBL_RNA_ALN_TBL,
    }

    def __init__(self, evidDb, genomeReader):
        self.evidDbConn = sqliteConnect(evidDb, readOnly=True)
        self.evidTbls = {evidSrc: EvidencePslDbFactory(self.evidDbConn, self.evidSrcToTbl[evidSrc], genomeReader)
                         for evidSrc in EvidenceSource}

    def get(self, evidSrc):
        "return table object to use for reading from appropriate table"
        return self.evidTbls[evidSrc]

    def close(self):
        self.evidDbConn.close()
        self.evidDbConn = None


class EvidenceCache(object):
    """Cache for a range of evidence, normally for one gene"""
    def __init__(self, evidReader, bounds):
        self.evidReader = evidReader
        self.bounds = bounds
        self.evidTbls = {evidSrc: None for evidSrc in EvidenceSource}

    def get(self, evidSrc):
        "get a type of evidence"
        if self.evidTbls[evidSrc] is None:
            # sort for testing purposes
            dbTbl = self.evidReader.get(evidSrc)
            self.evidTbls[evidSrc] = tuple(sorted(dbTbl.overlappingGen(self.bounds.name, self.bounds.start, self.bounds.end,
                                                                       rnaStrand=self.bounds.strand, minExons=2),
                                                  key=lambda a: (a.rna.name, a.chrom.name, a.chrom.start)))
        return self.evidTbls[evidSrc]


class AnnotationEvidenceEval(namedtuple("AnnotationEvidence",
                                        ("annotId", "evidSrc", "evidId", "evidEval"))):
    """Evaluation of an annotation against an evidence alignment."""
    slots = ()


class AnnotationEvidenceCollector(object):
    """Collection of evidence supporting a transcript annotation"""
    def __init__(self, annotTrans):
        self.annotTrans = annotTrans
        self.evidEvals = []

    def add(self, evidEval):
        self.evidEvals.append(evidEval)


class SupportClassifier(object):
    """TSL classifier for transcripts with evidence and annotations in sqlite3
    databases. Classification on a per-transcript basis, however it is grouped
    by gene to gives good locality of the evidence."""
    def __init__(self, evidDb, genomeReader):
        self.evidReader = EvidenceReader(evidDb, genomeReader)

    @staticmethod
    def writeTsvHeaders(tslFh, detailsFh=None):
        fileOps.prRowv(detailsFh, "transcriptId", "level")
        if detailsFh is not None:
            fileOps.prRowv(detailsFh, "transcriptId", "evidSrc", "evidId", "evidEval")

    def __writeDetails(self, detailsFh, annotTrans, evidSrc, evidTrans, evidEval):
        fileOps.prRowv(detailsFh, annotTrans.rna.name, evidSrc, evidTrans.rna.name, evidEval)

    def __compareWithEvidence(self, annotTrans, evidSrc, evidTrans, evidCollector, detailsFh):
        evidEval = compareMegWithEvidence(annotTrans, evidTrans)
        if detailsFh is not None:
            self.__writeDetails(detailsFh, annotTrans, evidSrc, evidTrans, evidEval)
        if evidEval < EvidenceEval.feat_mismatch:
            evidCollector.add(AnnotationEvidenceEval(annotTrans, evidSrc, evidTrans, evidEval))

    def __collectTransSupport(self, annotTrans, evidCache, detailsFh):
        evidCollector = AnnotationEvidenceCollector(annotTrans)
        for evidSrc in EvidenceSource:
            for evidTrans in evidCache.get(evidSrc):
                self.__compareWithEvidence(annotTrans, evidSrc, evidTrans, evidCollector, detailsFh)
        return evidCollector

    def __classifyTrans(self, annotTrans, evidCache, tslFh, detailsFh):
        evidCollector = self.__collectTransSupport(annotTrans, evidCache, detailsFh)

    def classifyGeneTranscripts(self, geneAnnotTranses, tslFh, detailsFh=None):
        """classify a list of transcripts, which must be all on the same chromosome
        and should be from the same gene."""
        evidCache = EvidenceCache(self.evidReader, _findAnnotBounds(geneAnnotTranses))
        for annotTrans in geneAnnotTranses:
            self.__classifyTrans(annotTrans, evidCache, tslFh, detailsFh)
