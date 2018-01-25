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
from gencode_icedb.general.transFeatures import ExonFeature, IntronFeature, ChromInsertFeature, RnaInsertFeature
from gencode_icedb.tsl.supportDefs import EvidenceComparison


# limits on size of as single indel in an exon.
exonPolymorhicSizeLimit = 12
# fraction of allowed total indel size relative exon length
exonPolymorhicFactionLimit = 0.05


def checkMegExonIndels(evidExon):
    "check for allowed indel polymorphism"

    def getIndelSize(aln):
        if isinstance(aln ChromInsertFeature):
            return len(aln.chrom)
        elif isinstance(aln RnaInsertFeature):
            return len(aln.rna)
        else:
            return 0  # not an indel

    totalIndelSize = 0
    for aln in evidExon.alignFeatures:
        indelSize = getIndelSize(aln)
        if indelSize > exonPolymorhicSizeLimit:
            return EvidenceComparison.large_indel_size
        totalIndelSize += indelSize
    if totalIndelSize > exonPolymorhicFactionLimit * len(evidExon.chrom):
        return EvidenceComparison.large_indel_content
    else:
        return EvidenceComparison.polymorphic

def compareMegExon(annotExon, evidExon):
    # boundaries are check with introns, so just check indels
    if len(evidExon.alignFeatures) > 1:
        return checkMegExonIndels(evidExon)
    else:
        return EvidenceComparison.good


def compareIntron(annotIntron, evidIntron):
    if annotIntron.chrom != evidIntron.chrom:
        return EvidenceComparison.exon_boundry_mismatch
    elif len(evidIntron.alignFeatures) > 1:
        return EvidenceComparison.internal_unaligned
    else:
        return EvidenceComparison.good


def compareFeature(annotFeat, evidFeat):
    assert type(annotFeat) == type(evidFeat)
    if isinstance(annotFeat, ExonFeature):
        return compareMegExon(annotFeat, evidFeat)
    else:
        return compareIntron(annotFeat, evidFeat)


def compareMegWithEvidence(annotTrans, evidTrans):
    """Compare a multi-exon annotation with a given piece of evidence"""
    if len(evidTrans.features) != len(annotTrans.features):
        return EvidenceComparison.featMismatch
    worstCmpr = EvidenceComparison.good
    for iFeat in range(len(annotTrans.features)):
        cmpr = compareFeature(annotTrans.features[iFeat], evidTrans[iFeat])
        if cmpr > worstCmpr:  # > is worse
            worstCmpr = cmpr
    return worstCmpr
