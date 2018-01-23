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
from gencode_icedb.tsl.supportDefs import EvidenceComparison


def compareMegExon(annotExon, evidExon):
    pass


def compareFeature(annotExon, evidExon):
    pass


def compareMegWithEvidence(annotTrans, evidTrans):
    """Compare a multi-exon annotation with a given piece of evidence"""
    if len(evidTrans.features) != len(annotTrans.features):
        return EvidenceComparison.featMismatch

    bestCmp = EvidenceComparison.good
    for iFeat in range(len(annotTrans.features)):
        pass
