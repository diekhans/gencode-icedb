"""
Types related to support levels as well as some common fuctions
"""
import re
from pycbio.sys.symEnum import SymEnum, auto
from gencode_icedb.general.transFeatures import ExonFeature

# FIXME: make organization consistent between TSL and RSL stuff


class Organism(SymEnum):
    hs = 1
    mm = 2


class GenbankProblemReason(SymEnum):
    nedo = 1
    athRage = 2
    orestes = 3


class EvidenceType(SymEnum):
    """Type of evidence used in support"""
    __slots__ = ()
    RNA = auto()
    EST = auto()
    NANOPORE_DRNA = auto()
    NANOPORE_CDNA = auto()
    ISOSEQ_CDNA = auto()


class EvidenceSupport(SymEnum):
    """One or more of these attributes describe the support provided to an annotation
    by a given piece of evidence"""
    __slots__ = ()
    # FIXME: rename SupportCategory
    good = 1                    # good support
    polymorphic = 2             # minor variation

    poor = 10                   # value greater than this a detail of poor
    extends_exons = 11          # provides support with additional exons
    large_indel_size = 12       # size of an indel exceeds a threshold (excludes initial/terminal)
    large_indel_content = 14    # indel context if a given exon exceeded (excludes initial/terminal)
    internal_unaligned = 15     # intron contains unaligned
    large_indel_ends = 16       # size of an initial or terminal indel exceeds a threshold

    exon_boundry_mismatch = 50  # exon boundaries differs
    feat_count_mismatch = 51    # different number of features
    feat_mismatch = 52          # mismatch of features
    not_useful = 90             # support deemed not useful for other reasons
    no_support = 100            # values >= to this are record so we know gene was analyzed
    no_eval = 200               # generic no_eval, any of the below
    no_eval_single_exon = 201   # not evaluated due to be a single exon gene
    no_eval_olfactory_receptor = 202  # olfactory receptor gene
    no_eval_immunoglobin = 203        # immunoglobin gene
    no_eval_tcell_receptor = 204      # tcell receptor
    no_eval_hla = 205                 # HLA gene
    no_eval_pseudogene = 205          # pseudogenes that were not transcribed
    worst = 1000                      # worse than anything, for min compares

class TrascriptionSupportLevel(SymEnum):
    """Transcription support levels"""
    tslNA = -1
    tsl1 = 1
    tsl2 = 2
    tsl3 = 3
    tsl4 = 4
    tsl5 = 5


# FIXME these are from the ccds2/modules/gencode/src/lib/gencode/data/gencodeGenes.py, migrate to new module
# FIXME move these somewhere else
def transIsSingleExon(transAnnot):
    return len(transAnnot.getFeaturesOfType(ExonFeature)) <= 1


def geneTypeEvidSupport(annot):
    """ Should this gene or transcript be evaulated? Returns and
    EvidenceSupport, with good if it should be analyzed. The annot maybe gene or transcript.
    A transTypeEvidSupport must be called to check for single exon.
    """
    geneBioType = annot.attrs.geneType
    geneName = annot.attrs.geneName
    # FIXME: this was part of ccds gencodeGenes module, we need to make that independent and use it here
    # trans.isPseudo()
    if (geneBioType.find("pseudogene") >= 0) and (geneBioType != "polymorphic_pseudogene") and (geneBioType.find("transcribed") < 0) and (geneBioType.find("translated") < 0):
        return EvidenceSupport.no_eval_pseudogene
    # trans.isOlfactoryReceptor()
    if re.match("^OR[0-9].+", geneName):
        return EvidenceSupport.no_eval_olfactory_receptor
    # trans.isHLA()
    if geneName.find("HLA-") == 0:
        return EvidenceSupport.no_eval_hla
    # trans.isImmunoglobin()
    if geneBioType.startswith("IG_"):
        return EvidenceSupport.no_eval_immunoglobin
    # trans.isTCellReceptor())
    if geneBioType.startswith("TR_"):
        return EvidenceSupport.no_eval_tcell_receptor
    return EvidenceSupport.good


def transTypeEvidSupport(transAnnot):
    if transIsSingleExon(transAnnot):
        return EvidenceSupport.no_eval_single_exon
    else:
        return geneTypeEvidSupport(transAnnot)


def geneTypeIsEvaulated(annot):
    """Should this gene or transcript be evaulated? The annot can be a gene or transcript.
    This does not detect single-exon
    """
    return geneTypeEvidSupport(annot) == EvidenceSupport.good


def transTypeIsEvaulated(annot):
    """Should this transcript be evaluated?
    """
    return transTypeEvidSupport(annot) == EvidenceSupport.good
