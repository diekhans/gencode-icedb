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


# FIXME: this is tmp
fakeGenbankUuid = {
    EvidenceType.RNA: "02c995e3-372c-4cde-b216-5d3376c51988",
    EvidenceType.EST: "fcfc5121-7e07-42f6-93a2-331a513eeb2c",
}


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
    no_support = 100


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


def isGeneIgnored(annot):
    "can be gene or transcript annotation"
    bioType = annot.attrs.geneType
    geneName = annot.attrs.geneName
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


def isTransIgnored(transAnnot):
    return transIsSingleExon(transAnnot) or isGeneIgnored(transAnnot)
