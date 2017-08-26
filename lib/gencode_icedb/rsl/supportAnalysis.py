"""
Common functions to support analysis.
"""

from gencode_icedb.rsl.rslModel import GencodeSupport, GencodeNovel


def isPseudo(suppRec):
    # FIXME: switch to using functions now in CCDS tree
    return (suppRec.transcriptType.find("pseudogene") >= 0) and (suppRec.transcriptType != "polymorphic_pseudogene")


def _intronSupportReaderFilter(isNovelDb, rec):
    """for novel, don't keep ones generated from GENCODE , for genes, don't
    use pseudo"""
    if isNovelDb:
        return rec.numExprs > 0
    else:
        return not isPseudo(rec)


def intronSupportReader(isNovelDb, chrom=None):
    """read with filtering introns support from database.  Should have already
    connected and bound ORM.  If chrom is supplied, only query that chrom for
    testing purposes."""
    ormCls = GencodeNovel if isNovelDb else GencodeSupport
    query = ormCls.select()
    if chrom is not None:
        query = query.where(ormCls.chrom == chrom)
    for rec in query:
        if _intronSupportReaderFilter(isNovelDb, rec):
            yield rec


def intronSupportAlreadyProcessed(rec, processed):
    """check in an intron has already been processed; handles duplication by gene"""
    key = (rec.chrom, rec.intronStart, rec.intronEnd, rec.strand)
    if key in processed:
        return True
    else:
        processed.add(key)
        return False
