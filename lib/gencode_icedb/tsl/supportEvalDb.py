"""
Results of support evaluation of transcripts against a single evidence sources.
These can be stored in a file or database.  Output is to a TSV which are then
loaded into a database.
"""
from collections import namedtuple
from uuid import UUID
from gencode_icedb.tsl.supportDefs import EvidenceSupport
from gencode_icedb.general.dataOps import outCnv, inCnv


class SupportEvidEvalResult(namedtuple("SupportEvidEvalEval",
                                       ("transcriptId", "evidSetUuid", "evidId", "support",
                                        "offset5", "offset3",
                                        "extend5Exons", "extend3Exons"))):
    """Evaluation of an annotation against an item of evidence.
    This is the output of collection jobs that will later be combined
    and stored in a database."""
    slots = ()

    def __new__(cls, transcriptId, evidSetUuid, evidId, support, offset5=0, offset3=0, extend5Exons=0, extend3Exons=0):
        return super(SupportEvidEvalResult, cls).__new__(cls, transcriptId, inCnv(evidSetUuid, UUID), evidId,
                                                         EvidenceSupport(support),
                                                         int(offset5), int(offset3),
                                                         int(extend5Exons), int(extend5Exons))

    @classmethod
    def tsvHeader(cls):
        return cls._fields

    def toRow(self):
        """convert to a TVS row"""
        return [outCnv(getattr(self, f)) for f in self._fields]


class SupportEvalResult(namedtuple("SupportEvalResult",
                                   ("transcriptId", "evidSetUuid", "support", "evidCount",
                                    "offset5", "offset3",
                                    "extend5Exons", "extend3Exons"))):
    """Results of evaluating a transcript against an evidence data set.  This
    records the best support that was found in the dataset.  Non-supporting
    overlapping evidence is not reported.  This is used to store results
    from a batch before loading into a database or other analysis.  The details
    of the specific support

    :param transcriptId: GENCODE transcript identifier
    :param evidSetUuid:' GUID of evidence set providing support
    :param support: EvidenceSupport best support by this evidence set.
    :param count: Number of evidence alignments at this support level.
    :param offset5:  Offset of the start of the first evidence exon to the start of
           annotation.  Positive indicates the evidence is longer than annotation,
           negative is shorter. This is in the direction of transcription.  This is
           the longest support of any supporting evidence.
    :param offset3: Offset of the end of the first evidence exon to the end of
           annotation.  Interpretation is the same as offset3.
    :para extend5Exons: Maximum extending 5' number of exons.
    :para extend3Exons: Maximum extending 3' number of exons.
    """
    __slots__ = ()

    def __new__(cls, transcriptId, evidSetUuid, support, evidCount, offset5, offset3, extend5Exons=0, extend3Exons=0):
        return super(SupportEvalResult, cls).__new__(cls, transcriptId,
                                                     inCnv(evidSetUuid, UUID),
                                                     EvidenceSupport(support),
                                                     int(evidCount),
                                                     int(offset5), int(offset3),
                                                     int(extend5Exons), int(extend5Exons))

    @classmethod
    def tsvHeader(cls):
        return cls._fields

    def toRow(self):
        """convert to a TVS row"""
        return [str(getattr(self, f)) for f in self._fields]
