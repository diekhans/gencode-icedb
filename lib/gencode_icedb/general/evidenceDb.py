"""
Read evidence alignments from a database.
"""
from pycbio.sys.symEnum import SymEnum
from pycbio.hgdata.hgLite import sqliteConnect, PslDbTable
from gencode_icedb.general.evidFeatures import EvidencePslFactory

# FIXME: is this too specific to TSL to be in general?


class EvidenceSource(SymEnum):
    """Source of evidence used in support"""
    __slots__ = ()
    UCSC_RNA = 1
    ENSEMBL_RNA = 2
    UCSC_EST = 3


# sqlite3 table names
UCSC_RNA_ALN_TBL = "ucsc_rna_aln"
ENSEMBL_RNA_ALN_TBL = "ensembl_rna_aln"
UCSC_EST_ALN_TBL = "ucsc_est_aln"

evidenceSourceTableMap = {
    EvidenceSource.UCSC_RNA: UCSC_RNA_ALN_TBL,
    EvidenceSource.ENSEMBL_RNA: ENSEMBL_RNA_ALN_TBL,
    EvidenceSource.UCSC_EST: UCSC_EST_ALN_TBL
}


class EvidenceReader(object):
    """Object for accessing alignment evidence data"""
    def __init__(self, evidDbFile, genomeReader=None):
        self.conn = sqliteConnect(evidDbFile)
        self.dbTables = {}  # by EvidenceSource
        for evidSrc in EvidenceSource:
            self.dbTables[evidSrc] = PslDbTable(self.conn, evidenceSourceTableMap[evidSrc])
        self.evidFactory = EvidencePslFactory(genomeReader)

    def close(self):
        self.conn.close()
        self.conn = None

    def __del__(self):
        if self.conn is not None:
            self.close()

    def genOverlapping(self, evidSrc, chrom, start, end, rnaStrand=None, minExons=0):
        """Generator of overlapping alignments as TranscriptFeatures.
        """
        dbTable = self.dbTables[evidSrc]
        strand = (None if rnaStrand is None
                  else ('+', '++') if rnaStrand == '+' else ('-', '+-'))
        extraWhere = "blockCount >= {}".format(minExons) if minExons > 0 else None
        for psl in dbTable.getTRangeOverlap(chrom, start, end, strand=strand, extraWhere=extraWhere):
            trans = self.evidFactory.fromPsl(psl)
            if len(trans.features) >= minExons + (minExons - 1):
                yield trans
