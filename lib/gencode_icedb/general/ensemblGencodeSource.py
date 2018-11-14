"""
Read GENCODE or other Ensembl gene build annotations from an Ensembl
database.
"""
from pycbio.hgdata import hgDb
from gencode_icedb.general.dataOps import ensureList, isChrYPar, ensemblIdSplit
from gencode_icedb.general.geneAnnot import geneAnnotGroup
from gencode_icedb.general.ensemblDbAnnotFeatures import EnsemblDbAnnotationFactory
from gencode_icedb.general.ensemblDbQuery import ensemblGeneIdRecQuery, ensemblTransQuery, ensemblGeneQuery


class EnsemblGencodeReader(object):
    """Object for accessing a GENCODE and other Ensembl annotations from the Ensembl MySql data base.
    """
    def __init__(self, ensemblDb, genomeReader=None, filterChrYPar=True,
                 transcriptTypes=None):
        # FIXME: need to remove UCSC assumption here
        self.conn = hgDb.connect(ensemblDb, useAutoSqlConv=False)
        self.filterChrYPar = filterChrYPar
        self.transcriptTypes = frozenset(transcriptTypes) if transcriptTypes is not None else None
        self.annotFactory = EnsemblDbAnnotationFactory(genomeReader)

    def close(self):
        self.conn.close()
        self.conn = None

    def _makeTransAnnot(self, ensTrans):
        "will return None if chrY PAR trans and these are being filtered"
        transAnnot = self.annotFactory.fromEnsemblDb(ensTrans)
        if self.filterChrYPar and isChrYPar(transAnnot):
            return None
        elif (self.transcriptTypes is not None) and (transAnnot.attrs.transcriptType not in self.transcriptTypes):
            return None
        else:
            return transAnnot

    def getGeneIds(self):
        return sorted(set([rec.geneId for rec in ensemblGeneIdRecQuery(self.coords)]))

    def _getByTranscriptId(self, transId):
        """Get annotations as TranscriptFeatures by transcript id.  For
        PAR regions, multiple transcript can be returned.
        """
        ensTranses = ensemblTransQuery(self.conn, transId)
        if len(ensTranses) == 0:
            raise Exception("no transcripts with id {}".format(transId))
        return [self._makeTransAnnot(ensTrans) for ensTrans in ensTranses]

    def getByTranscriptIds(self, transIds):
        """Get an annotation as a TranscriptFeatures by transcript id or ids.  For PAR regions,
        multiple transcript can be returned.
        """
        transIds = ensureList(transIds)
        transAnnots = []
        for transId in transIds:
            transAnnots.extend(self._getByTranscriptId(transId))
        return transAnnots

    def getByGeneId(self, geneId):
        "get all transcripts associated with a geneId"
        ensTranses = ensemblGeneQuery(self.conn, geneId)
        if len(ensTranses) == 0:
            raise Exception("no genes with id {}".format(geneId))
        return [self._makeTransAnnot(ensTrans) for ensTrans in ensTranses]

    def _getByGencodeIds(self, gencodeIds):
        """get annotations as TranscriptFeatures by gene or transcript id (or ids)"""
        gencodeIds = ensureList(gencodeIds)
        geneIds, transIds = ensemblIdSplit(gencodeIds)
        transAnnots = []
        for geneId in geneIds:
            transAnnots.extend(self.getByGeneId(geneId))
        for transId in transIds:
            transAnnots.extend(self._getByTranscriptId(transId))
        return transAnnots

    def getGenesByGencodeIds(self, gencodeIds):
        """get annotations as GeneAnnotation object, containing the gene's
        TranscriptFeatures.  This can be by gene or transcript id (or list of
        ids)"""
        return geneAnnotGroup(self._getByGencodeIds(gencodeIds))

    def getTranscriptsOverlapping(self, chrom, start, end, strand=None):
        """generator get overlapping annotations as TranscriptFeatures"""
        raise Exception("EnsemblDbAnnotationFactory.getTranscriptsOverlapping not yet implemented")
