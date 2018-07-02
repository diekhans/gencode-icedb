"""
Read GENCODE annotations from a database.
"""
import six
from collections import defaultdict
from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.hgLite import sqliteConnect, SqliteCursor, GenePredDbTable, GencodeAttrsDbTable, GencodeTagDbTable
from pycbio.hgdata.genePred import GenePred
from pycbio.hgdata.rangeFinder import Binner
from gencode_icedb.general.genePredAnnotFeatures import AnnotationGenePredFactory

# FIXME: use Coords here
# FIXME: rename module to gencodeReader.py
# FIXME: transcriptType filter should be done in sql

# tables in sqlite databases
GENCODE_ANN_TABLE = "gencode_ann"
GENCODE_ATTRS_TABLE = "gencode_attrs"
GENCODE_TRANSCRIPT_SOURCE_TABLE = "gencode_transcript_source"
GENCODE_TRANSCRIPTION_SUPPORT_LEVEL_TABLE = "gencode_transcription_support_level"
GENCODE_TAG_TABLE = "gencode_tag"


def _isChrYPar(annotTrans):
    return ("PAR" in annotTrans.attrs.tags) and (annotTrans.chrom.name == "chrY")


class UcscGencodeReader(object):
    """Object for accessing a GENCODE sqlite database with UCSC tables """
    def __init__(self, gencodeDbFile, genomeReader=None, filterChrYPar=True,
                 transcriptTypes=None):
        self.conn = None
        self.conn = sqliteConnect(gencodeDbFile)
        self.filterChrYPar = filterChrYPar
        self.transcriptTypes = frozenset(transcriptTypes) if transcriptTypes is not None else None
        self.genePredDbTable = GenePredDbTable(self.conn, GENCODE_ANN_TABLE)
        self.attrDbTable = GencodeAttrsDbTable(self.conn, GENCODE_ATTRS_TABLE)
        self.tagDbTable = GencodeTagDbTable(self.conn, GENCODE_TAG_TABLE)
        self.annotFactory = AnnotationGenePredFactory(genomeReader)

    def close(self):
        self.conn.close()
        self.conn = None

    def _getAttrs(self, gp):
        attrs = ObjDict()
        attrs.update(self.attrDbTable.getrByTranscriptId(gp.name)._asdict())
        attrs.tags = frozenset([t.tag for t in self.tagDbTable.getByTranscriptId(gp.name)])
        return attrs

    def _makeTransAnnot(self, gp):
        "will return None if chrY PAR trans and these are being filtered"
        transAnnot = self.annotFactory.fromGenePred(gp, self._getAttrs(gp))
        if self.filterChrYPar and _isChrYPar(transAnnot):
            return None
        elif (self.transcriptTypes is not None) and (transAnnot.attrs.transcriptType not in self.transcriptTypes):
            return None
        else:
            return transAnnot

    def getGeneIds(self):
        return self.attrDbTable.getGeneIds()

    def getByGeneId(self, geneId):
        transIds = self.attrDbTable.getGeneTranscriptIds(geneId)
        return self.getByTranscriptIds(transIds)

    def getByGeneIds(self, geneIds):
        "get transcripts for genes as a list of list of transcripts, sorted for reputability"
        return [self.getByGeneId(geneId) for geneId in sorted(geneIds)]

    def getStartingInBounds(self, chrom, start, end):
        """Get the annotations genes that *start* in the range. To get whole
        chrom, start and end maybe None.
        """
        columns = ",".join(GenePredDbTable.columnNames)
        if start is None:
            rangeWhere = "(chrom = '{}')".format(chrom)
        else:
            # note: end is at txStart+1
            rangeWhere = Binner.getOverlappingSqlExpr("bin", "chrom", "txStart", "txStart+1", chrom, start, end)
        sql = "SELECT {columns} FROM {table} WHERE {rangeWhere}".format(columns=columns, table=GENCODE_ANN_TABLE, rangeWhere=rangeWhere)
        transAnnots = []
        with SqliteCursor(self.conn, rowFactory=lambda cur, row: GenePred(row)) as cur:
            cur.execute(sql)
            for gp in cur:
                transAnnot = self._makeTransAnnot(gp)
                if transAnnot is not None:
                    transAnnots.append(transAnnot)
        return transAnnots

    def getByTranscriptIds(self, transIds):
        """get annotations as TranscriptFeatures by transcript id (or ids)"""
        if isinstance(transIds, six.string_types):
            transIds = [transIds]
        transAnnots = []
        for transId in transIds:
            for gp in self.genePredDbTable.getByName(transId):
                transAnnot = self._makeTransAnnot(gp)
                if transAnnot is not None:
                    transAnnots.append(transAnnot)
        return transAnnots

    def _getByGencodeId(self, gencodeId):
        """get annotations as TranscriptFeatures by a gene or transcript id"""
        attrs = self.attrDbTable.getByGeneId(gencodeId)
        if len(attrs) > 0:
            return self.getByGeneId(gencodeId)
        else:
            attrs = self.attrDbTable.getByTranscriptId(gencodeId)
            if attrs is None:
                raise Exception("Not a valid GENCODE gene or transcript id: {}".format(gencodeId))
            return self.getByTranscriptIds(gencodeId)

    def getByGencodeIds(self, gencodeIds):
        """get annotations as TranscriptFeatures by gene or transcript id (or ids)"""
        if isinstance(gencodeIds, six.string_types):
            gencodeIds = [gencodeIds]
        transAnnots = []
        for gencodeId in gencodeIds:
            transAnnots.extend(self._getByGencodeId(gencodeId))
        return transAnnots

    def getByGencodeIdsGrouped(self, gencodeIds):
        """get annotations as TranscriptFeatures by gene or transcript id (or ids),
        group into list or lists by geneId"""
        byGeneId = defaultdict(list)
        for transAnnot in self.getByGencodeIds(gencodeIds):
            byGeneId[transAnnot.attrs.geneId].append(transAnnot)
        return [[ta for ta in sorted(byGeneId[gi], key=lambda t: t.rna.name)]
                for gi in sorted(byGeneId.keys())]

    def getOverlapping(self, chrom, start, end, strand=None):
        """generator get overlapping annotations as TranscriptFeatures"""
        # FIXME: is this needed?
        transAnnots = []
        for gp in self.genePredDbTable.getRangeOverlap(chrom, start, end, strand=strand):
            transAnnot = self._makeTransAnnot(gp)
            if transAnnot is not None:
                transAnnots.append(transAnnot)
        return transAnnots

    def getAll(self):
        """Get all annotations as TranscriptFeatures"""
        transAnnots = []
        for gp in self.genePredDbTable.getAll():
            transAnnot = self._makeTransAnnot(gp)
            if transAnnot is not None:
                transAnnots.append(transAnnot)
        return transAnnots


def findAnnotationBounds(geneTranses):
    """find the bounds of a list of gene transcripts on the same chromosome"""
    annotId = geneTranses[0].rna.name
    name = geneTranses[0].chrom.name
    start = geneTranses[0].chrom.start
    end = geneTranses[0].chrom.end
    strand = geneTranses[0].rna.strand
    for geneTrans in geneTranses:
        if geneTrans.chrom.name != name:
            raise Exception("Bug: mix of chromosomes provided: {} and {}".format(geneTrans.rna.name, annotId))
        if geneTrans.chrom.strand != '+':
            raise Exception("Bug: assumes positive chromosome strand: {} and {}".format(geneTrans.rna.name, annotId))
        if geneTrans.rna.strand != strand:
            raise Exception("Bug: mix of RNA strand provided: {} and {}".format(geneTrans.rna.name, annotId))
        start = min(geneTrans.chrom.start, start)
        end = max(geneTrans.chrom.end, end)
    return Coords(name, start, end, strand)
