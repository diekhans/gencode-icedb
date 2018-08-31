"""
Query Ensembl's database for gene annotations.
"""
import sys
import re
from collections import namedtuple, defaultdict
from pycbio.db import mysqlOps
import MySQLdb   # mysqlclient is required for python 3
import MySQLdb.cursors
from pycbio.hgdata.coords import Coords

# max number of IN queries in one query
_maxSetSelect = 4096

# This does not use PeeWee, since the Ensembl database doesn't have foreign
# key constraints.  These would be needed to generated to PeeWee model objects
# with pwiz.

mysqlOps.mySqlSetErrorOnWarn()


def _query(conn, dataCls, sql, args=None):
    cur = conn.cursor(cursorclass=MySQLdb.cursors.DictCursor)
    #print("sql", sql)
    #print("arg", len(args), args)
    try:
        cur.execute(sql, args)
        for row in cur:
            yield dataCls(**row)
    finally:
        cur.close()


def _dump(fh, val, indent):
    "debugging print"
    fh.write((indent * '    ') + str(val) + '\n')


####################################################################################################
class EnsemblGeneIdRec(namedtuple("EnsemblGeneRec",
                                  ("geneDbId",
                                   "geneId"))):
    """Ensembl gene database id and gene id"""
    __slots__ = ()
    pass


# Select that generates EnsemblGeneIdRec objects
_ensemblGeneIdRecSelect = """
SELECT
    g.gene_id AS geneDbId,
    concat(g.stable_id, ".", g.version) AS geneId
FROM
    gene AS g
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = t.seq_region_id)
WHERE (g.source NOT IN ("LRG database")) {extrawhere};"""


def ensemblGeneIdRecQuery(conn, coords=None):
    """Select EnsemblGeneIdRec records, optionally limiting by location"""
    extrawhere = ""
    args = None
    if coords is not None:
        if coords.start is None:
            extrawhere = "AND (chrom = %s)"
            args = (coords.name,)
        else:
            extrawhere = "AND (chrom = %s) AND (t.seq_region_end > %s) AND (t.seq_region_start-1 < %s)"
            args = (coords.name, coords.start, coords.end)

    yield from _query(conn, EnsemblGeneIdRec,
                      _ensemblGeneIdRecSelect.format(extrawhere=extrawhere),
                      args)


####################################################################################################
class EnsemblTransRec(namedtuple("EnsemblGeneTrans",
                                 ("geneDbId", "geneId",
                                  "geneType", "geneSource", "geneName",
                                  "transcriptDbId", "transcriptId",
                                  "transcriptType", "transcriptSource",
                                  "chrom", "start", "end", "strand", "chromSize",
                                  "cdsStartExonDbId", "cdsStartOffset",
                                  "cdsEndExonDbId", "cdsEndOffset"))):
    """
    One row for a transcript, including gene information
    """
    __slots__ = ()
    pass


# Select that generates EnsemblTransRec objects
_ensemblTransRecSelect = """
SELECT
    g.gene_id AS geneDbId,
    concat(g.stable_id, ".", g.version) AS geneId,
    g.biotype AS geneType,
    g.source AS geneSource,
    gxr.display_label AS geneName,
    t.transcript_id AS transcriptDbId,
    concat(t.stable_id, ".", t.version) AS transcriptId,
    t.biotype AS transcriptType,
    t.source AS transcriptSource,
    sr.name AS chrom,
    t.seq_region_start-1 AS start,
    t.seq_region_end AS end,
    IF(t.seq_region_strand < 0, "-", "+") AS strand,
    sr.length AS chromSize,
    xlate.start_exon_id AS cdsStartExonDbId,
    xlate.seq_start-1 AS cdsStartOffset,
    xlate.end_exon_id AS cdsEndExonDbId,
    xlate.seq_end-1 AS cdsEndOffset
FROM
    gene AS g
    LEFT JOIN xref AS gxr ON (g.display_xref_id = gxr.xref_id)
    LEFT JOIN transcript AS t ON (t.gene_id = g.gene_id)
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = t.seq_region_id)
    LEFT JOIN translation AS xlate ON (xlate.transcript_id = t.transcript_id)
WHERE (g.source NOT IN ("LRG database")) {extrawhere};"""


def _splitEnsIds(ensIds):
    """split into gene and transcript ids"""
    # split by matching pattern
    if isinstance(ensIds, str):
        ensIds = [ensIds]
    geneIds = []
    transIds = []
    for ensId in ensIds:
        if re.match("E[A-Z]+G[0-9]+\\.[0-9]+$", ensId):
            geneIds.append(ensId)
        elif re.match("E[A-Z]+T[0-9]+\\.[0-9]+$", ensId):
            transIds.append(ensId)
        else:
            raise Exception("not a valid Ensembl id: {}".format(ensId))
    return geneIds, transIds


def _idSqlInClause(idField, ids):
    return "({} IN ({}))".format(idField, ','.join(len(ids) * ['%s']))


def ensemblTransRecQuery(conn, ensIds):
    """Select EnsemblGeneTrans records for a list of gene or transcript id.
    """
    geneIds, transIds = _splitEnsIds(ensIds)
    inClauses = []
    args = []
    if len(geneIds) > 0:
        inClauses.append("({})".format(_idSqlInClause('concat(g.stable_id, ".", g.version)', geneIds)))
        args.extend(geneIds)
    if len(transIds) > 0:
        inClauses.append("({})".format(_idSqlInClause('concat(t.stable_id, ".", t.version)', transIds)))
        args.extend(transIds)

    extrawhere = "AND " + " OR ".join(inClauses)
    yield from _query(conn, EnsemblTransRec,
                      _ensemblTransRecSelect.format(extrawhere=extrawhere),
                      tuple(args))


####################################################################################################
class EnsemblExonRec(namedtuple("EnsemblExonRec",
                                ("transcriptDbId", "exonDbId",
                                 "start", "end",
                                 "startPhase", "endPhase"))):
    """
    One row for a transcript exon
    """
    __slots__ = ()
    pass


# Select that generates EnsemblExonRec objects
_ensemblExonRecSelect = """
SELECT
    et.transcript_id AS transcriptDbId,
    e.exon_id AS exonDbId,
    e.seq_region_start-1 AS start,
    e.seq_region_end AS end,
    e.phase AS startPhase,
    e.end_phase AS endPhase
FROM
    exon_transcript AS et
    LEFT JOIN exon AS e ON (e.exon_id = et.exon_id)
WHERE {extrawhere};
"""


def ensemblExonRecQuery(conn, transcriptDbIds):
    """Select EnsemblExonRec records for a list of transcripts.
    """
    for i in range(0, len(transcriptDbIds), _maxSetSelect):
        args = tuple(transcriptDbIds[i:_maxSetSelect])
        extrawhere = "(et.transcript_id IN ({}))".format(",".join(len(args) * ["%s"]))
        yield from _query(conn, EnsemblExonRec,
                          _ensemblExonRecSelect.format(extrawhere=extrawhere),
                          args)


####################################################################################################
class EnsemblTransAttrRec(namedtuple("EnsemblTransAttrRec",
                                     ("transcriptDbId", "code", "value"))):
    """
    One row for a transcript attribute
    """
    __slots__ = ()
    pass


# Select that generates EnsemblTransAttrRec objects
_ensemblTransAttrRecSelect = """
SELECT
    ta.transcript_id AS transcriptDbId,
    atyp.code AS code,
    ta.value AS value
FROM
    transcript_attrib AS ta
    LEFT JOIN attrib_type as atyp ON  (atyp.attrib_type_id = ta.attrib_type_id)
WHERE (atyp.code IN ("Frameshift", "mRNA_start_NF", "mRNA_end",
                     "cds_start_NF", "cds_end_NF", "upstream_ATG",
                     "readthrough_tra", "gencode_basic", "appris",
                     "TSL", "ccds_transcript")) {extrawhere}
"""


def ensemblTransAttrRecQuery(conn, transDbIds):
    """Select EnsemblTransAttrRec records for a list of transcript.
    """
    for i in range(0, len(transDbIds), _maxSetSelect):
        args = tuple(transDbIds[i:_maxSetSelect])
        extrawhere = " AND ta.transcript_id IN ({})".format(",".join(len(args) * ["%s"]))
        yield from _query(conn, EnsemblTransAttrRec,
                          _ensemblTransAttrRecSelect.format(extrawhere=extrawhere),
                          args)


####################################################################################################
# Select that generates Coords objects of PAR
_ensemblParCoordsRecSelect = """
SELECT
    sr.name AS name,
    ae.seq_region_start-1 AS start,
    ae.seq_region_end AS end,
    sr.length AS size,
FROM
    assembly_exception AS ae
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = ae.seq_region_id)
WHERE ae.exc_type = "PAR";
"""


def ensemblParCoordsRecQuery(conn):
    """Select PAR regions as Coord objects from Ensembl Db"""
    yield from _query(conn, Coords, _ensemblParCoordsRecSelect)


####################################################################################################
class EnsemblGene(namedtuple("EnsemblGene",
                             ("geneDbId",
                              "geneId",
                              "transcripts"))):
    """Result records for a gene and it's transcripts.  If query is
    by-transcript, it may not have all the genes.  Multiple results are
    returned when the same gene is in different sequences (e.g. PAR).  Use
    geneDbId to distinguish."""
    __slots__ = ()

    def dump(self, fh=sys.stderr, indent=0):
        _dump(fh, self, indent)
        for t in self.transcripts:
            _dump(fh, t, indent + 1)


class EnsemblTranscript(namedtuple("EnsemblTranscript",
                                   ("transcriptDbId",
                                    "transcript", "exons", "attrs"))):
    """Result records for a transcripts."""
    __slots__ = ()

    def dump(self, fh=sys.stderr, indent=0):
        _dump(fh, self.transcript, indent)
        for e in self.exons:
            _dump(fh, e, indent + 1)
        for a in self.attrs:
            _dump(fh, a, indent + 1)


def ensemblGeneQuery(conn, ensIds):
    """Combined query to create EnsemblGene objects for a set of geneIds """
    transDbIds = set()
    # genes and transcripts
    transRecs = defaultdict(dict)  # by geneDbId then transcriptDbId
    for transRec in ensemblTransRecQuery(conn, ensIds):
        transRecs[transRec.geneDbId][transRec.transcriptDbId] = transRec
        transDbIds.add(transRec.transcriptDbId)
    transDbIds = tuple(transDbIds)

    # per-transcript data
    exonRecs = defaultdict(list)  # by transcript id
    for exonRec in ensemblExonRecQuery(conn, transDbIds):
        exonRecs[exonRec.transcriptDbId].append(exonRec)
    attrRecs = defaultdict(list)  # by transcript id
    for attrRec in ensemblTransAttrRecQuery(conn, transDbIds):
        attrRecs[attrRec.transcriptDbId].append(attrRec)

    # pull all together
    genes = []
    for geneDbId in sorted(transRecs.keys()):
        geneTrans = []
        for transDbId in sorted(transRecs[geneDbId].keys()):
            exons = tuple(sorted(exonRecs[transDbId], key=lambda e: e.start))
            attrs = tuple(sorted(attrRecs[transDbId], key=lambda a: (a.code, a.value)))
            geneTrans.append(EnsemblTranscript(transDbId,
                                               transRecs[geneDbId][transDbId],
                                               exons, attrs))
        genes.append(EnsemblGene(geneDbId, geneTrans[0].transcript.geneId, tuple(geneTrans)))

    return tuple(genes)
