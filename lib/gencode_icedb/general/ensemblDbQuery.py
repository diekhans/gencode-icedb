"""
Query Ensembl's database for gene annotations.
"""
import six
import sys
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
    try:
        cur.execute(sql, args)
        for row in cur:
            yield dataCls(**row)
    finally:
        cur.close()


def _dump(fh, val, indent):
    "debugging print"
    fh.write((indent * '    ') + str(val) + '\n')


def _ensureList(strOrList):
    if isinstance(strOrList, six.string_types):
        return [strOrList]
    else:
        return strOrList


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
    xlate.seq_end AS cdsEndOffset
FROM
    gene AS g
    LEFT JOIN xref AS gxr ON (g.display_xref_id = gxr.xref_id)
    LEFT JOIN transcript AS t ON (t.gene_id = g.gene_id)
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = t.seq_region_id)
    LEFT JOIN translation AS xlate ON (xlate.transcript_id = t.transcript_id)
WHERE (g.source NOT IN ("LRG database")) {extrawhere};"""


def _idSqlInClause(idField, ids):
    return "({} IN ({}))".format(idField, ','.join(len(ids) * ['%s']))


def ensemblTransRecByTransQuery(conn, transIds):
    """Select EnsemblGeneTrans records for a transcript ids or list of transcript ids.
    """
    transIds = _ensureList(transIds)
    clause = " AND ({})".format(_idSqlInClause('concat(t.stable_id, ".", t.version)', transIds))
    yield from _query(conn, EnsemblTransRec,
                      _ensemblTransRecSelect.format(extrawhere=clause),
                      tuple(transIds))


def ensemblTransRecByGeneQuery(conn, geneIds):
    """Select EnsemblGeneTrans records for a gene id or list of gene ids.
    """
    geneIds = _ensureList(geneIds)
    clause = " AND ({})".format(_idSqlInClause('concat(g.stable_id, ".", g.version)', geneIds))
    yield from _query(conn, EnsemblTransRec,
                      _ensemblTransRecSelect.format(extrawhere=clause),
                      tuple(geneIds))


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


def _makeTrans(transRec, exonRecs, attrRecs):
    return EnsemblTranscript(transRec.transcriptDbId, transRec,
                             tuple(sorted(exonRecs, key=lambda e: e.start)),
                             tuple(sorted(attrRecs, key=lambda a: (a.code, a.value))))


def _buildTrans(conn, transRecs):
    # per-transcript data
    transDbIds = [transRec.transcriptDbId for transRec in transRecs]
    exonRecs = defaultdict(list)  # by transcript id
    for exonRec in ensemblExonRecQuery(conn, transDbIds):
        exonRecs[exonRec.transcriptDbId].append(exonRec)
    attrRecs = defaultdict(list)  # by transcript id
    for attrRec in ensemblTransAttrRecQuery(conn, transDbIds):
        attrRecs[attrRec.transcriptDbId].append(attrRec)

    # pull all together
    return [_makeTrans(transRec, exonRecs[transRec.transcriptDbId], attrRecs[transRec.transcriptDbId])
            for transRec in transRecs]


def ensemblTransQuery(conn, transIds):
    """Combined query to a list of EnsemblTranscript for a set of transcript ids"""
    return _buildTrans(conn, list(ensemblTransRecByTransQuery(conn, transIds)))


def ensemblGeneQuery(conn, geneIds):
    """Combined query to a list of EnsemblTranscript for a set of gene ids"""
    return _buildTrans(conn, list(ensemblTransRecByGeneQuery(conn, geneIds)))
