"""
Query Ensembl's database for gene annotations.
"""
from collections import namedtuple
from pycbio.db import mysqlOps
import MySQLdb   # mysqlclient is required for python 3
import MySQLdb.cursors
from python.hgdata.coords import Coords

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


class EnsemblGeneTrans(namedtuple("EnsemblGeneTrans",
                                  ("geneDbId",
                                   "geneId",
                                   "geneType",
                                   "geneSource",
                                   "geneName",
                                   "transcriptDbId",
                                   "transcriptId",
                                   "transcriptType",
                                   "transcriptSource",
                                   "transcriptStart",
                                   "transcriptEnd",
                                   "transcriptStrand"))):
    """
    One row for a transcript, including gene information
    """
    __slots__ = ()
    pass


# Select that generates EnsemblGeneTrans objects
_ensemblGeneTransSelect = """SELECT
    g.gene_id AS geneDbId,
    concat(g.stable_id, ".", g.version) AS geneId,
    g.biotype AS geneType,
    g.source AS geneSource,
    gxr.display_label AS geneName,
    t.transcript_id AS transcriptDbId,
    concat(t.stable_id, ".", t.version) AS transcriptId,
    t.biotype AS transcriptType,
    t.source AS transcriptSource,
    t.seq_region_start-1 AS transcriptStart,
    t.seq_region_end AS transcriptEnd,
    IF(t.seq_region_strand < 0, "-", "+") AS transcriptStrand
FROM
    gene AS g
    LEFT JOIN xref AS gxr ON (g.display_xref_id = gxr.xref_id)
    LEFT JOIN transcript AS t ON (t.gene_id = g.gene_id)
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = t.seq_region_id)
WHERE (g.source NOT IN ("LRG database")) {extrawhere};"""


def ensemblGeneTransQuery(conn, coords=None):
    """Select EnsemblGeneTrans records, optionally limiting by location
    """
    extrawhere = ""
    args = None
    if coords is not None:
        if coords.start is None:
            extrawhere = "AND (chrom = ?)"
            args = (coords.name,)
        else:
            extrawhere = "AND (chrom = ?) AND (t.seq_region_end > ?) AND (t.seq_region_start-1 < ?)"
            args = (coords.name, coords.start, coords.end)

    yield from _query(conn, EnsemblGeneTrans,
                      _ensemblGeneTransSelect.format(extrawhere=extrawhere),
                      args)


class EnsemblTransExon(namedtuple("EnsemblTransExon",
                                  ("geneDbId",
                                   "transcriptDbId",
                                   "chrom",
                                   "transStart",
                                   "transEnd",
                                   "transStrand",
                                   "exonStart",
                                   "exonEnd",
                                   "startPhase",
                                   "endPhase"))):
    """
    One row for a transcript exon
    """
    __slots__ = ()
    pass


# Select that generates EnsemblTransExon objects
_ensemblTransExonsSelect = """SELECT
    t.gene_id AS geneDbId,
    t.transcript_id AS transcriptDbId,
    sr.name AS chrom,
    t.seq_region_start-1 AS transStart,
    t.seq_region_end AS transEnd,
    t.seq_region_strand AS transStrand,
    e.seq_region_start-1 AS exonStart,
    e.seq_region_end AS exonEnd,
    e.phase AS startPhase,
    e.end_phase AS endPhase
FROM
    transcript AS t
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = t.seq_region_id)
    LEFT JOIN exon_transcript AS et ON (et.transcript_id = t.transcript_id)
    LEFT JOIN exon AS e ON (e.exon_id = et.exon_id);
WHERE WHERE (t.source NOT IN ("LRG database")) {extrawhere};
"""


def ensemblTransExonQuery(conn, geneDbIds):
    """Select EnsemblGeneExon records for a list of genes.
    """
    for i in range(0, len(geneDbIds), _maxSetSelect):
        args = tuple(geneDbIds[i:_maxSetSelect])
        extrawhere = " AND t.gene_id IN ({})".format(",".join(len(args) * ["?"]))
        yield from _query(conn, EnsemblTransExon,
                          _ensemblTransExonsSelect.format(extrawhere=extrawhere),
                          args)


# Select that generates EnsemblTransAttr objects
_ensemblTransAttrSelect = """SELECT
    t.gene_id AS geneDbId,
    t.transcript_id AS transcriptDbId,
    concat(t.stable_id, ".", t.version) AS transcriptId,
    attrt.code AS transAttrCode,
    tattr.value AS transAttrValue
FROM
    transcript AS t
    LEFT JOIN transcript_attrib AS tattr ON (tattr.transcript_id = t.transcript_id)
    LEFT JOIN attrib_type as attrt ON  (attrt.attrib_type_id = tattr.attrib_type_id)
WHERE (t.source NOT IN ("LRG database"))
      and (attrt.code IN ("Frameshift",
                          "mRNA_start_NF",
                          "mRNA_end",
                          "cds_start_NF",
                          "cds_end_NF",
                          "upstream_ATG",
                          "readthrough_tra",
                          "gencode_basic",
                          "appris",
                          "TSL",
                          "ccds_transcript"));
"""


class EnsemblTransAttr(namedtuple("EnsemblTransAttr",
                                  ("geneDbId",
                                   "transcriptDbId",
                                   "transcriptId",
                                   "transAttrCode",
                                   "transAttrValue"))):
    """
    One row for a transcript attribute
    """
    __slots__ = ()
    pass


def ensemblTransAttrQuery(conn, geneDbIds):
    """Select EnsemblGeneAttr records for a list of genes.
    """
    for i in range(0, len(geneDbIds), _maxSetSelect):
        args = tuple(geneDbIds[i:_maxSetSelect])
        extrawhere = " AND t.gene_id IN ({})".format(",".join(len(args) * ["?"]))
        yield from _query(conn, EnsemblTransAttr,
                          _ensemblTransAttrSelect.format(extrawhere=extrawhere),
                          args)


# Select that generates Coords objects of PAR
_ensemblParCoordsSelect = """SELECT
    sr.name AS name,
    ae.seq_region_start-1 AS start,
    ae.seq_region_end AS end
FROM
    assembly_exception as ae
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = ae.seq_region_id)
WHERE ae.exc_type = "PAR";
"""


def ensemblParQuery(conn):
    """Select PAR regions as Coord objects from Ensembl Db"""
    yield from _query(conn, Coords, _ensemblParCoordsSelect)
