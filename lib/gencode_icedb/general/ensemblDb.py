"""
Query Ensembl's database for gene annotations.
"""
from collections import namedtuple

# This does not use PeeWee, since the Ensembl database doesn't have foreign
# key constraints.  These would be needed to generated to PeeWee model objects
# with pwiz.


class EnsemblGeneTrans(namedtuple("EnsemblGeneTrans",
                                  ("geneDbId",
                                   "geneId",
                                   "geneType",
                                   "geneSource",
                                   "geneName",
                                   "transcriptDbId",
                                   "transcriptId",
                                   "transcriptType",
                                   "transcriptScoure",
                                   "transcriptStart",
                                   "transcriptEnd",
                                   "transcriptStrand"))):
    """
    One row for a transcript, including gene information
    """
    __slots__ = ()
    pass


# Select that generate EnsemblGeneTrans fields
_ensemblGeneTransSelect = """SELECT
    g.gene_id AS geneDbId,
    concat(g.stable_id, ".", g.version) AS geneId,
    g.biotype AS geneType,
    g.source AS geneSource,
    gxr.display_label AS geneName,
    t.transcript_id AS transcriptDbId,
    concat(t.stable_id, ".", t.version) AS transcriptId,
    t.biotype AS transcriptType,
    t.source AS transcriptScoure,
    t.seq_region_start-1 AS transcriptStart,
    t.seq_region_end AS transcriptEnd,
    IF(t.seq_region_strand < 0, "-", "+") AS transcriptStrand
FROM
    gene AS g
    LEFT JOIN xref AS gxr ON (g.display_xref_id = gxr.xref_id)
    LEFT JOIN transcript AS t ON (t.gene_id = g.gene_id)
    LEFT JOIN seq_region AS sr ON (sr.seq_region_id = t.seq_region_id)
WHERE (g.source NOT IN ("LRG database")) {extrawhere};"""
