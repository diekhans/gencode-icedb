"""
GenBank problem case accession and reason table.
"""
from collections import namedtuple
from pycbio.hgdata.hgSqlite import HgSqliteTable
from gencode_icedb.tsl.supportDefs import EvidenceType, GenbankProblemReason, Organism

# IMPORTANT: this does not use peewee because the rest of the evidence database
# does not use peewee.

# FIXME: merge with evidenceDb or with genbankProblemCase.py

GENBANK_PROBLEM_CASE_TBL = "genbank_problem_case"


class GenbankProblemCase(namedtuple("GenbankProblemCase",
                                    ("organism", "etype", "startAcc", "endAcc", "reason"))):
    """a range of problem case"""
    __slots__ = ()

    def __new__(cls, organism, etype, startAcc, endAcc, reason):
        """converts types when loading from strings"""
        return super(GenbankProblemCase, cls).__new__(cls, Organism(organism), EvidenceType(etype),
                                                      startAcc, endAcc, GenbankProblemReason(reason))


class GenbankProblemCaseSqliteTable(HgSqliteTable):
    """
    Storage for problem cases
    """
    _createSql = """CREATE TABLE {table} (
            organism text not null,
            etype text not null,
            startAcc text not null,
            endAcc text not null,
            reason text not null);"""
    _insertSql = """INSERT INTO {table} ({columns}) VALUES ({values});"""
    _indexSql = """CREATE UNIQUE INDEX {table}_startAcc on {table} (startAcc);
                    CREATE UNIQUE INDEX {table}_endAcc on {table} (endAcc);"""
    columnNames = ("organism", "etype", "startAcc", "endAcc", "reason")

    def __init__(self, conn, table=GENBANK_PROBLEM_CASE_TBL, create=False):
        super(GenbankProblemCaseSqliteTable, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self._createSql)

    def index(self):
        """create index after loading"""
        self._index(self._indexSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or GenbankProblemCase"""
        self._inserts(self._insertSql, self.columnNames, rows)

    def genAll(self):
        "generator to retrieve all records"
        sql = "SELECT {columns} FROM {table}"
        return self.queryRows(sql, self.columnNames,
                              lambda cur, row: GenbankProblemCase(*row))
