"""
GenBank problem case accession and reason table.
"""
from collections import namedtuple
import sqlite3
from pycbio.hgdata.hgLite import HgLiteTable
from pycbio.sys import symEnum


class GenbankProblemReason(symEnum.SymEnum):
    nedo = 1
    athRage = 2
    orestes = 3


class GenbankProblemCase(namedtuple("GenbankProblemCase",
                                    ("startAcc", "endAcc", "reason"))):
    """a range of problem case"""
    __slots__ = ()
    pass


sqlite3.register_adapter(GenbankProblemReason, str)
sqlite3.register_converter("reason", GenbankProblemReason)


class GenbankProblemCaseDbTable(HgLiteTable):
    """
    Storage for problem cases
    """
    _createSql = """CREATE TABLE {table} (
            startAcc text not null,
            endAcc text not null,
            reason text not null);"""
    _insertSql = """INSERT INTO {table} ({columns}) VALUES ({values});"""
    _indexSql = """CREATE UNIQUE INDEX {table}_startAcc on {table} (startAcc);
                    CREATE UNIQUE INDEX {table}_endAcc on {table} (endAcc);"""
    columnNames = ("startAcc", "endAcc", "reason")

    def __init__(self, conn, table, create=False):
        super(GenbankProblemCaseDbTable, self).__init__(conn, table)
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

    def get(self, acc):
        "retrieve a record by accession, or None"
        sql = "SELECT {columns}, reason FROM {table} WHERE (startAcc >= ?) and (endAcc <= ?);"
        row = next(self.query(sql, self.columnNames, acc, acc), None)
        if row is None:
            return None
        else:
            return GenbankProblemCase(*row)
