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
    __createSql = """CREATE TABLE {table} (
            startAcc text not null,
            endAcc text not null,
            reason text not null);"""
    __insertSql = """INSERT INTO {table} (startAcc, endAcc, reason) VALUES (?, ?, ?);"""
    __indexSql = """CREATE UNIQUE INDEX {table}_startAcc on {table} (startAcc);
                    CREATE UNIQUE INDEX {table}_endAcc on {table} (endAcc);"""

    def __init__(self, conn, table, create=False):
        super(GenbankProblemCaseDbTable, self).__init__(conn, table)
        if create:
            self.create()

    def create(self):
        """create table"""
        self._create(self.__createSql)

    def index(self):
        """create index after loading"""
        self._index(self.__indexSql)

    def loads(self, rows):
        """load rows into table.  Each element of row is a list, tuple, or GenbankProblemCase"""
        self._inserts(self.__insertSql, rows)

    def get(self, acc):
        "retrieve a record by accession, or None"
        sql = "SELECT startAcc, endAcc, reason FROM {table} WHERE (startAcc >= ?) and (endAcc <= ?);"
        row = next(self.query(sql, acc, acc), None)
        if row is None:
            return None
        else:
            return GenbankProblemCase(*row)
