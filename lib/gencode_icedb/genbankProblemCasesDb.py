"""
GenBank problem case accession and reason table.
"""
from collections import namedtuple
import sqlite3
from pycbio.hgdata.hgLite import HgLiteTable
from pycbio.sys import symEnum

GenbankProblemReason = symEnum.SymEnum("GenbankProblemReason",
                                       ("nedo", "athRage", "orestes"))


class GenbankProblemCase(namedtuple("GenbankProblemCase",
                                    ("acc", "reason"))):
    """a problem case"""
    __slots__ = ()
    pass

sqlite3.register_adapter(GenbankProblemReason, str)
sqlite3.register_converter("reason", GenbankProblemReason)


class GenbankProblemCaseDbTable(HgLiteTable):
    """
    Storage for problem cases
    """
    __createSql = """CREATE TABLE {table} (
            acc text not null,
            reason text not null);"""
    __insertSql = """INSERT INTO {table} (acc, reason) VALUES (?, ?);"""
    __indexSql = """CREATE UNIQUE INDEX {table}_acc on {table} (acc);"""

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
        """load rows into table.  Each element of row is a list, tuple, or GenbankProblemCase object of name and seq"""
        self._inserts(self.__insertSql, rows)

    def get(self, acc):
        "retrieve a record by accession, or None"
        sql = "SELECT acc, reason FROM {table} WHERE acc = ?;"
        row = next(self.query(sql, acc), None)
        if row is None:
            return None
        else:
            return GenbankProblemCase(*row)
