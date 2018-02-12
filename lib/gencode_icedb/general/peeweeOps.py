"""
Functions to support common peewee operations.
"""
import re
import apsw
from playhouse.apsw_ext import APSWDatabase
from pycbio.db.sqliteOps import sqliteSetSynchronous


def peeweeConnect(sqliteDb, bindModelFunc, create=False, readonly=True, timeout=None, synchronous=None):
    """connect to sqlite3 database and bind to model,  If sqliteDb is None, create an in-memory database.
    bindModelFunc is called to bind connection to peeweemodel """
    if create:
        readonly = False
    if sqliteDb is None:
        sqliteDb = ":memory:"
        create = True
        readonly = False
    flags = apsw.SQLITE_OPEN_READONLY if readonly else apsw.SQLITE_OPEN_READWRITE
    if create:
        flags |= apsw.SQLITE_OPEN_CREATE
    kwargs = {"flags": flags}
    if timeout is not None:
        kwargs["timeout"] = timeout
    conn = APSWDatabase(sqliteDb, **kwargs)
    bindModelFunc(conn)
    if synchronous is not None:
        sqliteSetSynchronous(conn, synchronous)
    return conn


def peeweeClose(conn):
    "close database"
    # not sure why it might be in closed state even after open, maybe lazy open?
    if not conn.is_closed():
        conn.close()


def peeweeClassToTableName(cls):
    """Change tables names from came camel case to lower-case, underscore
    (e.g. BowWow to bow_wow).  To enable this in class Meta, do
    table_func = peeweeClassToTableName
    """
    return re.sub("([a-z])([A-Z])", "\\1_\\2", cls.__name__).lower()


class PeeweeModelMixins(object):
    """Adds some useful functions to a PeeWee model"""

    @classmethod
    def get_table_name(cls):
        """name of underlying database table"""
        return cls._meta.table_name

    @classmethod
    def get_fields(cls):
        "get field objects for this model"
        return list(cls._meta.fields.values())
