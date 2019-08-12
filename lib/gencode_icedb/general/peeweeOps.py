"""
Functions to support common peewee operations.
"""
import os
import re
import apsw
import configparser
import urllib.parse as urlparse
from collections import namedtuple
from playhouse.apsw_ext import APSWDatabase
from peewee import SqliteDatabase, MySQLDatabase, DatabaseError, CharField
from pycbio.sys.symEnum import SymEnum
from pycbio.db import sqliteOps, mysqlOps

# URL code inspired by https://github.com/kennethreitz/dj-database-url
urlparse.uses_netloc.append('mysql')
urlparse.uses_netloc.append('sqlite')


class DbType(SymEnum):
    "type of database"
    mysql = 1
    sqlite = 2


class DbUrl(namedtuple("DbUrl",
                       ("url", "scheme", "netloc", "user", "passspec", "host", "port", "database"))):
    """Parsed JDBC-like URL use to specify the database, host, and user, etc.
    """
    __slots__ = ()

    @classmethod
    def parse(cls, url):
        def emptyIsNone(s):
            return None if s in ('', None) else urlparse.unquote(s)

        p = urlparse.urlparse(url)
        if p.netloc == '':
            scheme = p.scheme if p.scheme != '' else 'file'
            return DbUrl(url=url,
                         scheme=scheme,
                         netloc=None,
                         user=None,
                         passspec=None,
                         host=None,
                         port=None,
                         database=emptyIsNone(p.path))
        else:
            # parse netloc
            m = re.match("^(((?P<user>[^:@]+)?" "(:(?P<passspec>[^@]+))?" ")@)?"
                         "(?P<host>[^:@]+)" "(:(?P<port>[0-9]+))?$", p.netloc)
            if m is None:
                raise Exception("can't parser database URL: {}".format(url))

            # only drop first `/' from database, second indicates sqlite absolute path
            return DbUrl(url=url,
                         scheme=p.scheme,
                         netloc=emptyIsNone(p.netloc),
                         user=emptyIsNone(m.group("user")),
                         passspec=emptyIsNone(m.group("passspec")),
                         host=emptyIsNone(m.group("host")),
                         port=emptyIsNone(m.group("port")),
                         database=emptyIsNone(p.path[1:]))


def _sqliteConnect(dburl, bindModelFunc, create=False, readonly=True, timeout=None, synchronous=None):
    """connect to sqlite3 database and bind to model.  If dburl has no database name,
    or the name is :memory:, an in memory database is created.
    bindModelFunc is called to bind connection to the peewee models."""
    if dburl.netloc is not None:
        raise Exception("Can't specify network location in sqlite database URL, found '{}': {}", dburl.netloc, dburl.url)
    if create:
        readonly = False
    if dburl.database is None:
        dbFile = ":memory:"
        create = True
        readonly = False
    else:
        dbFile = os.path.abspath(os.path.expanduser(dburl.database))
    flags = apsw.SQLITE_OPEN_READONLY if readonly else apsw.SQLITE_OPEN_READWRITE
    if create:
        flags |= apsw.SQLITE_OPEN_CREATE
    kwargs = {"flags": flags}
    if timeout is not None:
        kwargs["timeout"] = timeout
    conn = APSWDatabase(dbFile, **kwargs)
    bindModelFunc(conn)
    if synchronous is not None:
        sqliteOps.setSynchronous(conn, synchronous)
    try:
        conn.connect()  # connect now to get error messages
    except Exception as ex:
        raise DatabaseError("can't open database: {}".format(dbFile)) from ex
    return conn


def _loadMyCnf(dburl):
    """load my.cnf based on information in parsed URL, or an empty dict if file or
    sections is not found"""
    if dburl.passspec is None:
        return {}
    # MYCNF:SECTION or SECTION
    if dburl.passspec.find(':') >= 0:
        mycnf, section = dburl.passspec.rsplit(':', 1)
    else:
        mycnf, section = ("~/.my.cnf", dburl.passspec)
    cnf = configparser.ConfigParser().read(os.path.expanduser(mycnf))
    if cnf is None:
        return {}
    sec = cnf.get(section)
    if sec is None:
        return {}
    return sec


def _mysqlConnect(dburl, bindModelFunc):
    """connect to mysql database using the information in the dburl"
    """
    def addparam(name, kwargs, dburl, mycnf, argname=None):
        """add a parameter to kwargs if defined in dburl or mycnf.  argname allows for different name than
        in input. Arguments not found are not added to kwargs, as some can't be None."""
        v = dburl.get(name)
        if v is None:
            v = mycnf.get(name)
        if v is not None:
            kwargs[argname if argname is not None else name] = v

    mycnf = _loadMyCnf(dburl)
    kwargs = {}
    addparam("host", kwargs, dburl, mycnf)
    addparam("port", kwargs, dburl, mycnf)
    addparam("user", kwargs, dburl, mycnf)
    addparam("password", kwargs, dburl, mycnf, "passwd")
    addparam("database", kwargs, dburl, mycnf)
    mysqlOps.mySqlSetErrorOnWarn()
    conn = MySQLDatabase(dburl, **kwargs)
    conn.execute_sql("SET default_storage_engine=InnoDb;")
    bindModelFunc(conn)
    return conn


def peeweeConnect(url, bindModelFunc, create=False, readonly=True, timeout=None, synchronous=None):
    """connect to database using a bind to model. URL is described in
    configuration.md. The bindModelFunc is called to bind connection to
    peeweemodel """
    dburl = DbUrl.parse(url)
    if dburl.scheme in ("sqlite", "file"):
        return _sqliteConnect(dburl, bindModelFunc, create=create, readonly=readonly, timeout=readonly, synchronous=synchronous)
    elif dburl.scheme == "mysql":
        return _mysqlConnect(dburl, bindModelFunc)
    else:
        raise Exception("invalid scheme '{}', expect one of 'sqlite' or 'mysql': {}".format(dburl.scheme, url))


def peeweeDbType(conn):
    "return DbType of database"
    if isinstance(conn, SqliteDatabase):
        return DbType.sqlite
    elif isinstance(conn, MySQLDatabase):
        return DbType.mysql
    else:
        raise Exception("unknown peewee database type: {}".format(type(conn)))


def peeweeClose(conn):
    "close database"
    # not sure why it might be in closed state even after open, maybe lazy open?
    if not conn.is_closed():
        if peeweeDbType(conn) == DbType.sqlite:
            conn.execute_sql("PRAGMA optimize;")
        conn.close()


def peeweeClassToTableName(cls):
    """Change tables names from came camel case to lower-case, underscore
    (e.g. BowWow to bow_wow).  To enable this in class Meta, do
    table_func = peeweeClassToTableName
    """
    return re.sub("([a-z])([A-Z])", "\\1_\\2", cls.__name__).lower()


class SymEnumField(CharField):
    """A field storing an SymEnum derived field as a string in the database."""
    def __init__(self, symEnumCls, **kwargs):
        self.symEnumCls = symEnumCls
        super().__init__(**kwargs)

    def db_value(self, value):
        return str(value)

    def python_value(self, value):
        return self.symEnumCls(value)


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


def peeweeBulkLoadSetup(conn):
    conn.cache_size = -32 * 1024 * 1024
    conn.synchronous = 0
    conn.journal_mode = 0
    conn.locking_mode = "EXCLUSIVE"
    conn.count_changes = 0
    conn.temp_store = "MEMORY"
    conn.auto_vacuum = 0
