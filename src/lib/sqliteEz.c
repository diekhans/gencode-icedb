#include "sqliteEz.h"
#include "common.h"

/* open or abort */
sqlite3 *sqliteEzOpen(char *sqliteDb, int flags) {
    sqlite3 *conn = NULL;

    if (sqlite3_open_v2(sqliteDb, &conn, flags, NULL) != SQLITE_OK) {
        errAbort("can't open sqlite3 db \"%s\": %s", sqliteDb, sqlite3_errmsg(conn));
    }
    return conn;
}

/* close or abort */
void sqliteEzClose(sqlite3 *conn) {
    if (sqlite3_close_v2(conn) != SQLITE_OK) {
        errAbort("error closing sqlite3 db: %s", sqlite3_errmsg(conn));
    }
}

/* execute a sql statement that doesn't return a results */
void sqliteEzExec(sqlite3 *conn, char *sql) {
    char *errmsg = NULL;
    if (sqlite3_exec(conn, sql, NULL, NULL, &errmsg) != SQLITE_OK) {
        errAbort("sql exec failed: %s: %s", errmsg, sql);
    }
}

/* execute a sql statement that doesn't return a results, replace
* string {table} with table name in sql */
void sqliteEzExecTable(sqlite3 *conn, char *table, char *sql) {
    char *sqlTbl = sqliteEzSubTable(table, sql);
    sqliteEzExec(conn, sqlTbl);
    freeMem(sqlTbl);
}

/* create an table, optionally dropping if it exists.
 * All occurrences of the string `{table}' will be replaced
 * with table name in SQL */
void sqliteEzCreateTable(sqlite3 *conn,
                         char *table,
                         char *createSql,
                         bool drop) {
    if (drop) {
        sqliteEzExecTable(conn, table, "DROP TABLE IF EXISTS {table};");
    }
    sqliteEzExecTable(conn, table, createSql);
}
