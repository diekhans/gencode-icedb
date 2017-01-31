/* functions to warp sqlite3 calls */
#ifndef sqliteEz_h
#define sqliteEz_h
#include "common.h"
#include "sqlite3.h"

/* open or abort */
sqlite3 *sqliteEzOpen(char *sqliteDb, int flags);

/* close or abort */
void sqliteEzClose(sqlite3 *conn);

/* substitute `{table}' with table name in sql template. */
INLINE char *sqliteEzSubTable(char *table, char *sqlTemplate) {
    return replaceChars(sqlTemplate, "{table}", table);
}

/* execute a sql statement that doesn't return a results */
void sqliteEzExec(sqlite3 *conn, char *sql);

/* execute a sql statement that doesn't return a results, replace
* string {table} with table name in sql */
void sqliteEzExecTable(sqlite3 *conn, char *table, char *sql);

/* create an table, optionally dropping if it exists.
 * All occurrences of the string `{table}' will be replaced
 * with table name in SQL */
void sqliteEzCreateTable(sqlite3 *conn,
                         char *table,
                         char *createSql,
                         bool drop);

/* abort on a bind error */
INLINE void sqliteExBindErrAbort(sqlite3_stmt* stmt, int paramIdx) {
    errAbort("sqlite3 error binding parameter %d: %s", paramIdx, sqlite3_errmsg(sqlite3_db_handle(stmt)));
}

/* bind_blob that aborts */
INLINE void sqliteExBindBlob(sqlite3_stmt* stmt, int paramIdx,
                             const void* value, int size) {
    if (sqlite3_bind_blob(stmt, paramIdx, value, size, SQLITE_STATIC) != SQLITE_OK) {
        sqliteExBindErrAbort(stmt, paramIdx);
    }
}

/* bind_double that aborts */
INLINE void sqliteExBindDouble(sqlite3_stmt* stmt, int paramIdx, double value) {
    if (sqlite3_bind_double(stmt, paramIdx, value) != SQLITE_OK) {
        sqliteExBindErrAbort(stmt, paramIdx);
    }
}

/* bind_ intthat aborts */
INLINE void sqliteExBindInt(sqlite3_stmt* stmt, int paramIdx, int value) {
    if (sqlite3_bind_int(stmt, paramIdx, value) != SQLITE_OK) {
        sqliteExBindErrAbort(stmt, paramIdx);
    }
}

/* bind_text that aborts */
INLINE void sqliteExBindText(sqlite3_stmt* stmt, int paramIdx,
                             const char* value) {
    if (sqlite3_bind_text(stmt, paramIdx, value, -1, SQLITE_STATIC) != SQLITE_OK) {
        sqliteExBindErrAbort(stmt, paramIdx);
    }
}
#endif
