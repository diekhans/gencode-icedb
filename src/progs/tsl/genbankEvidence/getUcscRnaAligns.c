#include "common.h"
#include "options.h"
#include "psl.h"
#include "hash.h"
#include "hdb.h"
#include "jksql.h"
#include "estOrientInfo.h"
#include "verbose.h"
#include "dystring.h"
#include "sqliteEz.h"

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "getUcscRnaAligns ucscDb type sqliteDb sqliteTable\n"
        "\n"
        "Load PSL alignments from UCSC all_mrna or all_est tables into an SQLite\n"
        "database.  EST PSLs will be reverse-complement if estOrientInfo table\n"
        "indicates.  Type is `rna' or `est'.\n"
        "\n"
        "Options:\n"
        "  -verbose=n\n"
        "  -chrom=chrom - restrict to this chrom, for testing\n"
        ;
    errAbort("%s:\n%s", msg, usageMsg);
}
static struct optionSpec optionSpecs[] = {
    {"chrom", OPTION_STRING},
    {NULL, 0}
};

static char *pslCreateSqliteTbl =
    "CREATE TABLE {table} ("
    "bin int unsigned not null,"
    "matches int unsigned not null,"
    "misMatches int unsigned not null,"
    "repMatches int unsigned not null,"
    "nCount int unsigned not null,"
    "qNumInsert int unsigned not null,"
    "qBaseInsert int unsigned not null,"
    "tNumInsert int unsigned not null,"
    "tBaseInsert int unsigned not null,"
    "strand text not null,"
    "qName text not null,"
    "qSize int unsigned not null,"
    "qStart int unsigned not null,"
    "qEnd int unsigned not null,"
    "tName text not null,"
    "tSize int unsigned not null,"
    "tStart int unsigned not null,"
    "tEnd int unsigned not null,"
    "blockCount int unsigned not null,"
    "blockSizes blob not null,"
    "qStarts text not null,"
    "tStarts text not null);";
static char *pslInsertSqlite =
    "INSERT INTO {table} VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16, ?17, ?18, ?19, ?20, ?21, ?22);";
static char *pslCreateSqliteBinIndex =
    "CREATE INDEX {table}_tName_bin on {table} (tName, bin);";
static char *pslCreateSqliteQnameIndex =
    "CREATE INDEX {table}_qname on {table} (qName);";

/* load PSLs and sort by target */
static struct psl* loadPsls(struct sqlConnection *hgConn, char *table,
                            char *restrictChrom) {
    char *sqlTemplate =
        "SELECT matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, "
        "tBaseInsert, strand, concat(qName,\".\",version), qSize, qStart, qEnd, tName, "
        "tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts "
        "FROM %s, hgFixed.gbCdnaInfo where (qName = acc) %s;";
    char where[256];
    where[0] = '\0';
    if (restrictChrom != NULL) {
        safef(where, sizeof(where), " AND (tName = \"%s\")", restrictChrom);
    }
    char sql[2048];
    safef(sql, sizeof(sql), sqlTemplate, table, where);
    struct sqlResult *sr = sqlGetResult(hgConn, sql);
    char **row;
    struct psl *psls = NULL;
    while ((row = sqlNextRow(sr)) != NULL) {
        slAddHead(&psls, pslLoad(row));
    }
    slSort(&psls, pslCmpTarget);
    return psls;
}

/* get hash key for estOrientInfo records. WARNING static return */
static char *getOrientInfoKey(char *name, char *chrom,
                              unsigned chromStart, unsigned chromEnd) {
    static char key[1024];
    safef(key, sizeof(key), "%s@%s:%d-%d", name, chrom, chromStart, chromEnd);
    return key;
}

/* load EST orient info, hashed by name, chrom, chromStart, chromEnd.
 * record are in hash local memory. */
struct hash* loadEstOrientInfos(struct sqlConnection *hgConn, char *table,
                                char *restrictChrom) {
    struct hash* orientInfoMap = hashNew(18);
    char sql[1024], where[256];
    where[0] = '\0';
    if (restrictChrom != NULL) {
        safef(where, sizeof(where), " WHERE (chrom = \"%s\")", restrictChrom);
    }
    safef(sql, sizeof(sql), "SELECT * FROM %s%s", table, where);
    struct sqlResult *sr = sqlGetResult(hgConn, sql);
    char **row;
    while ((row = sqlNextRow(sr)) != NULL) {
        struct estOrientInfo *eoi = estOrientInfoLoadLm(row+1, orientInfoMap->lm);
        char *key = getOrientInfoKey(eoi->name, eoi->chrom, eoi->chromStart, eoi->chromEnd);
        hashAdd(orientInfoMap, key, eoi);
    }
    sqlFreeResult(&sr);
    return orientInfoMap;
}

/* check is a psl looks reversed */
static boolean isPslReversed(struct hash* orientInfoMap, struct psl *psl) {
    char *key = getOrientInfoKey(psl->qName, psl->tName, psl->tStart, psl->tEnd);
    struct estOrientInfo *orientInfo = hashFindVal(orientInfoMap, key);
    return (orientInfo != NULL) && (orientInfo->intronOrientation < 0);
}

/* orient ESTs when possible */
static void orientEstPsls(struct sqlConnection *hgConn, struct psl *psls, char *restrictChrom) {
    struct hash* orientInfoMap = loadEstOrientInfos(hgConn, "estOrientInfo", restrictChrom);
    for (struct psl *psl = psls; psl != NULL; psl = psl->next) {
        if (isPslReversed(orientInfoMap, psl)) {
            pslRc(psl);
        }
    }
    hashFree(&orientInfoMap);
}

/* load and orient alignments */
static struct psl *loadAligns(char *ucscDb, char *type, char *restrictChrom) {
    struct sqlConnection *hgConn = sqlConnect(ucscDb);
    bool isEst = sameString(type, "est");
    struct psl *psls = loadPsls(hgConn, (isEst? "all_est" : "all_mrna"), restrictChrom);
    if (isEst) {
        orientEstPsls(hgConn, psls, restrictChrom);
    }
    sqlDisconnect(&hgConn);
    return psls;
}

/* autoSql pack an unsigned array into a string */
static void strPackUnsignedArray(struct dyString *strBuf,
                                 int count,
                                 unsigned *values) {
    dyStringClear(strBuf);
    for (int i = 0; i < count; i++) {
        dyStringPrintf(strBuf, "%d,", values[i]);
    }
}

/* prepare insert statement */
static sqlite3_stmt* prepPslInsert(sqlite3 *conn,
                                   char *table) {
    char *sql = sqliteEzSubTable(table, pslInsertSqlite);
    sqlite3_stmt* stmt = NULL;
    if (sqlite3_prepare_v2(conn, sql, -1, &stmt, NULL) != SQLITE_OK) {
        errAbort("error preparing open sqlite3 statement %s: \"%s\"", sqlite3_errmsg(conn), sql);
    }
    freeMem(sql);
    return stmt;
}

/* write a single PSL to the database */
static void writePslToDb(struct psl *psl,
                         sqlite3_stmt* stmt) {
    static struct dyString *blockSizesBuf = NULL, *qStartsBuf = NULL, *tStartsBuf = NULL;
    if (blockSizesBuf == NULL) {
        blockSizesBuf = dyStringNew(512);
        qStartsBuf = dyStringNew(512);
        tStartsBuf = dyStringNew(512);
    }
    strPackUnsignedArray(blockSizesBuf, psl->blockCount, psl->blockSizes);
    strPackUnsignedArray(qStartsBuf, psl->blockCount, psl->qStarts);
    strPackUnsignedArray(tStartsBuf, psl->blockCount, psl->tStarts);

    sqliteExBindInt(stmt, 1, hFindBin(psl->tStart, psl->tEnd));
    sqliteExBindInt(stmt, 2, psl->match);
    sqliteExBindInt(stmt, 3, psl->misMatch);
    sqliteExBindInt(stmt, 4, psl->repMatch);
    sqliteExBindInt(stmt, 5, psl->nCount);
    sqliteExBindInt(stmt, 6, psl->qNumInsert);
    sqliteExBindInt(stmt, 7, psl->qBaseInsert);
    sqliteExBindInt(stmt, 8, psl->tNumInsert);
    sqliteExBindInt(stmt, 9, psl->tBaseInsert);
    sqliteExBindText(stmt, 10, psl->strand);
    sqliteExBindText(stmt, 11, psl->qName);
    sqliteExBindInt(stmt, 12, psl->qSize);
    sqliteExBindInt(stmt, 13, psl->qStart);
    sqliteExBindInt(stmt, 14, psl->qEnd);
    sqliteExBindText(stmt, 15, psl->tName);
    sqliteExBindInt(stmt, 16, psl->tSize);
    sqliteExBindInt(stmt, 17, psl->tStart);
    sqliteExBindInt(stmt, 18, psl->tEnd);
    sqliteExBindInt(stmt, 19, psl->blockCount);
    sqliteExBindBlob(stmt, 20, dyStringContents(blockSizesBuf), dyStringLen(blockSizesBuf));
    sqliteExBindBlob(stmt, 21, dyStringContents(qStartsBuf), dyStringLen(qStartsBuf));
    sqliteExBindBlob(stmt, 22, dyStringContents(tStartsBuf), dyStringLen(tStartsBuf));

    if (sqlite3_step(stmt) != SQLITE_DONE) {
        errAbort("PSL insert failed %s", sqlite3_errmsg(sqlite3_db_handle(stmt)));
    }
    sqlite3_reset(stmt);
}

/* write a PSLs to the database */
static void writePslsToDb(struct psl *psls, sqlite3 *conn, char *table) {
    sqlite3_stmt* stmt = prepPslInsert(conn, table);
    for (struct psl *psl = psls; psl != NULL; psl = psl->next) {
        writePslToDb(psl, stmt);
    }
    if (sqlite3_finalize(stmt) != SQLITE_OK) {
        errAbort("error finalizing sqlite3 statement %s", sqlite3_errmsg(conn));
    }
}

/* load PSLs into sqlite database, adding bin */
static void storeSqliteDb(struct psl *psls, char *sqliteDb, char *table) {
    sqlite3 *conn = sqliteEzOpen(sqliteDb, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
    sqliteEzExec(conn, "BEGIN TRANSACTION;");
    sqliteEzCreateTable(conn, table, pslCreateSqliteTbl, TRUE);
    writePslsToDb(psls, conn, table);
    sqliteEzExec(conn, "COMMIT TRANSACTION;");
    sqliteEzExecTable(conn, table, pslCreateSqliteBinIndex);
    sqliteEzExecTable(conn, table, pslCreateSqliteQnameIndex);
    sqliteEzClose(conn);
}

/* get rna or est alignments */
static void getUcscRnaAligns(char *ucscDb, char *type, char *sqliteDb, char *sqliteTable,
                             char *restrictChrom) {
    struct psl *psls = loadAligns(ucscDb, type, restrictChrom);
    storeSqliteDb(psls, sqliteDb, sqliteTable);
    // Don't bother: pslFreeList(&psls);
}

/* entry */
int main(int argc, char** argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5)
        usage("wrong # args");
    if (!(sameString(argv[2], "rna") || sameString(argv[2], "est"))) {
        usage("expected type of `rna' or `est'");
    }
    char *restrictChrom = optionVal("chrom", NULL);
    getUcscRnaAligns(argv[1], argv[2], argv[3], argv[4], restrictChrom);
    return 0;
}
