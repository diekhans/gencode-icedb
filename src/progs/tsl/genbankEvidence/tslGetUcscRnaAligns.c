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
    static char* usageMsg = "tslGetUcscRnaAligns ucscDb type sqliteDb\n"
        "\n"
        "Load PSL alignments from UCSC all_mrna or all_est tables into an SQLite\n"
        "database.  EST PSLs will be reverse-complement if estOrientInfo table\n"
        "indicates.  Type is `rna' or `est'.\n"
        "\n"
        "Options:\n"
        "  -verbose=n\n"
        "  -table=tbl - load this table, defaults to ucsc_rna_aln or ucsc_est_aln.\n"
        "  -chromSpec=spec - restrict to this chrom or chrom range, for testing.\n"
        "   Maybe repeated. Duplicates caused alignments being in multiple ranges\n"
        "   are discarded.\n"
        ;
    errAbort("%s:\n%s", msg, usageMsg);
}
static struct optionSpec optionSpecs[] = {
    {"chromSpec", OPTION_STRING|OPTION_MULTI},
    {"table", OPTION_STRING},
    {NULL, 0}
};

static char *UCSC_RNA_ALN_TBL = "ucsc_rna_aln";
static char *UCSC_EST_ALN_TBL = "ucsc_est_aln";


static char *pslCreateSqliteTbl =
    "CREATE TABLE {table} ("
    "bin INT UNSIGNED NOT NULL,"
    "matches INT UNSIGNED NOT NULL,"
    "misMatches INT UNSIGNED NOT NULL,"
    "repMatches INT UNSIGNED NOT NULL,"
    "nCount INT UNSIGNED NOT NULL,"
    "qNumInsert INT UNSIGNED NOT NULL,"
    "qBaseInsert INT UNSIGNED NOT NULL,"
    "tNumInsert INT UNSIGNED NOT NULL,"
    "tBaseInsert INT UNSIGNED NOT NULL,"
    "strand TEXT NOT NULL,"
    "qName TEXT NOT NULL,"
    "qSize INT UNSIGNED NOT NULL,"
    "qStart INT UNSIGNED NOT NULL,"
    "qEnd INT UNSIGNED NOT NULL,"
    "tName TEXT NOT NULL,"
    "tSize INT UNSIGNED NOT NULL,"
    "tStart INT UNSIGNED NOT NULL,"
    "tEnd INT UNSIGNED NOT NULL,"
    "blockCount INT UNSIGNED NOT NULL,"
    "blockSizes TEXT NOT NULL,"
    "qStarts TEXT NOT NULL,"
    "tStarts TEXT NOT NULL);";
static char *pslInsertSqlite =
    "INSERT INTO {table} VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16, ?17, ?18, ?19, ?20, ?21, ?22);";
static char *pslCreateSqliteBinIndex =
    "CREATE INDEX {table}_tName_bin ON {table} (tName, bin);";
static char *pslCreateSqliteQnameIndex =
    "CREATE INDEX {table}_qname ON {table} (qName);";

/* chromosome spec */
struct ChromSpec {
    struct ChromSpec *next;
    char* chrom;    // NULL if not restriction
    int start;      // zero length for whole chromosome
    int end;
};

/* a verbose message about a psl */
static void pslVerb(int level, char* msg, struct psl* psl) {
    verbose(level, "%s: %s:%d-%d <=> %s:%d-%d (%s)\n", msg,
            psl->qName, psl->qStart, psl->qEnd,
            psl->tName, psl->tStart, psl->tEnd,
            psl->strand);
}

/* parse chrome spec */
static struct ChromSpec* parseChromSpec(char* db,
                                        char* chromSpecStr) {
    struct ChromSpec *chromSpec;
    AllocVar(chromSpec);
    if (!hgParseChromRange(db, chromSpecStr,
                           &chromSpec->chrom, &chromSpec->start, &chromSpec->end)) {
        errAbort("invalid chromSpec: %s", chromSpecStr);
    }
    if (chromSpec->start < chromSpec->end) {
        chromSpec->start++; // hgParseChromRange assumes one-based
    }
    return chromSpec;
}

/* parse list chrome spec */
static struct ChromSpec* parseChromSpecs(char* db,
                                         struct slName* chromSpecStrs) {
    struct ChromSpec* chromSpecs = NULL;
    for (struct slName *chromSpecStr = chromSpecStrs; chromSpecStr != NULL; chromSpecStr = chromSpecStr->next) {
        slAddHead(&chromSpecs, parseChromSpec(db, chromSpecStr->name));
    }
    return chromSpecs;
}

/* common code to build sql query with chromSpec */
static struct dyString* makeQuery(char *sqlTemplate, char *table,
                                  char *chromCol, char *startCol, char *endCol,
                                  struct ChromSpec *chromSpec) {
    struct dyString *query = dyStringNew(256);
    dyStringPrintf(query, sqlTemplate, table);
    if (chromSpec != NULL) {
        dyStringPrintf(query, " WHERE (%s = \"%s\")", chromCol, chromSpec->chrom);
        if (chromSpec->start < chromSpec->end) {
            dyStringPrintf(query, " AND ");
            hAddBinToQueryGeneral("bin", chromSpec->start, chromSpec->end, query);
            dyStringPrintf(query, " (%s < %d) AND (%s > %d)", startCol, chromSpec->end, endCol, chromSpec->start);
        }
    }
    return query;
}

/* load PSLs for a range */
static struct psl* loadPslsRange(struct sqlConnection *hgConn, char *table,
                                 struct ChromSpec *chromSpec) {
    char *sqlTemplate =
        "NOSQLINJ SELECT matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, "
        "tBaseInsert, strand, concat(qName,\".\",version), qSize, qStart, qEnd, tName, "
        "tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts "
        "FROM %s LEFT JOIN hgFixed.gbCdnaInfo ON (qName = acc)";
    struct dyString *query = makeQuery(sqlTemplate, table, "tName", "tStart", "tEnd", chromSpec);
    verbose(3, "%s: %s\n", sqlGetDatabase(hgConn), dyStringContents(query));

    struct sqlResult *sr = sqlGetResult(hgConn, dyStringContents(query));
    char **row;
    struct psl *psls = NULL;
    while ((row = sqlNextRow(sr)) != NULL) {
        struct psl* psl = pslLoad(row); 
        slAddHead(&psls, psl);
        pslVerb(3, "load", psl);
    }
    dyStringFree(&query);
    return psls;
}

/* Compare to sort based on target start, but check all fields before
 * declaring unique */
static int pslFullCmpTarget(const void *va, const void *vb) {
#define cmpNumMacro(f1, f2) {int dif = f1 - f2; if (dif != 0) {return dif;}}
#define cmpStrMacro(f1, f2) {int dif = strcmp(f1, f2); if (dif != 0) {return dif;}}
    const struct psl *a = *((struct psl **)va);
    const struct psl *b = *((struct psl **)vb);
    cmpStrMacro(a->tName, b->tName);
    cmpNumMacro(a->tStart, b->tStart);
    cmpNumMacro(a->tEnd, b->tEnd);
    cmpNumMacro(a->tSize, b->tSize);
    cmpStrMacro(a->qName, b->qName);
    cmpNumMacro(a->qSize, b->qSize);
    cmpNumMacro(a->qStart, b->qStart);
    cmpNumMacro(a->qEnd, b->qEnd);
    cmpStrMacro(a->strand, b->strand);
    cmpNumMacro(a->match, b->match);
    cmpNumMacro(a->misMatch, b->misMatch);
    cmpNumMacro(a->repMatch, b->repMatch);
    cmpNumMacro(a->nCount, b->nCount);
    cmpNumMacro(a->qNumInsert, b->qNumInsert);
    cmpNumMacro(a->qBaseInsert, b->qBaseInsert);
    cmpNumMacro(a->tNumInsert, b->tNumInsert);
    cmpNumMacro(a->tBaseInsert, b->tBaseInsert);
    cmpNumMacro(a->blockCount, b->blockCount);
    for (int i = 0; i < a->blockCount; i++) {
        cmpNumMacro(a->blockSizes[i], b->blockSizes[i]);
        cmpNumMacro(a->qStarts[i], b->qStarts[i]);
        cmpNumMacro(a->tStarts[i], b->tStarts[i]);
    }
    return 0;
#undef cmpNumMacro
#undef cmpStrMacro
}

/* load PSLs for all ranges and sort by target */
static struct psl* loadPsls(struct sqlConnection *hgConn, char *table,
                            struct ChromSpec *chromSpecs) {
    struct psl* psls = NULL;
    if (chromSpecs == NULL) {
        psls = loadPslsRange(hgConn, table, NULL);
    } else {
        for (struct ChromSpec *chromSpec = chromSpecs; chromSpec != NULL; chromSpec = chromSpec->next) {
            psls = slCat(psls, loadPslsRange(hgConn, table, chromSpec));
        }
    }
    // both sort and unique to avoid duplication due to an alignment being in two
    // ranges
    slUniqify(&psls, pslFullCmpTarget, pslFree);
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
static void loadEstOrientInfosRange(struct sqlConnection *hgConn, char *table,
                                    struct ChromSpec *chromSpec,
                                    struct hash* orientInfoMap) {
    struct dyString *query = makeQuery("NOSQLINJ SELECT * FROM %s", table, "chrom", "chromStart", "chromEnd", chromSpec);
    struct sqlResult *sr = sqlGetResult(hgConn, dyStringContents(query));
    char **row;
    while ((row = sqlNextRow(sr)) != NULL) {
        struct estOrientInfo *eoi = estOrientInfoLoadLm(row+1, orientInfoMap->lm);
        char *key = getOrientInfoKey(eoi->name, eoi->chrom, eoi->chromStart, eoi->chromEnd);
        hashAdd(orientInfoMap, key, eoi);
    }
    sqlFreeResult(&sr);
    dyStringFree(&query);
}

/* load EST orient info, hashed by name, chrom, chromStart, chromEnd.
 * record are in hash local memory. */
static struct hash* loadEstOrientInfos(struct sqlConnection *hgConn, char *table,
                                struct ChromSpec *chromSpecs) {
    struct hash* orientInfoMap = hashNew(18);
    for (struct ChromSpec *chromSpec = chromSpecs; chromSpec != NULL; chromSpec = chromSpec->next) {
        loadEstOrientInfosRange(hgConn, table, chromSpec, orientInfoMap);
    }
    return orientInfoMap;
}

/* check is a psl looks reversed */
static boolean isPslReversed(struct hash* orientInfoMap, struct psl *psl) {
    char *key = getOrientInfoKey(psl->qName, psl->tName, psl->tStart, psl->tEnd);
    struct estOrientInfo *orientInfo = hashFindVal(orientInfoMap, key);
    return (orientInfo != NULL) && (orientInfo->intronOrientation < 0);
}

/* orient ESTs when possible */
static void orientEstPsls(struct sqlConnection *hgConn, struct psl *psls, struct ChromSpec *chromSpecs) {
    struct hash* orientInfoMap = loadEstOrientInfos(hgConn, "estOrientInfo", chromSpecs);
    for (struct psl *psl = psls; psl != NULL; psl = psl->next) {
        if (isPslReversed(orientInfoMap, psl)) {
            pslRc(psl);
        }
    }
    hashFree(&orientInfoMap);
}

/* load and orient alignments */
static struct psl *loadAligns(char *ucscDb, char *type, struct ChromSpec *chromSpecs) {
    struct sqlConnection *hgConn = sqlConnect(ucscDb);
    bool isEst = sameString(type, "est");
    struct psl *psls = loadPsls(hgConn, (isEst? "all_est" : "all_mrna"), chromSpecs);
    if (isEst) {
        orientEstPsls(hgConn, psls, chromSpecs);
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
    pslVerb(3, "store", psl);
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
    sqliteExBindText(stmt, 20, dyStringContents(blockSizesBuf));
    sqliteExBindText(stmt, 21, dyStringContents(qStartsBuf));
    sqliteExBindText(stmt, 22, dyStringContents(tStartsBuf));

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
static void tslGetUcscRnaAligns(char *ucscDb, char *type, char *sqliteDb, char *sqliteTable,
                                struct ChromSpec *chromSpecs) {
    struct psl *psls = loadAligns(ucscDb, type, chromSpecs);
    storeSqliteDb(psls, sqliteDb, sqliteTable);
    // Don't bother: pslFreeList(&psls);
}

/* entry */
int main(int argc, char** argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 4)
        usage("wrong # args");
    if (!(sameString(argv[2], "rna") || sameString(argv[2], "est"))) {
        usage("expected type of `rna' or `est'");
    }
    struct ChromSpec *chromSpecs = NULL;
    struct slName *chromSpecStrs = optionMultiVal("chromSpec", NULL);
    if (chromSpecStrs != NULL) {
        chromSpecs = parseChromSpecs(argv[1], chromSpecStrs);
    }
    char *sqliteTable = optionVal("table",
                                  (sameString(argv[2], "rna") ? UCSC_RNA_ALN_TBL : UCSC_EST_ALN_TBL));
    tslGetUcscRnaAligns(argv[1], argv[2], argv[3], sqliteTable, chromSpecs);
    return 0;
}
