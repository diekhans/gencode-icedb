#include "common.h"
#include "options.h"
#include "psl.h"
#include "hash.h"
#include "hdb.h"
#include "jksql.h"
#include "estOrientInfo.h"
#include "verbose.h"
#include "dystring.h"

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "tslGetUcscRnaAligns ucscDb type pslFile\n"
        "\n"
        "Load PSL alignments from UCSC all_mrna or all_est tables and write to\n"
        "sorted PSL file for indexing.  EST PSLs will be reverse-complement if\n"
        "the estOrientInfo table indicates it is a 3' EST.  Type is `rna' or `est'.\n"
        "\n"
        "Options:\n"
        "  -verbose=n\n"
        "  -table=tbl - load this table, defaults to ucsc_rna_aln or ucsc_est_aln.\n"
        "  -chromSpec=spec - restrict to this chrom or chrom range, for testing.\n"
        "   Maybe repeated. Duplicates caused by alignments being in multiple ranges\n"
        "   are discarded.\n"
        ;
    errAbort("%s:\n%s", msg, usageMsg);
}
static struct optionSpec optionSpecs[] = {
    {"chromSpec", OPTION_STRING|OPTION_MULTI},
    {NULL, 0}
};

static int noOrientInfoCnt = 0;

/* chromosome spec */
struct ChromSpec {
    struct ChromSpec *next;
    char* chrom;    // NULL if not restriction
    int start;      // zero length for whole chromosome
    int end;
};

/* a verbose message about a psl */
static void pslVerb(int level, char* msg, struct psl* psl) {
    verbose(level, "%s: %s:%d-%d <=> %s:%d-%d (%s) blks: %d\n", msg,
            psl->qName, psl->qStart, psl->qEnd,
            psl->tName, psl->tStart, psl->tEnd,
            psl->strand, psl->blockCount);
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
    char nameNoVer[64];
    safecpy(nameNoVer, sizeof(nameNoVer), name);
    char *dot = strchr(nameNoVer, '.');
    if (dot != NULL) {
        *dot = '\0';
    }
    safef(key, sizeof(key), "%s@%s:%d-%d", nameNoVer, chrom, chromStart, chromEnd);
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
        verbose(4, "load orientInfo: %s\n", key);
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
    if (chromSpecs != NULL) {
        for (struct ChromSpec *chromSpec = chromSpecs; chromSpec != NULL; chromSpec = chromSpec->next) {
            loadEstOrientInfosRange(hgConn, table, chromSpec, orientInfoMap);
        }
    } else {
        loadEstOrientInfosRange(hgConn, table, NULL, orientInfoMap);
    }
    return orientInfoMap;
}

/* check is a psl looks reversed */
static boolean isPslReversed(struct hash* orientInfoMap, struct psl *psl) {
    char *key = getOrientInfoKey(psl->qName, psl->tName, psl->tStart, psl->tEnd);
    struct estOrientInfo *orientInfo = hashFindVal(orientInfoMap, key);
    if (orientInfo != NULL) {
        verbose(3, "isPslReversed: %s %d\n", psl->qName, orientInfo->intronOrientation);
    } else {
        noOrientInfoCnt++;
        verbose(3, "isPslReversed: %s no orientInfo (%s)\n", psl->qName, key);
    }
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

/* get rna or est alignments */
static void tslGetUcscRnaAligns(char *ucscDb, char *type, char *pslFile,
                                struct ChromSpec *chromSpecs) {
    struct psl *psls = loadAligns(ucscDb, type, chromSpecs);
    pslWriteAll(psls, pslFile, FALSE);
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
    tslGetUcscRnaAligns(argv[1], argv[2], argv[3], chromSpecs);
    if (noOrientInfoCnt > 0) {
        fprintf(stderr, "WARNING: %d orientInfo records not found", noOrientInfoCnt);
    }
    return 0;
}
