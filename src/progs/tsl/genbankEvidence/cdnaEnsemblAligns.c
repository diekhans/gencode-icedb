#include "common.h"
#include "options.h"
#include "dnautil.h"
#include "psl.h"
#include "pslTransMap.h"
#include "genomeRangeTree.h"
#include "jksql.h"
#include "hash.h"
#include "linefile.h"
#include "dystring.h"
#include "verbose.h"
#include <stdbool.h>


static char *ensDbHost = "ensembldb.ensembl.org"; 
static char *ensDbUser = "anonymous";
static char *ensDbPass = NULL;
static char* ensDbPort = "5306";

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "cdnaEnsemblAligns [options] ensemblDb cdnaLengthCacheDb mappingPsl outPsl\n"
        "\n"
        "Fetch cDNA alignments from Ensembl and write as PSLs.\n"
        "\n"
        "The Ensembl alignments don't include the poly-A and lack the lengths,\n"
        "The length of cDNA sequences are obtained from either the UCSC browser\n"
        "database or the cdnaExtraLengths file. If lengths can be found, the\n"
        "accessions are added to missingAccFile.\n"
        "\n"
        "The mappingPsl is used to map Ensembl chromosome names and non-reference\n"
        "assembled haplotype UCSC chromosomes and deal with the difference chrM\n"
        "sequences. These are produced by ensToUcscChromMap and the UCSC chrM mappings.\n"
        "\n"
        "Options:\n"
        "  -genomeDb=db - obtain sizes from this genome database\n"
        "  -ensemblPsl=ens.psl - write Ensembl PSL records to this file without mapping to UCSC\n"
        "  -cdnaExtraLengths=tsvFile - obtain sizes from this table file of accesionVersion\n"
        "   and cDNA length\n"
        "  -accverList=file - list of accession/versions to process, used for testing.\n"
        "  -chrom=chr - restrict to this UCSC chrom name, maybe repeated, used for testing\n"
        "  -verbose=n - verbose level\n";
    errAbort("%s:\n%s", msg, usageMsg);
}

static struct optionSpec optionSpecs[] = {
    {"genomeDb", OPTION_STRING},
    {"ensemblPsl", OPTION_STRING},
    {"cdnaExtraLengths", OPTION_STRING},
    {"chrom", OPTION_STRING|OPTION_MULTI},
    {"accverList", OPTION_STRING},
    {NULL, 0}
};
static struct slName *restrictEnsChroms = NULL;

/* sql query to load cDNA alignments from ensembl database
 * Using only alignments that are in the transcript_supporting_feature
 * removes per-exon features.
 */
/* FIXME: ensembl 85-86 did not have external_db set in dna_align_feature */
#define EXTERNAL_DB_HACK 0
static char *ensCDnaAlnQuery =
    "select "
    "  daf.hit_name, daf.hit_strand, daf.hit_start, daf.hit_end, "
    "  sr.name, sr.length, daf.seq_region_strand, daf.seq_region_start, daf.seq_region_end, "
    "  daf.cigar_line "
    "from "
#if EXTERNAL_DB_HACK
    "   dna_align_feature daf, seq_region sr, transcript_supporting_feature tsf "
#else
    "   dna_align_feature daf, seq_region sr, external_db ed, transcript_supporting_feature tsf "
#endif
    "where "
    "((daf.seq_region_id = sr.seq_region_id) "
#if EXTERNAL_DB_HACK
    " and (daf.hit_name not like \"N%\\_%\")"
#else
    " and (daf.external_db_id = ed.external_db_id) "
    " and (ed.db_name = \"EMBL\") "
#endif
    " and (daf.dna_align_feature_id = tsf.feature_id) "
    " and (tsf.feature_type = \"dna_align_feature\"))";


/* load chromosome mappings psl */
static struct genomeRangeTree *loadEnsUcscMappings(char *mappingsPslFile) {
    struct genomeRangeTree *ensUcscMappings = genomeRangeTreeNew();
    struct psl *psl, *psls = pslLoadAll(mappingsPslFile);
    while ((psl = slPopHead(&psls)) != NULL) {
        genomeRangeTreeAddValList(ensUcscMappings, psl->qName, psl->qStart, psl->qEnd, psl);
    }
    return ensUcscMappings;
}

/* convert one UCSC chromosome name to one or more Ensemble names */
static struct slName *convertUcscToEnsChromName(struct genomeRangeTree *ensUcscMappings, char *ucscChrom) {
    struct slName *ensChroms = NULL;
    // must scan all, not indexed by UCSC chrom
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(ensUcscMappings->hash);
    while ((hel = hashNext(&cookie)) != NULL) {
        for (struct range *range = genomeRangeTreeList(ensUcscMappings, hel->name); range != NULL; range = range->next) {
            for (struct psl *psl = range->val; psl != NULL; psl = psl->next) {
                if (sameString(psl->tName, ucscChrom)) {
                    slAddHead(&ensChroms, slNameNew(psl->qName));
                }
            }
        }
    }
    if (ensChroms == NULL) {
        errAbort("can't convert UCSC chrom \"%s\" to an Ensembl chrom", ucscChrom);
    }
    return ensChroms;
}

/* convert a list of ENSEMBL chrom names into UCSC chrom names.  This
 * is slow and only for -chrom argument parsing*/
static struct slName *convertUcscToEnsChromNames(struct genomeRangeTree *ensUcscMappings, struct slName *ucscChroms) {
    struct slName *ensChroms = NULL;
    for (struct slName *ucscChrom = ucscChroms; ucscChrom != NULL; ucscChrom = ucscChrom->next) {
        ensChroms = slCat(ensChroms, convertUcscToEnsChromName(ensUcscMappings, ucscChrom->name));
    }
    return ensChroms;
}

/* map an Ensembl alignment to one a UCSC alignments */
static struct psl *mapEnsToUcsc(struct genomeRangeTree *ensUcscMappings, struct psl *ensPsl) {
    /* assume one only will map, will need to change code if multiple mappings occur */
    struct psl *ucscPsls = NULL;
    for (struct range *range = genomeRangeTreeAllOverlapping(ensUcscMappings, ensPsl->tName, ensPsl->tStart, ensPsl->tEnd); range != NULL; range = range->next) {
        for (struct psl *mapPsl = range->val; mapPsl != NULL; mapPsl = mapPsl->next) {
            slAddHead(&ucscPsls, pslTransMap(pslTransMapNoOpts, ensPsl, mapPsl));
        }
    }
    if (slCount(ucscPsls) > 1) {
        errAbort("BUG: multiple UCSC psls produced by mapping %s:%d-%d", ensPsl->tName, ensPsl->tStart, ensPsl->tEnd);
    }
    return ucscPsls;
}

static char orientToStrand(int orient) {
    return (orient < 0) ? '-' : '+';
}

/* parse the next cigar op */
static bool getNextCigarOp(char **ptr, char *op, int *size) {
    char *str = *ptr;

    if (str == NULL) {
        return false;
    }
    str = skipLeadingSpaces(str);
    if (*str == '\0') {
        return false;
    }
    char *end;
    *size = strtoul(str, &end, 10);
    if (*size == 0) {
        *size = 1;
    }
    str = skipLeadingSpaces(end);
    if (*str == '\0') {
        errAbort("invalid CIGAR: missing op");
    }
    *op = *str++;
    *ptr = str;
    return true;
}

/* create a PSL from an ENSEMBL-style cigar formatted alignment */
static struct psl* pslFromCigar(char *qName, int qSize, int qStart, int qEnd,
                                char *tName, int tSize, int tStart, int tEnd,
                                char *strand, char *cigar) {
    int blocksAlloced = 4;
    struct psl *psl = pslNew(qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, strand, blocksAlloced, 0);

    char cigarSpec[strlen(cigar+1)];  // copy since parsing is destructive
    strcpy(cigarSpec, cigar);
    char *cigarNext = cigarSpec;
    char op;
    int size;
    int qNext = qStart, qBlkEnd = qEnd;
    if (strand[0] == '-') {
        reverseIntRange(&qNext, &qBlkEnd, qSize);
    }
    int tNext = tStart, tBlkEnd = tEnd;
    if (strand[1] == '-') {
        reverseIntRange(&tNext, &tBlkEnd, tSize);
    }

    while (getNextCigarOp(&cigarNext, &op, &size)) {
        switch (op) {
        case 'M': // match or mismatch (gapless aligned block)
            if (psl->blockCount == blocksAlloced)
                pslGrow(psl, &blocksAlloced);
            
            psl->blockSizes[psl->blockCount] = size;
            psl->qStarts[psl->blockCount] = qNext;
            psl->tStarts[psl->blockCount] = tNext;
            psl->blockCount++;
            psl->match += size;
            tNext += size;
            qNext += size;
            break;
        case 'I': // inserted in target
            tNext += size;
            psl->tNumInsert++;
            psl->tBaseInsert += size;
            break;
        case 'D': // deleted from target
            qNext += size;
            psl->qNumInsert++;
            psl->qBaseInsert += size;
            break;

        default:
            errAbort("invalid CIGAR op %c in %s", op, cigar);
        }
    }
    if (qNext != qBlkEnd) {
        errAbort("CIGAR length does not match aligned query range: %s %s", qName, cigar);
    }
    if (tNext != tBlkEnd) {
        errAbort("CIGAR length does not match aligned target range: %s %s", qName, cigar);
    }
    if (psl->strand[1] == '-') {
        pslRc(psl);
    }
    psl->strand[1] = '\0';
    return psl;
}

/* generate list of SQL queries for all cDNAs based on restrictions */
static struct dyString* mkEnsCDnaAlnQuery(char *accver) {
    struct dyString *sql = dyStringNew(0);
    dyStringAppend(sql, ensCDnaAlnQuery);
    if (accver != NULL) {
        dyStringPrintf(sql, " and (daf.hit_name = \"%s\")", accver);
    }
    if (restrictEnsChroms != NULL) {
        dyStringAppend(sql, " and (sr.name in (");
        for (struct slName *ensChrom = restrictEnsChroms; ensChrom != NULL; ensChrom = ensChrom->next) {
            if (ensChrom != restrictEnsChroms) {
                dyStringAppend(sql, ", ");  // not first
            }
            dyStringPrintf(sql, "\"%s\"", ensChrom->name);
        }
        dyStringAppend(sql, "))");
    }
    return sql;

}

/* load sizes into a hash */
static void loadCdnaInfo(struct hash *cdnaSizes, char *cdnaExtraLengths) {
    // columns: accver gi size moddate seqType
    struct lineFile *lf = lineFileOpen(cdnaExtraLengths, TRUE);
    char *row[5];
    lineFileNextRowTab(lf, row, ArraySize(row));  // skip header
    while (lineFileNextRowTab(lf, row, ArraySize(row))) {
        hashAddInt(cdnaSizes, row[0], sqlUnsigned(row[2]));
    }
    lineFileClose(&lf);
}

/* load sizes into a hash */
static void loadFaSizeTab(struct hash *cdnaSizes, char *faSizeTab) {
    // columns: accver size
    struct lineFile *lf = lineFileOpen(faSizeTab, TRUE);
    char *row[2];
    while (lineFileNextRowTab(lf, row, ArraySize(row))) {
        hashAddInt(cdnaSizes, row[0], sqlUnsigned(row[1]));
    }
    lineFileClose(&lf);
}

/* get size of a cDNA or 0 if not in database */
static int getCDnaSizeFromDb(char *accver, struct sqlConnection *genomeDbConn) {
    char acc[strlen(accver)+1];
    strcpy(acc, accver);
    char *dot = strchr(acc, '.');
    if (dot == NULL) {
        errAbort("invalid cDNA identifier, expected acc.version: \"%s\"", accver);
    }
    *dot = '\0';
    char *ver = dot+1;
    char sql[128];
    safef(sql, sizeof(sql), "select size from gbSeq where (acc = \"%s\") and (version = %s)", acc, ver);
    return sqlQuickNum(genomeDbConn, sql);
}

/* get size of a cDNA or 0 if in database or table */
static int getCDnaSize(char *accver, struct sqlConnection *genomeDbConn, struct hash *cdnaSizes) {
    int size = 0;
    if (genomeDbConn != NULL) {
        size = getCDnaSizeFromDb(accver, genomeDbConn);
    }
    if ((size == 0) && (cdnaSizes != NULL)) {
        size = hashIntValDefault(cdnaSizes, accver, 0);
    }
    return size;
}

/* create an Ensembl PSL from the result */
static struct psl *convertResultToEnsPsl(char *hitName, int hitLength, int hitStrand, int hitStart, int hitEnd,
                                         char *seqName, int seqLength, int seqStrand, int seqStart, int seqEnd,
                                         char *cigar) {
    char strand[3] = {orientToStrand(hitStrand), orientToStrand(seqStrand), '\0'};
    return pslFromCigar(hitName, hitLength, hitStart, hitEnd,
                        seqName, seqLength, seqStart, seqEnd,
                        strand, cigar);
}

/* convert one row */
static void convertResult(char *hitName, int hitStrand, int hitStart, int hitEnd,
                          char *seqName, int seqLength, int seqStrand, int seqStart, int seqEnd, char *cigar,
                          struct genomeRangeTree *ensUcscMappings, struct sqlConnection *genomeDbConn, struct hash *cdnaSizes, FILE *pslFh, FILE *missingFh, FILE *ensemblPslFh) {
    int cdnaSize = getCDnaSize(hitName, genomeDbConn, cdnaSizes);
    if (cdnaSize == 0) {
        fprintf(missingFh, "%s\n", hitName);
    } else {
        struct psl *ensPsl = convertResultToEnsPsl(hitName, cdnaSize, hitStrand, hitStart, hitEnd,
                                                   seqName, seqLength, seqStrand, seqStart, seqEnd,
                                                   cigar);
        struct psl *ucscPsl = mapEnsToUcsc(ensUcscMappings, ensPsl);
        if (ucscPsl == NULL) {
            verbose(2, "can't map Ensembl range %s:%d-%d for %s to UCSC\n", seqName, seqStart, seqEnd, hitName);
        } else {
            pslTabOut(ucscPsl, pslFh);
            pslFree(&ucscPsl);
        }
        if (ensemblPslFh != NULL) {
            pslTabOut(ensPsl, ensemblPslFh);
        }
        pslFree(&ensPsl);
    }
}

/* convert one row from the ensembl query */
static void convertRow(char **row, struct genomeRangeTree *ensUcscMappings, struct sqlConnection *genomeDbConn, struct hash *cdnaSizes, FILE *pslFh, FILE *missingFh, FILE *ensemblPslFh) {
    convertResult(row[0], sqlSigned(row[1]), sqlSigned(row[2])-1, sqlSigned(row[3]), 
                  row[4], sqlSigned(row[5]), sqlSigned(row[6]), sqlSigned(row[7])-1, sqlSigned(row[8]), row[9],
                  ensUcscMappings, genomeDbConn, cdnaSizes, pslFh, missingFh, ensemblPslFh);
}

/* convert one */
static void convertOne(char *accver, struct sqlConnection *ensConn, struct genomeRangeTree *ensUcscMappings, struct sqlConnection *genomeDbConn, struct hash *cdnaSizes, FILE *pslFh, FILE *missingFh, FILE *ensemblPslFh) {
    struct dyString *sql = mkEnsCDnaAlnQuery(accver);
    struct sqlResult *sr = sqlGetResult(ensConn, sql->string);
    int cnt = 0;
    char **row;
    while ((row = sqlNextRow(sr)) != NULL) {
        convertRow(row, ensUcscMappings, genomeDbConn, cdnaSizes, pslFh, missingFh, ensemblPslFh);
        cnt++;
    }
    if (cnt == 0) {
        errAbort("no cDNA alignments found for %s", accver);
    }
    sqlFreeResult(&sr);
    dyStringFree(&sql);
}

/* convert some */
static void convertSome(struct sqlConnection *ensConn, struct genomeRangeTree *ensUcscMappings, struct sqlConnection *genomeDbConn, struct hash *cdnaSizes, FILE *pslFh, FILE *missingFh, FILE *ensemblPslFh, struct slName *accvers) {
    for (struct slName *accver = accvers; accver != NULL; accver = accver->next) {
        convertOne(accver->name, ensConn, ensUcscMappings, genomeDbConn, cdnaSizes, pslFh, missingFh, ensemblPslFh);
    }
}

/* convert all */
static void convertAll(struct sqlConnection *ensConn, struct genomeRangeTree *ensUcscMappings, struct sqlConnection *genomeDbConn, struct hash *cdnaSizes, FILE *pslFh, FILE *missingFh, FILE *ensemblPslFh) {
    struct dyString *sql = mkEnsCDnaAlnQuery(NULL);
    verbose(3, "%s\n", sql->string);
    struct sqlResult *sr = sqlGetResult(ensConn, sql->string);
    int cnt = 0;
    char **row;
    while ((row = sqlNextRow(sr)) != NULL) {
        convertRow(row, ensUcscMappings, genomeDbConn, cdnaSizes, pslFh, missingFh, ensemblPslFh);
        cnt++;
    }
    if (cnt == 0) {
        errAbort("no cDNA alignments found");
    }
    sqlFreeResult(&sr);
    dyStringFree(&sql);
}

/* connect to ensembl db */
static struct sqlConnection *ensemblDbConn(char *ensemblDb) {
    // must go through profiles, since connect function doesn't have port
    static char *profileName = "ensemblDb";
    struct slPair* params = NULL;
    slAddHead(&params, slPairNew("name", profileName));
    slAddHead(&params, slPairNew("host", ensDbHost));
    slAddHead(&params, slPairNew("port", ensDbPort));
    slAddHead(&params, slPairNew("user", ensDbUser));
    slAddHead(&params, slPairNew("pass", ensDbPass));
    sqlProfileConfig(params);
    slPairFreeList(&params);
    return sqlConnectProfile(profileName, ensemblDb);
}

/* retrieve the PSLs */
static void cdnaGetEnsemblAligns(char *ensemblDb, struct genomeRangeTree *ensUcscMappings, char *outPsl, char *missingAccFile,
                                 char *genomeDb, char *cdnaExtraLengths, char *faSizeTab, char *accverList, char *ensemblPsl) {
    struct slName *accvers = (accverList != NULL) ? slNameLoadReal(accverList) : NULL;
    struct sqlConnection *ensConn = ensemblDbConn(ensemblDb);
    struct sqlConnection *genomeDbConn = (genomeDb != NULL) ? sqlConnect(genomeDb) : NULL;
    struct hash *cdnaSizes = NULL;
    if ((cdnaExtraLengths != NULL) || (faSizeTab != NULL)) {
        cdnaSizes = hashNew(0);
        if (cdnaExtraLengths != NULL) {
            loadCdnaInfo(cdnaSizes, cdnaExtraLengths);
        }
        if (faSizeTab != NULL) {
            loadFaSizeTab(cdnaSizes, faSizeTab);
        }
    }
    FILE *pslFh = mustOpen(outPsl, "w");
    FILE *missingFh = mustOpen(missingAccFile, "w");
    FILE *ensemblPslFh = (ensemblPsl != NULL) ? mustOpen(ensemblPsl, "w") : NULL;

    if (accvers != NULL) {
        convertSome(ensConn, ensUcscMappings, genomeDbConn, cdnaSizes, pslFh, missingFh, ensemblPslFh, accvers);
    } else {
        convertAll(ensConn, ensUcscMappings, genomeDbConn, cdnaSizes, pslFh, missingFh, ensemblPslFh);
    }

    carefulClose(&ensemblPslFh);
    carefulClose(&missingFh);
    carefulClose(&pslFh);
    sqlDisconnect(&genomeDbConn);
    sqlDisconnect(&ensConn);
}

/* entry */
int main(int argc, char** argv) {
    dnaUtilOpen();
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5)
        usage("wrong # args");
    struct genomeRangeTree *ensUcscMappings = loadEnsUcscMappings(argv[2]);
    restrictEnsChroms = convertUcscToEnsChromNames(ensUcscMappings, optionMultiVal("chrom", NULL));
    if (!(optionExists("genomeDb") || optionExists("cdnaExtraLengths") || optionExists("faSizeTab"))) {
        errAbort("must specify at least one of -genomeDb, -cdnaExtraLengths, or -faSizeTab");
    }
    cdnaGetEnsemblAligns(argv[1], ensUcscMappings, argv[3], argv[4],
                         optionVal("genomeDb", NULL),
                         optionVal("cdnaExtraLengths", NULL),
                         optionVal("faSizeTab", NULL),
                         optionVal("accverList", NULL),
                         optionVal("ensemblPsl", NULL));
    return 0;
}
