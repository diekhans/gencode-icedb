#include "common.h"
#include "options.h"
#include "intronMap.h"
#include "starSpliceJunction.h"
#include "genePred.h"
#include "filePath.h"
#include "intronCounts.h"

static struct optionSpec optionSpecs[] = {
    {"countsReport", OPTION_BOOLEAN},
    {"minOverhang", OPTION_INT},
    {"minNumUniqueMapReads", OPTION_INT},
    {"minNumMultiMapReads", OPTION_INT},
    {NULL, 0}
};

static int gMinOverhang = 0;
static int gMinNumUniqueMapReads = 0;
static int gMinNumMultiMapReads = 0;

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "spliceJunctionCollectEvidence gencodeGenePred gencodeSpliceTsv starSpliceJunctionList reportTsv\n\n"
        "Collect splice junction supporting evidence\n"
        "\n"
        "  o starSpliceJunctionList is list of splice junction files,\n"
        "    with file names relative to location of list file.\n"
        "Options:\n"
        "   -minOverhang=n - minimum overhang for a STAR splice junction call.\n"
        "        records with less than this maximum overhang are discards\n"
        "        and not considered more.\n"
        "   -minNumUniqueMapReads=n \n"
        "   -minNumMultiMapReads=n \n";
    errAbort("%s:\n%s", msg, usageMsg);
}

/* load intron map from files */
static struct intronMap* loadIntronMap(char* gencodeGenePred,
                                       char* gencodeSpliceTsv,
                                       char* starSpliceJunctionList) {
    struct intronMap* intronMap = intronMapNew();
    struct slName *spliceJuncFiles = slNameLoadReal(starSpliceJunctionList);
    intronMapLoadTranscripts(intronMap, gencodeGenePred);
    for (struct slName *spliceJuncFile = spliceJuncFiles; spliceJuncFile != NULL; spliceJuncFile = spliceJuncFile->next) {
        char* juncPath = pathRelativeToFile(starSpliceJunctionList, spliceJuncFile->name);
        intronMapLoadStarJuncs(intronMap, juncPath, gMinOverhang);
        freeMem(juncPath);
    }
    slFreeList(&spliceJuncFiles);
    intronMapLoadTranscriptSpliceSites(intronMap, gencodeSpliceTsv);
    return intronMap;
}

/* get annotation strand (WARNING: static return) */
static char* getAnnotStrand(struct intronInfo* intronInfo) {
    static char buf[64];
    buf[0] = '\0';
    // report multi strands if annotations conflict
    for (struct intronTransLink* intronTrans = intronInfo->intronTranses;
         intronTrans != NULL; intronTrans = intronTrans->next) {
        if (strstr(buf, intronTrans->transcript->strand) == NULL) {
            if (strlen(buf) > 0) {
                safecat(buf, sizeof(buf), "/");
            }
            safecat(buf, sizeof(buf), intronTrans->transcript->strand);
        }
    }
    return buf;
}

/* get RNA-Seq strand (WARNING: static return) */
static char* getRnaSeqStrand(struct intronInfo* intronInfo) {
    if (intronInfo->mappingsSum == NULL) {
        return "";
    }
    switch (intronInfo->mappingsSum->strand) {
        case 0: return "?";
        case 1: return "+";
        case 2: return "-";
        default:
            errAbort("invalid RNA-Seq strand code: %d", intronInfo->mappingsSum->strand);
            return NULL;
    }
}

/* write report header */
static void reportEvidenceHeader(FILE* reportFh) {
    static char* header = "chrom\t" "intronStart\t" "intronEnd\t" "novel\t"
        "annotStrand\t" "rnaSeqStrand\t" "intronMotif\t"
        "numUniqueMapReads\t" "numMultiMapReads\t"
        "transcripts\n";
    fputs(header, reportFh);
}

/* report on an intron */
static void reportEvidenceIntron(struct intronInfo* intronInfo,
                                 FILE* reportFh) {
    int numUniqueMapReads = intronInfo->mappingsSum != NULL ? intronInfo->mappingsSum->numUniqueMapReads : 0;
    int numMultiMapReads = intronInfo->mappingsSum != NULL ? intronInfo->mappingsSum->numMultiMapReads : 0;
    fprintf(reportFh, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t",
            intronInfo->chrom, intronInfo->chromStart,
            intronInfo->chromEnd,
            intronInfoIsNovel(intronInfo),
            getAnnotStrand(intronInfo),
            getRnaSeqStrand(intronInfo),
            intronInfoMotifStr(intronInfo),
            numUniqueMapReads, numMultiMapReads);
    for (struct intronTransLink* intronTrans = intronInfo->intronTranses;
         intronTrans != NULL; intronTrans = intronTrans->next) {
        if (intronTrans != intronInfo->intronTranses) {
            fputc(',', reportFh);
        }
        fputs(intronTrans->transcript->name, reportFh);
    }
    fputc('\n', reportFh);
}

/* report on results */
static void reportEvidence(struct intronMap* intronMap,
                           FILE* reportFh) {
    reportEvidenceHeader(reportFh);
    for (struct intronInfo* intronInfo = intronMapGetSorted(intronMap); intronInfo != NULL; intronInfo = intronInfo->next) {
        reportEvidenceIntron(intronInfo, reportFh);
    }
}

/* main */
static void spliceJunctionCollectEvidence(char* gencodeGenePred,
                                          char* gencodeSpliceTsv,
                                          char* starSpliceJunctionList,
                                          char* reportTsv) {
    struct intronMap* intronMap = loadIntronMap(gencodeGenePred, gencodeSpliceTsv, starSpliceJunctionList);
    FILE* reportFh = mustOpen(reportTsv, "w");
    reportEvidence(intronMap, reportFh);
    carefulClose(&reportFh);
    intronMapFree(&intronMap);
}

/* entry */
int main(int argc, char** argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 5) {
        usage("wrong # args");
    }
    gMinOverhang = optionInt("minOverhang", gMinOverhang);
    gMinNumUniqueMapReads = optionInt("minNumUniqueMapReads", gMinNumUniqueMapReads);
    gMinNumMultiMapReads = optionInt("minNumMultiMapReads", gMinNumMultiMapReads);

    spliceJunctionCollectEvidence(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}

