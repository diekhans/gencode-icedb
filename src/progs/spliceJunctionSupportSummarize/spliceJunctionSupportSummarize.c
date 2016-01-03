#include "common.h"
#include "options.h"
#include "intronMap.h"
#include "starSpliceJunction.h"
#include "genePred.h"

static struct optionSpec optionSpecs[] = {
    {"countsReport", OPTION_BOOLEAN},
    {"minOverhang", OPTION_INT},
    {"minNumUniqueMapReads", OPTION_INT},
    {"minNumMultiMapReads", OPTION_INT},
    {NULL, 0}
};

static bool gCountsReport = FALSE;
static int gMinOverhang = 0;
static int gMinNumUniqueMapReads = 0;
static int gMinNumMultiMapReads = 0;

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "spliceJunctionSupportSummarize gencodeGenePred starSpliceJunctionList reportTsv\n\n"
        "Summarize splice junction support\n"
        "\n"
        "Options:\n"
        "   -countsReport - output counts rather than support levels. The -minOverhang\n"
        "        filter is applied.\n"
        "   -minOverhang=n - minimum overhang for a STAR splice junction call.\n"
        "        records with less than this maximum overhang are discards\n"
        "        and not considered more.\n"
        "   -minNumUniqueMapReads=n \n"
        "   -minNumMultiMapReads=n \n";
    errAbort("%s:\n%s", msg, usageMsg);
}

/* load intron map from files */
static struct intronMap* loadIntronMap(char* gencodeGenePred,
                                       char* starSpliceJunctionList) {
    struct intronMap* intronMap = intronMapNew();
    struct slName *spliceJuncFiles = slNameLoadReal(starSpliceJunctionList);
    intronMapLoadTranscripts(intronMap, gencodeGenePred);
    for (struct slName *spliceJuncFile = spliceJuncFiles; spliceJuncFile != NULL; spliceJuncFile = spliceJuncFile->next) {
        intronMapLoadStarJuncs(intronMap, spliceJuncFile->name, gMinOverhang);
    }
    slFreeList(&spliceJuncFiles);
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

/* check if the is novel? */
static int getIsNovel(struct intronInfo* intronInfo) {
    return (intronInfo->mappingsSum != NULL)
        && (intronInfo->mappingsSum->annotated == 0);
}

/* convert intron motif code to string */
static char* getIntronMotifStr(struct intronInfo* intronInfo) {
    if (intronInfo->mappingsSum == NULL) {
        return "unk"; 
    }
    switch (intronInfo->mappingsSum->intronMotif) {
        case 0: return "ncon"; 
        case 1: return "GT/AG";
        case 2: return "CT/AC";
        case 3: return "GC/AG";
        case 4: return "CT/GC";
        case 5: return "AT/AC";
        case 6: return "GT/AT";
        default:
            errAbort("unknown intron motif code: %d", intronInfo->mappingsSum->intronMotif);
            return 0;
    }
}

/* determine the support level */
static int getSupportLevel(struct intronInfo* intronInfo) {
    // level 1: strongest support
    // level 2: at least 
    // level 3: at least 
    // level 4: at least 
    // level 5: no support
    return 0;
}

/* write support report header */
static void reportSupportHeader(FILE* reportFh) {
    static char* header = "chrom\t" "intronStart\t" "intronEnd\t" "novel\t"
        "annotStrand\t" "rnaSeqStrand\t" "intronMotif\t" "supportLevel\t"
        "transcripts\n";
    fputs(header, reportFh);
}

/* report on an intron */
static void reportSupportIntron(struct intronInfo* intronInfo,
                                FILE* reportFh) {
    fprintf(reportFh, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t",
            intronInfo->chrom, intronInfo->chromStart,
            intronInfo->chromEnd,
            getIsNovel(intronInfo),
            getAnnotStrand(intronInfo),
            getRnaSeqStrand(intronInfo),
            getIntronMotifStr(intronInfo),
            getSupportLevel(intronInfo));
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
static void reportSupport(struct intronMap* intronMap,
                          FILE* reportFh) {
    reportSupportHeader(reportFh);
    struct intronInfo* intronInfos = intronMapGetSorted(intronMap);
    for (struct intronInfo* intronInfo = intronInfos; intronInfo != NULL; intronInfo = intronInfo->next) {
        reportSupportIntron(intronInfo, reportFh);
    }
}
#ifdef FIXME
/* write counts report header */
static void reportCountsHeader(FILE* reportFh) {
    static char* header =  "novel\t" "intronMotif\t"
        "intronCount\t" "numUniqueMapReads\t" "numMultiMapReads\t"
        "transcriptCount\n";
    fputs(header, reportFh);
}
#endif

/* main */
static void spliceJunctionSupportSummarize(char* gencodeGenePred,
                                           char* starSpliceJunctionList,
                                           char* reportTsv) {
    struct intronMap* intronMap = loadIntronMap(gencodeGenePred, starSpliceJunctionList);
    FILE* reportFh = mustOpen(reportTsv, "w");
    if (gCountsReport) {
    } else {
        reportSupport(intronMap, reportFh);
    }
    carefulClose(&reportFh);
    intronMapFree(&intronMap);
}

/* entry */
int main(int argc, char** argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 4) {
        usage("wrong # args");
    }
    gCountsReport = optionExists("countsReport");
    gMinOverhang = optionInt("minOverhang", gMinOverhang);
    gMinNumUniqueMapReads = optionInt("minNumUniqueMapReads", gMinNumUniqueMapReads);
    gMinNumMultiMapReads = optionInt("minNumMultiMapReads", gMinNumMultiMapReads);

    spliceJunctionSupportSummarize(argv[1], argv[2], argv[3]);
    return 0;
}

