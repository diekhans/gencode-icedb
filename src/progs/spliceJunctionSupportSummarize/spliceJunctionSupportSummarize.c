#include "common.h"
#include "options.h"
#include "intronMap.h"

static struct optionSpec optionSpecs[] = {
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
    static char* usageMsg = "spliceJunctionSupportSummarize gencodeGenePred starSpliceJunctionList summaryTsv\n\n"
        "Summarize splice junction support\n"
        "\n"
        "\n"
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

/* report on results */
static void report(struct intronMap* intronMap,
                   FILE* summaryFh) {
    static char* header = "transcript\t"
        "chrom\t" "intronStart\t" "intronEnd\t" "strand\t"
        ;
    fputs(header, summaryFh);
}

/* main */
static void spliceJunctionSupportSummarize(char* gencodeGenePred,
                                           char* starSpliceJunctionList,
                                           char* summaryTsv) {
    struct intronMap* intronMap = loadIntronMap(gencodeGenePred, starSpliceJunctionList);
    FILE* summaryFh = mustOpen(summaryTsv, "w");
    report(intronMap, summaryFh);
    carefulClose(&summaryFh);
    intronMapFree(&intronMap);
}

/* entry */
int main(int argc, char** argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 4) {
        usage("wrong # args");
    }
    gMinOverhang = optionInt("minOverhang", gMinOverhang);
    gMinNumUniqueMapReads = optionInt("minNumUniqueMapReads", gMinNumUniqueMapReads);
    gMinNumMultiMapReads = optionInt("minNumMultiMapReads", gMinNumMultiMapReads);

    spliceJunctionSupportSummarize(argv[1], argv[2], argv[3]);
    return 0;
}

