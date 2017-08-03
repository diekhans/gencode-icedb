#include "common.h"
#include "options.h"
#include "starResultsDir.h"
#include "starSpliceJunction.h"
#include "linefile.h"
#include "portable.h"

#define CHROM_MAX 512

static struct optionSpec optionSpecs[] = {
    {"startDirLine", OPTION_INT},
    {"endDirLine", OPTION_INT},
    {"minOverhang", OPTION_INT},
    {NULL, 0}
};

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "rslStarSjOutSplit starResultsDirTsv chromOutDir\n\n"
        "Split STAR sjout files into per-chromosome files, converting to zero-based\n"
        "coordinates and adding the mapping_symid. Output files are named in the form:\n"
        "  chromOutDir/chrom.mapping_symid.sjsup\n"
        "\n"
        "The output files must not exist.\n"
        "\n"
        "  o starResultsDirTsv is TSV file with the columns:\n"
        "       run_acc mapping_param_symid mapping_symid sjout\n"
        "    with file paths relative to location of list file.\n"
        "    A file header is skipped, but not used; columns must be in this order\n"
        "Options:\n"
        "   -startDirLine=n - Zero based line number of the first line in starResultsDirTsv\n"
        "    to process.  The TSV header is line zero.  If not specified, process all lines.\n"
        "    Used to parallelize splitting.\n"
        "   -endDirLine=n - Half-open line number of the last line in starResultsDirTsv\n"
        "    to process.\n"
        "   -minOverhang=n - minimum overhang for a STAR splice junction to include.\n"
        ;
    errAbort("%s:\n%s", msg, usageMsg);
}

/* open new chrom file, error if it exists. Close current one if open.  */
static void chromSplitFile(char *chromOutDir,
                           char* mapping_symid,
                           char* chrom,
                           FILE **chromOutFhPtr) {
    char chromDir[PATH_LEN], path[PATH_LEN];
    safef(chromDir, sizeof(chromDir), "%s/%s", chromOutDir, chrom);
    safef(path, sizeof(path), "%s/%s.sjsup", chromDir, mapping_symid);
    if (fileExists(path)) {
        errAbort("chrom split output file already exists, either output directory was not empty or rows not sorted by chrom: %s",
                 path);
    }
    carefulClose(chromOutFhPtr);
    makeDir(chromDir);
    *chromOutFhPtr = mustOpen(path, "w");
}

/* write one record, opening new chromOut if needed */
static void starSjSupportWrite(char *chromOutDir,
                               char* mapping_symid,
                               char currentChrom[CHROM_MAX],
                               FILE **chromOutFhPtr,
                               struct starSpliceJunction *sjout) {
    if (!sameString(sjout->chrom, currentChrom)) {
        chromSplitFile(chromOutDir, mapping_symid, sjout->chrom, chromOutFhPtr);
        safecpy(currentChrom, CHROM_MAX, sjout->chrom);
    }
    fprintf(*chromOutFhPtr,
            "%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\n",
            sjout->chrom,
            sjout->chromStart,
            sjout->chromEnd,
            sjout->strand,
            sjout->intronMotif,
            sjout->annotated,
            sjout->numUniqueMapReads,
            sjout->numMultiMapReads,
            sjout->maxOverhang,
            mapping_symid);
}

/* split a sjout file, assumes grouped by chrom  */
static void starSjOutSplit(struct starResults *starResults,
                            int minOverhang,
                            char *chromOutDir) {
    FILE *chromOutFh = NULL;
    char currentChrom[CHROM_MAX];
    ZeroVar(currentChrom);
    char* sjoutRow[STARSPLICEJUNCTION_NUM_COLS];

    struct lineFile *lf = lineFileOpen(starResults->sjout, TRUE);
    while (lineFileNextRowTab(lf, sjoutRow, STARSPLICEJUNCTION_NUM_COLS)) {
        struct starSpliceJunction *sjout = starSpliceJunctionLoad(sjoutRow);
        if (sjout->maxOverhang >= minOverhang) {
            starSjSupportWrite(chromOutDir, starResults->mapping_symid,
                               currentChrom, &chromOutFh, sjout);
        }
        starSpliceJunctionFree(&sjout);
    }
    carefulClose(&chromOutFh);
}

/* main */
static void rslStarSjOutSplit(char* starResultsDirTsv,
                              char* chromOutDir,
                              int minOverhang,
                              int startDirLine,
                              int endDirLine) {
    struct starResults *starResultsDir = starResultsDirLoad(starResultsDirTsv);
    if (startDirLine < 0) {
        // default to whole file
        startDirLine = 1;  // skip header
        endDirLine = slCount(starResultsDir) + 1;
    }
    makeDir(chromOutDir);

    int iDirLine = 1; // header already skipped
    for (struct starResults *starResults = starResultsDir; starResults != NULL; starResults = starResults->next) {
        if ((startDirLine <= iDirLine) && (iDirLine < endDirLine)) {
            starSjOutSplit(starResults, minOverhang, chromOutDir);
        }
        iDirLine++;
    }
    starResultsDirFree(starResultsDir);
}

/* entry */
int main(int argc, char** argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 3) {
        usage("wrong # args");
    }
    int minOverhang = optionInt("minOverhang", 0);
    if (optionExists("startDirLine") != optionExists("endDirLine")) {
        errAbort("must specify either both or neither of -startDirLine and -endDirLine");
    }
    int startDirLine = optionInt("startDirLine", -1);
    int endDirLine = optionInt("endDirLine", -1);

    rslStarSjOutSplit(argv[1], argv[2], minOverhang, startDirLine, endDirLine);
    return 0;
}

