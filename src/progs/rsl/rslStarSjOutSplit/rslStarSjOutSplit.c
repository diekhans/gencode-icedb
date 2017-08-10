#include "common.h"
#include "options.h"
#include "starResultsDir.h"
#include "starSpliceJunction.h"
#include "starOps.h"
#include "linefile.h"
#include "portable.h"

#define CHROM_MAX 512

static struct optionSpec optionSpecs[] = {
    {"startDirRec", OPTION_INT},
    {"endDirRec", OPTION_INT},
    {"minOverhang", OPTION_INT},
    {NULL, 0}
};

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "rslStarSjOutSplit starResultsDirTsv chromOutDir\n\n"
        "Split STAR sjout files into per-chromosome files, converting to zero-based\n"
        "coordinates and adding the mapping_symid. Output files are named in the form:\n"
        "  chromOutDir/chrom/mapping_symid.sjsup\n"
        "\n"
        "The output files must not exist.\n"
        "\n"
        "  o starResultsDirTsv is TSV file with the columns:\n"
        "       run_acc mapping_param_symid mapping_symid sjout\n"
        "    with file paths relative to location of list file.\n"
        "    A file header is skipped, but not used; columns must be in this order\n"
        "Options:\n"
        "   -startDirRec=n - Zero based record number of the first line in starResultsDirTsv\n"
        "    to process.  The TSV header is not include.  If not specified, process all lines.\n"
        "    Used to parallelize splitting.\n"
        "   -endDirRec=n - Half-open record number of the last line in starResultsDirTsv\n"
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
            "%s\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%u\t%s\n",
            sjout->chrom,
            sjout->chromStart,
            sjout->chromEnd,
            starStrandCodeToChar(sjout->strand),
            starMotifCodeToStr(sjout->intronMotif),
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
    lineFileClose(&lf);
    carefulClose(&chromOutFh);
}

/* main */
static void rslStarSjOutSplit(char* starResultsDirTsv,
                              char* chromOutDir,
                              int minOverhang,
                              int startDirRec,
                              int endDirRec) {
    struct starResults *starResultsDir = starResultsDirLoad(starResultsDirTsv);
    if (startDirRec < 0) {
        // default to whole file
        startDirRec = 0;
        endDirRec = slCount(starResultsDir);
    }
    makeDirsOnPath(chromOutDir);

    int iDirLine = 0;
    for (struct starResults *starResults = starResultsDir; starResults != NULL; starResults = starResults->next) {
        if ((startDirRec <= iDirLine) && (iDirLine < endDirRec)) {
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
    if (optionExists("startDirRec") != optionExists("endDirRec")) {
        errAbort("must specify either both or neither of -startDirRec and -endDirRec");
    }
    int startDirRec = optionInt("startDirRec", -1);
    int endDirRec = optionInt("endDirRec", -1);

    rslStarSjOutSplit(argv[1], argv[2], minOverhang, startDirRec, endDirRec);
    return 0;
}

