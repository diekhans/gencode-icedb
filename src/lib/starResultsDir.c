#include "common.h"
#include "starResultsDir.h"
#include "filePath.h"
#include "linefile.h"

/* construct start results dir object */
static struct starResults *starResultsNew(char *starResultsDirTsv,
                                          char *run_acc,
                                          char *mapping_param_symid,
                                          char *mapping_symid,
                                          char *sjout) {
    struct starResults *starResults;
    AllocVar(starResults);
    starResults->run_acc = cloneString(run_acc);
    starResults->mapping_param_symid = cloneString(mapping_param_symid);
    starResults->mapping_symid = cloneString(mapping_symid);
    starResults->sjout = pathRelativeToFile(starResultsDirTsv, sjout);
    return starResults;
}

/* freee an analysis object */
static void starResultsFree(struct starResults *starResults) {
    freeMem(starResults->run_acc);
    freeMem(starResults->mapping_param_symid);
    freeMem(starResults->mapping_symid);
    freeMem(starResults->sjout);
    freeMem(starResults);
}

/* load an STAR results directory TSV files */
struct starResults *starResultsDirLoad(char *starResultsDirTsv) {
    struct starResults* starResultsDir = NULL;
    static const int numCols = 4;
    char *row[numCols];
    struct lineFile *lf = lineFileOpen(starResultsDirTsv, TRUE);
    lineFileNextRowTab(lf, row, numCols);     // skip header

    while (lineFileNextRowTab(lf, row, numCols)) {
        slAddHead(&starResultsDir,
                  starResultsNew(starResultsDirTsv, row[0], row[1], row[2], row[3]));
    }
    slReverse(&starResultsDir);
    lineFileClose(&lf);
    return starResultsDir;
}

/* free STAR results dir list */
void starResultsDirFree(struct starResults *starResultsDir) {
    struct starResults* starResults;
    while((starResults = slPopHead(&starResultsDir)) != NULL) {
        starResultsFree(starResults);
    }
}
