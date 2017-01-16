#include "common.h"
#include "rslAnalysisSet.h"
#include "filePath.h"
#include "lineFile.h"

/* construct new analysis object */
static struct rslAnalysis *rslAnalysisNew(char *runname, char *tissue, char *sjPath) {
    struct rslAnalysis *rslAnalysis;
    AllocVar(rslAnalysis);
    rslAnalysis->runname = cloneString(runname);
    rslAnalysis->tissue = cloneString(tissue);
    rslAnalysis->sjPath = cloneString(sjPath);
    return rslAnalysis;
}

/* freee an analysis object */
static void rslAnalysisFree(struct rslAnalysis *rslAnalysis) {
    freeMem(rslAnalysis->runname);
    freeMem(rslAnalysis->tissue);
    freeMem(rslAnalysis->sjPath);
    freeMem(rslAnalysis);
}

/* construct an RSL analysis for a geneset */
static struct rslAnalysisSet *rslAnalysisSetNew(char *setname) {
    struct rslAnalysisSet *rslAnalysisSet;
    AllocVar(rslAnalysisSet);
    rslAnalysisSet->setname = cloneString(setname);
    return rslAnalysisSet;
}

/* free an RSL analysis geneset */
void rslAnalysisSetFree(struct rslAnalysisSet *rslAnalysisSet) {
    struct rslAnalysis *rslAnalysis;
    while ((rslAnalysis = slPopHead(&rslAnalysisSet->analyses))) {
        rslAnalysisFree(rslAnalysis);
    }
    freeMem(rslAnalysisSet->setname);
    freeMem(rslAnalysisSet);
}

/* load one RSL analysis */
static void loadRslAnalysis(struct rslAnalysisSet *rslAnalysisSet, char *tsvFile,
                            char *runname, char *tissue, char *relSjPath) {
    char *sjPath = pathRelativeToFile(tsvFile, relSjPath);
    struct rslAnalysis *rslAnalysis= rslAnalysisNew(runname, tissue, sjPath);
    freeMem(sjPath);  // was copied
    slAddHead(&rslAnalysisSet->analyses, rslAnalysis);
}

/* load an RSL analysis set from a sqlite3 database plus output files */
struct rslAnalysisSet *rslAnalysisSetLoad(char *tsvFile,
                                          char *setname) {
    struct rslAnalysisSet *rslAnalysisSet = rslAnalysisSetNew(setname);
    struct lineFile *lf = lineFileOpen(tsvFile, TRUE);
    // skip header
    static const int numCols = 3;
    char *row[numCols];
    lineFileNextRowTab(lf, row, numCols);
    while (lineFileNextRowTab(lf, row, numCols)) {
        loadRslAnalysis(rslAnalysisSet, tsvFile, row[0], row[1], row[2]);
    }
    slReverse(&rslAnalysisSet->analyses);
    return rslAnalysisSet;
}
