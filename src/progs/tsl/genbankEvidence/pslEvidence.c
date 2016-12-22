/* various functions related to psl evidence */
#include "common.h"
#include "pslEvidence.h"
#include "psl.h"
#include "dnautil.h"
#include <stdbool.h>

struct Range {
    int start;
    int end;
};

/* Fill in a range array of Ranges, in chromosome coordinates, with small gaps
 * removed. Return the number of entries used in the array, which should be
 * large enough to hold all PSL blocks. */
static int getBlockRanges(struct psl *psl, int maxIgnoreTGapSize, struct Range *ranges) {
    if (pslTStrand(psl) == '-') {
        errAbort("can't handle negative target strand");
    }
    int nRanges = 0;
    for (int iBlk = 0; iBlk < psl->blockCount; iBlk++) {
        if ((nRanges == 0) || ((ranges[nRanges-1].end - psl->tStarts[iBlk]) >= maxIgnoreTGapSize)) {
            ranges[nRanges].start = psl->tStarts[iBlk];
            ranges[nRanges].end = pslTEnd(psl, iBlk);
            nRanges++;
        } else {
            ranges[nRanges].end = pslTEnd(psl, iBlk);
        }
    }
    return nRanges;
}

/* compare two range arrays of the same size, allowing for fuzzy ends */
static bool cmpRanges(struct Range *ranges1, struct Range *ranges2, int nRanges) {
    for (int iRange = 0; iRange < nRanges; iRange++) {
        if (((iRange > 0) && (ranges1[iRange].start != ranges2[iRange].start))
            || ((iRange < nRanges-1) && (ranges1[iRange].end != ranges2[iRange].end))) {
            return false;
        }
    }
    return true;
}

/* Do the target gaps in two psl alignments similar, regardless of strand?
 * this ignores target gaps up to the specified size. */
bool pslTGapsSimilar(struct psl *psl1, struct psl *psl2, int maxIgnoreTGapSize) {
    if (!sameString(psl1->tName, psl2->tName)) {
        return false;
    }
    struct Range *ranges1 = alloca(psl1->blockCount*sizeof(struct Range));
    int nRanges1 =  getBlockRanges(psl1, maxIgnoreTGapSize, ranges1);
    struct Range *ranges2 = alloca(psl2->blockCount*sizeof(struct Range));
    int nRanges2 = getBlockRanges(psl2, maxIgnoreTGapSize, ranges2);
    if (nRanges1 != nRanges2) {
        return false;
    }

    return cmpRanges(ranges1, ranges2, nRanges1);
}

/* Compare to sort based on target min start and max end. */
static int pslCmpTargetMax(const void *va, const void *vb) {
    const struct psl *a = *((struct psl **)va);
    const struct psl *b = *((struct psl **)vb);
    int dif = strcmp(a->tName, b->tName);
    if (dif == 0)
        dif = a->tStart - b->tStart;
    if (dif == 0)
        dif = b->tEnd - a->tEnd;
    return dif;
}

/* should this psl be used? */
static bool inclPsl(struct psl *psl, pslFilterFunc filterFunc) {
    return (((filterFunc == NULL) || filterFunc(psl)))
        && (pslCheck(NULL, NULL, psl) == 0);
}

/* utility to load all psl, dropping invalid ones and sort by target */
struct psl* pslEvidenceLoad(char *pslFile, pslFilterFunc filterFunc) {
    struct psl *psls = NULL, *psl;
    struct lineFile *fh = pslFileOpen(pslFile);
    while ((psl = pslNext(fh)) != NULL) {
        if (inclPsl(psl, filterFunc)) {
            slAddHead(&psls, psl);
        } else {
            pslFree(&psl);
        }
    }
    lineFileClose(&fh);
    slSort(&psls, pslCmpTargetMax);
    return psls;
}

