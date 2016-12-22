/* various functions related to psl evidence */
#ifndef pslEvidence_h
#define pslEvidence_h
#include "psl.h"
#include <stdbool.h>

/* assumed minimum size of an intron */
static const int pslEvidenceMinIntronSize = 30;

/* get size of target gap before block */
INLINE unsigned pslTGapSize(struct psl *psl, int blkIdx) {
    return pslTStart(psl, blkIdx) - pslTEnd(psl, blkIdx-1);
}

/* does a psl overlap the genomic range? */
INLINE bool pslTOverlap(struct psl *psl1, struct psl *psl2) {
    return sameString(psl1->tName, psl2->tName)
        && (psl1->tStart < psl2->tEnd) && (psl1->tEnd > psl2->tStart);
}

/* Do the target gaps in two psl alignments similar, regardless of strand?
 * this ignores target gaps up to the specified size. */
bool pslTGapsSimilar(struct psl *psl1, struct psl *psl2, int maxIgnoreTGapSize);

typedef bool(*pslFilterFunc)(struct psl *psl);

/* utility to load all psl, dropping invalid ones and sort by target */
struct psl* pslEvidenceLoad(char *pslFile, pslFilterFunc filterFunc);

#endif
