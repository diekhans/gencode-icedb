#include "common.h"
#include "intronMap.h"
#include "starSpliceJunction.h"
#include "hash.h"
#include "genePred.h"

/* constructor */
static struct intronTransLink* intronTransLinkNew(struct genePred* transcript,
                                                  int intronIdx) {
    struct intronTransLink* intronTrans; 
    AllocVar(intronTrans);
    intronTrans->transcript = transcript;
    intronTrans->intronIdx = intronIdx;
    return intronTrans;
}

/* destructor */
static void intronTransLinkFree(struct intronTransLink** intronTransPtr) {
    freez(intronTransPtr);
}

/* list destructor */
static void intronTransLinkFreeList(struct intronTransLink** intronTransList) {
    struct intronTransLink* intronTrans;
    while ((intronTrans = slPopHead(intronTransList)) != NULL) {
        intronTransLinkFree(&intronTrans);
    }
    intronTransList = NULL;
}

/* constructor */
static struct intronInfo* intronInfoNew(void) {
    struct intronInfo* intronInfo;
    AllocVar(intronInfo);
    return intronInfo;
}

/* destructor */
static void intronInfoFree(struct intronInfo** intronInfoPtr) {
    intronTransLinkFreeList(&((*intronInfoPtr)->intronTranses));
    starSpliceJunctionFreeList(&(*intronInfoPtr)->starMappings);
    starSpliceJunctionFree(&((*intronInfoPtr)->mappingsSum));
    freez(intronInfoPtr);
}

/* report splice junction summary check difference */
static void intronInfoSumCheckErr(struct starSpliceJunction* sum,
                                  char* fieldName,
                                  unsigned sumField,
                                  unsigned starField) {
    errAbort("star splice junction %s difference for %s %d-%d, %d vs %d",
             fieldName, sum->chrom, sum->chromStart, sum->chromEnd,
             sumField, starField);
}

/* check start splice junction info compatibility with summary  */
static void intronInfoSumCheck(struct starSpliceJunction* sum,
                               struct starSpliceJunction* starJunc) {
    if (sum->strand != starJunc->strand) {
        intronInfoSumCheckErr(sum, "strand", sum->strand, starJunc->strand);
    }
    if (sum->intronMotif != starJunc->intronMotif) {
        intronInfoSumCheckErr(sum, "intronMotif", sum->intronMotif, starJunc->intronMotif);
    }
    if (sum->annotated != starJunc->annotated) {
        intronInfoSumCheckErr(sum, "annotated", sum->annotated, starJunc->annotated);
    }
}


/* summarize splice info */
static void intronInfoSum(struct intronInfo* intronInfo,
                          struct starSpliceJunction* starJunc) {
    struct starSpliceJunction* sum = intronInfo->mappingsSum;
    if (sum == NULL) {
        AllocVar(sum);
        intronInfo->mappingsSum = sum;
        *sum = *starJunc;
        sum->chrom = cloneString(starJunc->chrom);
    } else {
        intronInfoSumCheck(sum, starJunc);
        sum->numUniqueMapReads += starJunc->numUniqueMapReads;
        sum->numMultiMapReads += starJunc->numMultiMapReads;
        sum->maxOverhang = max(sum->maxOverhang, starJunc->maxOverhang);
    }
}

/* construct a new object */
struct intronMap* intronMapNew(void) {
    struct intronMap* intronMap;
    AllocVar(intronMap);
    intronMap->intronHash = hashNew(0);
    return intronMap;
}

/* free object */
void intronMapFree(struct intronMap** intronMapPtr) {
    hashFreeWithVals(&(*intronMapPtr)->intronHash, intronInfoFree);
    genePredFreeList(&(*intronMapPtr)->transcripts);
    freez(intronMapPtr);
}

/* Get a key for intron.  STATIC RETURN */
static char* getIntronKey(char* chrom, int chromStart, int chromEnd) {
    static char key[256];
    safef(key, sizeof(key), "%s:%d-%d", chrom, chromStart, chromEnd);
    return key;
}

/* get or create an introInfo object */
static struct intronInfo* intronMapObtainIntronInfo(struct intronMap* intronMap,
                                                    char* chrom,
                                                    int chromStart,
                                                    int chromEnd) {
    char* key = getIntronKey(chrom, chromStart, chromEnd);
    struct hashEl* hel = hashStore(intronMap->intronHash, key);
    if (hel->val == NULL) {
        hel->val = intronInfoNew();
    }
    return hel->val;
}

/* add a star junction object */
static void intronMapAddStarJunc(struct intronMap* intronMap,
                                 struct starSpliceJunction* starJunc) {
    struct intronInfo* intronInfo
        = intronMapObtainIntronInfo(intronMap, starJunc->chrom,
                                    starJunc->chromStart, starJunc->chromEnd);
    slAddHead(&intronInfo->starMappings, starJunc);
    intronInfoSum(intronInfo, starJunc);
}

/* load a star junction file */
void intronMapLoadStarJuncs(struct intronMap* intronMap,
                            char* starJuncFile,
                            int minOverhang) {
    struct starSpliceJunction* starJuncs = starSpliceJunctionLoadAllByTab(starJuncFile);
    struct starSpliceJunction* starJunc;
    while ((starJunc = slPopHead(&starJuncs)) != NULL) {
        if (starJuncs->maxOverhang >= minOverhang) {
            intronMapAddStarJunc(intronMap, starJunc);
        } else {
            starSpliceJunctionFree(&starJunc);
        }
    }
}

/* add a transcript intron */
static void intronMapAddTransIntron(struct intronMap* intronMap,
                                    struct genePred* transcript,
                                    int intronIdx) {
    struct intronInfo* intronInfo
        = intronMapObtainIntronInfo(intronMap, transcript->chrom,
                                    transcript->exonEnds[intronIdx],
                                    transcript->exonStarts[intronIdx]);
    struct intronTransLink* intronTrans
        = intronTransLinkNew(transcript, intronIdx);
    slAddHead(&intronInfo->intronTranses, intronTrans);
}

/* should an transcript intron be included? */
static bool shouldIncludeTransIntron(struct genePred* transcript,
                                     int intronIdx) {
    static const int minIntronSize = 30;
    return (transcript->exonStarts[intronIdx] - transcript->exonEnds[intronIdx]) >= minIntronSize;
}

/* add a transcript object */
static void intronMapAddTranscript(struct intronMap* intronMap,
                                   struct genePred* transcript) {
    for (int intronIdx = 0; intronIdx < transcript->exonCount-1; intronIdx++) {
        if (shouldIncludeTransIntron(transcript, intronIdx)) {
            intronMapAddTransIntron(intronMap, transcript, intronIdx);
        }
    }
}

/* load a transcript file */
void intronMapLoadTranscripts(struct intronMap* intronMap,
                              char* transcriptFile) {
    // FIXME: maybe drop single exon or save in another list??
    intronMap->transcripts = genePredLoadAllByTab(transcriptFile);
    for (struct genePred* transcript = intronMap->transcripts; transcript != NULL; transcript = transcript->next) {
        intronMapAddTranscript(intronMap, transcript);
    }
}
