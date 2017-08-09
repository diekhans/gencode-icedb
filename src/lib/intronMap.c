#include "common.h"
#include "intronMap.h"
#include "starSpliceJunction.h"
#include "starOps.h"
#include "hash.h"
#include "sqlNum.h"
#include "genePred.h"
#include "linefile.h"
#include "rslAnalysisSet.h"

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
static struct intronInfo* intronInfoNew(char* chrom, int chromStart, int chromEnd) {
    struct intronInfo* intronInfo;
    AllocVar(intronInfo);
    intronInfo->chrom = cloneString(chrom);
    intronInfo->chromStart = chromStart;
    intronInfo->chromEnd = chromEnd;
    return intronInfo;
}

/* destructor */
static void intronInfoFree(struct intronInfo** intronInfoPtr) {
    intronTransLinkFreeList(&((*intronInfoPtr)->intronTranses));
    starSpliceJunctionFreeList(&(*intronInfoPtr)->starMappings);
    starSpliceJunctionFree(&((*intronInfoPtr)->mappingsSum));
    freez(&(*intronInfoPtr)->chrom);
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

/* collect summarize splice info for the first time */
static void intronInfoSumFirst(struct intronInfo* intronInfo,
                               struct starSpliceJunction* starJunc) {
    struct starSpliceJunction* sum;
    AllocVar(sum);
    intronInfo->mappingsSum = sum;
    *sum = *starJunc;  // byte copy
    // fix up dynamic pointers, srcAnalyses should only be one
    sum->next = NULL;
    sum->chrom = cloneString(starJunc->chrom);
    sum->srcAnalyses = rslAnalysisLinkCloneList(starJunc->srcAnalyses);  // normally one
}

/* collect summarize splice info for the subsequence times */
static void intronInfoSumRest(struct intronInfo* intronInfo,
                              struct starSpliceJunction* starJunc) {
    struct starSpliceJunction* sum = intronInfo->mappingsSum;
    intronInfoSumCheck(sum, starJunc);
    sum->numUniqueMapReads += starJunc->numUniqueMapReads;
    sum->numMultiMapReads += starJunc->numMultiMapReads;
    sum->maxOverhang = max(sum->maxOverhang, starJunc->maxOverhang);
    sum->srcAnalyses = slCat(sum->srcAnalyses,
                             rslAnalysisLinkCloneList(starJunc->srcAnalyses));
}

/* summarize splice info */
static void intronInfoSum(struct intronInfo* intronInfo,
                          struct starSpliceJunction* starJunc) {
    if (intronInfo->mappingsSum == NULL) {
        intronInfoSumFirst(intronInfo, starJunc);
    } else {
        intronInfoSumRest(intronInfo, starJunc);
    }
}

/* get the intron motif, either from the transcript or
 * the STAR record.  WARNING: static return */
char* intronInfoMotifStr(struct intronInfo* intronInfo) {
    if (strlen(intronInfo->transDonor) > 0) {
        static char motif[6];
        safef(motif, sizeof(motif), "%s/%s", intronInfo->transDonor, intronInfo->transAcceptor);
        return motif;
    } else if (intronInfo->mappingsSum != NULL) {
        return starMotifCodeToStr(intronInfo->mappingsSum->intronMotif);
    } else {
        return "\?\?/\?\?"; 
    }
}

/* intron info is novel */
bool intronInfoIsNovel(struct intronInfo* intronInfo) {
    if (intronInfo->mappingsSum != NULL) {
        if ((intronInfo->mappingsSum->annotated == 0) != (intronInfo->intronTranses == NULL)) {
            errAbort("intron %s:%d-%d: STAR annotated state (%d) not the same as transcript state (%d)",
                     intronInfo->chrom, intronInfo->chromStart, intronInfo->chromEnd,
                     (intronInfo->mappingsSum->annotated == 0),
                     (intronInfo->intronTranses == NULL));
        }
        return (intronInfo->mappingsSum->annotated == 0);
    } else {
        return (intronInfo->intronTranses == NULL);
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
        hel->val = intronInfoNew(chrom, chromStart, chromEnd);
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

/* read splice junctions from file and link if source analysis information */
static struct starSpliceJunction* spliceJunctionsLoad(struct rslAnalysis* rslAnalysis) {
    struct starSpliceJunction* starJuncs = starSpliceJunctionLoadAllByTab(rslAnalysis->sjPath);
    for (struct starSpliceJunction* starJunc = starJuncs; starJunc != NULL; starJunc = starJunc->next) {
        slAddHead(&starJunc->srcAnalyses, rslAnalysisLinkNew(rslAnalysis));
    }
    return starJuncs;
}

/* load a star junction file */
void intronMapLoadStarJuncs(struct intronMap* intronMap,
                            struct rslAnalysis* rslAnalysis,
                            int minOverhang) {
    struct starSpliceJunction* starJuncs = spliceJunctionsLoad(rslAnalysis);
    struct starSpliceJunction* starJunc;
    while ((starJunc = slPopHead(&starJuncs)) != NULL) {
        if (starJunc->maxOverhang >= minOverhang) {
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
                                    transcript->exonStarts[intronIdx+1]);
    struct intronTransLink* intronTrans
        = intronTransLinkNew(transcript, intronIdx);
    slAddHead(&intronInfo->intronTranses, intronTrans);
    safecpy(intronInfo->transStrand, sizeof(intronInfo->transStrand), transcript->strand);
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
                              char* transcriptGenePred) {
    intronMap->transcripts = genePredLoadAllByTab(transcriptGenePred);
    for (struct genePred* transcript = intronMap->transcripts; transcript != NULL; transcript = transcript->next) {
        intronMapAddTranscript(intronMap, transcript);
    }
}

/* compare by location */
static int intronInfoLocCmp(const void *va, const void *vb) {
    const struct intronInfo *a = *((struct intronInfo **)va);
    const struct intronInfo *b = *((struct intronInfo **)vb);
    int diff = strcmp(a->chrom, b->chrom);
    if (diff != 0) {
        return diff;
    }
    diff = a->chromStart - b->chromStart;
    if (diff != 0) {
        return diff;
    }
    return a->chromEnd - b->chromEnd;
}

/* get list of intronInfo objects (DON'T FREE) */
struct intronInfo* intronMapGet(struct intronMap* intronMap) {
    struct intronInfo* intronInfos = NULL;
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(intronMap->intronHash);
    while ((hel = hashNext(&cookie)) != NULL) {
        struct intronInfo* intronInfo = hel->val;
        slAddHead(&intronInfos, intronInfo);
    }
    return intronInfos;
}

/* get location-sorted list of intronInfo objects (DON'T FREE) */
struct intronInfo* intronMapGetSorted(struct intronMap* intronMap) {
    struct intronInfo* intronInfos = intronMapGet(intronMap);
    slSort(&intronInfos, intronInfoLocCmp);
    return intronInfos;
}

/* header for splice sites */
static const int spliceTsvWidth = 6;
static char* spliceTsvHeader[] = {
    "chrom", "chromStart", "chromEnd",
    "strand", "donor", "acceptor", "transcripts", NULL
};

/* write intron splice site TSV header */
static void intronMapWriteSpliceTsvHeader(FILE* spliceTsvFh) {
    for (int i = 0; spliceTsvHeader[i] != NULL; i++) {
        if (i > 0) {
            fputc('\t', spliceTsvFh);
        }
        fputs(spliceTsvHeader[i], spliceTsvFh);
    }
    fputc('\n', spliceTsvFh);
}

/* write a row of intron splice site TSV */
void intronMapWriteSpliceTsvRow(FILE* spliceTsvFh,
                                struct intronInfo* intronInfo) {
    fprintf(spliceTsvFh, "%s\t%d\t%d\t%s\t%s\t%s\t",
            intronInfo->chrom, 
            intronInfo->chromStart, intronInfo->chromEnd,
            intronInfo->transStrand,
            intronInfo->transDonor,
            intronInfo->transAcceptor);
    for (struct intronTransLink *transLink = intronInfo->intronTranses; transLink != NULL; transLink = transLink->next) {
        if (transLink != intronInfo->intronTranses) {
            fputc(',', spliceTsvFh);
        }
        fputs(transLink->transcript->name, spliceTsvFh);
    }
    fputc('\n', spliceTsvFh);
}

/* save splice sites obtained from transcripts to a TSV */
void intronMapSaveTranscriptSpliceSites(struct intronMap* intronMap,
                                        char* spliceTsv) {
    FILE* spliceTsvFh = mustOpen(spliceTsv, "w");
    intronMapWriteSpliceTsvHeader(spliceTsvFh);
    struct intronInfo* intronInfos = intronMapGetSorted(intronMap);
    for (struct intronInfo* intronInfo = intronInfos; intronInfo != NULL; intronInfo = intronInfo->next) {
        if (strlen(intronInfo->transDonor) > 0) {
            intronMapWriteSpliceTsvRow(spliceTsvFh, intronInfo);
        }
    }
    carefulClose(&spliceTsvFh);
}

/* check header */
static void intronMapLoadCheckHeader(struct lineFile *spliceTsvLf) {
    char* row[spliceTsvWidth];
    if (!lineFileNextRowTab(spliceTsvLf, row, spliceTsvWidth)) {
        errAbort("premature EOF on splice TSV"); 
    }
    for (int i = 0; spliceTsvHeader[i] != NULL; i++) {
        if (!sameString(row[i], spliceTsvHeader[i])) {
            errAbort("unexpected splice TSV column header \"%s\", expected \"%s\"",
                     row[i], spliceTsvHeader[i]);
        }
    }
}

/* load one transcript splice site */
static void intronMapLoadTranscriptSpliceSite(struct intronMap* intronMap,
                                              char* chrom, int chromStart, int chromEnd,
                                              char* strand, char* donor, char* acceptor) {
    struct intronInfo* intronInfo =  intronMapObtainIntronInfo(intronMap, chrom, chromStart, chromEnd);
    safecpy(intronInfo->transStrand, sizeof(intronInfo->transStrand), strand);
    safecpy(intronInfo->transDonor, sizeof(intronInfo->transDonor), donor);
    safecpy(intronInfo->transAcceptor, sizeof(intronInfo->transAcceptor), acceptor);
}

/* load splice sites obtained from transcripts to a TSV */
void intronMapLoadTranscriptSpliceSites(struct intronMap* intronMap,
                                        char* spliceTsv) {
    char* row[spliceTsvWidth];
    struct lineFile *spliceTsvLf = lineFileOpen(spliceTsv, TRUE);
    intronMapLoadCheckHeader(spliceTsvLf) ;

    while (lineFileNextRowTab(spliceTsvLf, row, spliceTsvWidth)) {
        intronMapLoadTranscriptSpliceSite(intronMap,
                                          row[0], 
                                          sqlSigned(row[1]),
                                          sqlSigned(row[2]),
                                          row[3], row[4], row[5]);
    }
    lineFileClose(&spliceTsvLf);
}
