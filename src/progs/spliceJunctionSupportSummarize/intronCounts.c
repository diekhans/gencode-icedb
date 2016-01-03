#include "common.h"
#include "intronCounts.h"
#include "intronMap.h"
#include "starSpliceJunction.h"
#include "hash.h"

/* constructor */
static struct intronCounts* intronCountsNew(bool annotated,
                                            char* intronMotif) {
    struct intronCounts* intronCounts;
    AllocVar(intronCounts);
    intronCounts->annotated = annotated;
    safecpy(intronCounts->intronMotif, sizeof(intronCounts->intronMotif), intronMotif);
    return intronCounts;
}

/* get key to hash data (WARNING: static return */
static char* intronCountsKey(bool annotated,
                             char* intronMotif) {
    static char key[32];
    safef(key, sizeof(key), "%d+%s", annotated, intronMotif);
    return key;
}

/* get an existing intron counts object or create a new one */
static struct intronCounts* intronCountsObtain(struct hash* intronCountsMap,
                                               bool annotated,
                                               char* intronMotif) {
    char* key = intronCountsKey(annotated, intronMotif);
    struct hashEl* hel = hashStore(intronCountsMap, key);
    struct intronCounts* intronCounts = hel->val;
    if (intronCounts == NULL) {
        intronCounts = hel->val = intronCountsNew(annotated, intronMotif);
    }
    return intronCounts;
}

/* count an intron map entry */
static void countIntronInfo(struct hash* intronCountsMap,
                            struct intronInfo* intronInfo,
                            struct intronCounts** intronCountsList) {
    struct intronCounts* intronCounts
        = intronCountsObtain(intronCountsMap,
                             intronInfoIsAnnotated(intronInfo),
                             intronInfoMotifStr(intronInfo));
    if (intronCounts->count == 0) {
        // first time
        slAddHead(intronCountsList, intronCounts);
    }
    intronCounts->count++;
    if (intronInfo->mappingsSum != NULL) {
        intronCounts->numUniqueMapReads += intronInfo->mappingsSum->numUniqueMapReads;
        intronCounts->numMultiMapReads += intronInfo->mappingsSum->numMultiMapReads;
    }
    intronCounts->transcriptCount += slCount(intronInfo->intronTranses);
}

/* collect intron counts from data */
struct intronCounts* intronCountsCollect(struct intronMap* intronMap) {
    // process sorted
    struct intronCounts* intronCountsList = NULL;
    struct hash* intronCountsMap = hashNew(0);
    for (struct intronInfo* intronInfo = intronMapGetSorted(intronMap); intronInfo != NULL; intronInfo = intronInfo->next) {
        countIntronInfo(intronCountsMap, intronInfo, &intronCountsList);
    }
    hashFree(&intronCountsMap);
    slReverse(&intronCountsList);
    return intronCountsList;
}
