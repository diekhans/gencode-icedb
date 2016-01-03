#include "common.h"
#include "intronCounts.h"
#include "intronMap.h"

/* constructor */
static struct intronCounts* intronCountsNew(bool annotated,
                                            int intronMotif) {
    struct intronCounts* intronCounts;
    AllocVar(intronCounts);
    intronCounts->annotated = annotated;
    intronCounts->intronMotif = intronMotif;
    return intronCounts;
}

/* get key to hash data (WARNING: static return */
static char* intronCountsKey(bool annotated,
                             int intronMotif) {
    static char key[32];
    safef(key, sizeof(key), "%d+%d", annotated, intronMotif);
    return key;
}

/* get an existing intron counts object or create a new one */
static struct intronCounts* intronCountsObtain(struct hash* intronCountsMap,
                                               bool annotated,
                                               int intronMotif) {
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
                            struct intronInfo* intronInfo) {
    struct intronCounts* intronCounts = intronCountsObtain(intronCountsMap,
                                                           annotated,
                                                           int intronMotif) {
}

/* collect intron counts from data */
struct intronCounts* intronCountsCollect(struct intronMap* intronMap) {
    struct hash* intronCountsMap = hashNew(0);
    
    hashFree(&intronCountsMap);
}
