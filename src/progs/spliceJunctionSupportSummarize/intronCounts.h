#ifndef intronCounts_h
#define intronCounts_h
struct intronMap;

/* count of support for an intron of a certain category */
struct intronCounts {
    struct intronCounts* next;
    bool annotated;    // is it annotated?
    char intronMotif[6];
    int count;
    int numUniqueMapReads;
    int numMultiMapReads;
    int transcriptCount;
};

/* collect intron counts from data */
struct intronCounts* intronCountsCollect(struct intronMap* intronMap);

#endif
