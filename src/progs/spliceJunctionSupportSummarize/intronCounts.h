#ifndef supportCounts_h
#define supportCounts_h

/* count of support for an intron of a certain category */
struct intronCounts {
    struct intronCounts* next;
    bool annotated;    // is it annotated?
    int intronMotif;   // motif code, -1 if no support
    int intronCount;
    int numUniqueMapReads;
    int numMultiMapReads;
    int transcriptCount;
};

/* collect intron counts from data */
struct intronCounts* intronCountsCollect(struct intronMap* intronMap);

#endif
