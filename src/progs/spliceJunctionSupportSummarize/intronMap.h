#ifndef intronMap_h
#define intronMap_h

/* linking of an intron to a transcript */
struct intronTransLink {
    struct intronTransLink* next;
    struct genePred* transcript;  // pointer to transcript (not owned)
    int intronIdx;                // index of intron (follows this exn)
};

/* information about an intron */
struct intronInfo {
    struct starSpliceJunction* starMappings;  // list of start mappings
    struct starSpliceJunction* mappingsSum;   // sum of mappings
    struct intronTransLink* intronTranses;     // links to transcripts
};

/* map to collect splice junctions */
struct intronMap {
    struct hash* intronHash;         // intronInfo, indexed by chr:start-end
    struct genePred *transcripts;  // all transcripts
};

/* construct a new object */
struct intronMap* intronMapNew(void);

/* free object */
void intronMapFree(struct intronMap** intronMapPtr);

/* load a star junction file */
void intronMapLoadStarJuncs(struct intronMap* intronMap,
                            char* starJuncFile,
                            int minOverhang);

/* load a transcript file */
void intronMapLoadTranscripts(struct intronMap* intronMap,
                              char* transcriptFile);

#endif
