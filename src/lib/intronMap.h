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
    struct intronInfo* next;  // for sorting
    char* chrom;              // location information
    int chromStart;           // zero-based, half-open
    int chromEnd;
    char transStrand[3];      // transcript strand
    char transDonor[3];       // obtain from transcript and genome
    char transAcceptor[3];
    struct starSpliceJunction* starMappings;  // list of start mappings
    struct starSpliceJunction* mappingsSum;   // sum of mappings
    struct intronTransLink* intronTranses;     // links to transcripts
};

/* map to collect splice junctions */
struct intronMap {
    struct hash* intronHash;         // intronInfo, indexed by chr:start-end
    struct genePred *transcripts;  // all transcripts
};

/* get the intron motif, either from the transcript or
 * the STAR record.  WARNING: static return  */
char* intronInfoMotifStr(struct intronInfo* intronInfo);

/* intron info is novel */
bool intronInfoIsNovel(struct intronInfo* intronInfo);

/* intron info is annotated */
static inline bool intronInfoIsAnnotated(struct intronInfo* intronInfo) {
    return !intronInfoIsNovel(intronInfo);
}

/* construct a new object */
struct intronMap* intronMapNew(void);

/* free object */
void intronMapFree(struct intronMap** intronMapPtr);

/* load a star junction file */
void intronMapLoadStarJuncs(struct intronMap* intronMap,
                            char* starJuncFile,
                            int minOverhang);

/* load a transcript genepred file */
void intronMapLoadTranscripts(struct intronMap* intronMap,
                              char* transcriptGenePred);

/* get list of intronInfo objects (DON'T FREE) */
struct intronInfo* intronMapGet(struct intronMap* intronMap);

/* get location-sorted list of intronInfo objects (DON'T FREE) */
struct intronInfo* intronMapGetSorted(struct intronMap* intronMap);

/* save splice sites obtained from transcripts to a TSV */
void intronMapSaveTranscriptSpliceSites(struct intronMap* intronMap,
                                        char* spliceTsv);

/* load splice sites obtained from transcripts to a TSV */
void intronMapLoadTranscriptSpliceSites(struct intronMap* intronMap,
                                        char* spliceTsv);

#endif
