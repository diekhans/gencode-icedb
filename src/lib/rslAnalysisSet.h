#ifndef rslAnalysisSets_h
#define rslAnalysisSets_h

/* An RSL analysis for a geneset using one star run */
struct rslAnalysis {
    struct rslAnalysis *next;
    char *runname;
    char *tissue;
    char *sjPath;
};

/* linked list of rslAnalysis object.  These don't own the objects */
struct rslAnalysisLink {
    struct rslAnalysisLink *next;
    struct rslAnalysis *rslAnalysis;
};

/* allocate a link */
INLINE struct rslAnalysisLink* rslAnalysisLinkNew(struct rslAnalysis *rslAnalysis) {
    struct rslAnalysisLink *ral;
    AllocVar(ral);
    ral->rslAnalysis = rslAnalysis;
    return ral;
}

/* clone a list of links (but not objects) */
INLINE struct rslAnalysisLink *rslAnalysisLinkCloneList(struct rslAnalysisLink *head) {
    struct rslAnalysisLink *head2 = NULL;;
    for (struct rslAnalysisLink *ral = head; ral != NULL; ral = ral->next) {
        slAddHead(&head2, rslAnalysisLinkNew(ral->rslAnalysis));
    }
    slReverse(&head2);
    return head2;
}

/* free list of links (but not objects) */
INLINE void rslAnalysisLinkFreeList(struct rslAnalysisLink **head) {
    struct rslAnalysisLink *ral;
    while ((ral = slPopHead(head)) != NULL) {
        freeMem(ral);
    }
}
    
/* a set of RSL analysis for a geneset */
struct rslAnalysisSet {
    char *setname;
    struct rslAnalysis *analyses;
};

/* load an RSL analysis set from a sqlite3 database plus output files */
struct rslAnalysisSet *rslAnalysisSetLoad(char *tsvFile,
                                          char *setname);

/* free an RSL analysis geneset */
void rslAnalysisSetFree(struct rslAnalysisSet *rslAnalysisSet);


#endif
