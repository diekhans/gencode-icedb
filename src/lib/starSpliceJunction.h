/* starSpliceJunction.h was originally generated by the autoSql program, which also 
 * generated starSpliceJunction.c and starSpliceJunction.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef STARSPLICEJUNCTION_H
#define STARSPLICEJUNCTION_H

#define STARSPLICEJUNCTION_NUM_COLS 9

extern char *starSpliceJunctionCommaSepFieldNames;

struct starSpliceJunction
/* Splice junction output from STAR */
    {
    struct starSpliceJunction *next;  /* Next in singly linked list. */
    char *chrom;	/* chromsome */
    unsigned chromStart;	/* first base of the intron (1-based in file, converted in memory) */
    unsigned chromEnd;	/* last base of the intron */
    unsigned strand;	/* strand 0: undefined, 1: +, 2: - */
    unsigned intronMotif;	/* 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT */
    unsigned annotated;	/* 0: unannotated, 1: annotated (only if splice junctions database is used) */
    unsigned numUniqueMapReads;	/* number of uniquely mapping reads crossing the junction */
    unsigned numMultiMapReads;	/* number of multi-mapping reads crossing the junction */
    unsigned maxOverhang;	/* maximum spliced alignment overhang */
    /* extra data */ 
    struct rslAnalysisLink* srcAnalyses;  /* links to analysis that were the source of this object.
                                           * it maybe a list when summed.  This is not in the file,
                                           * will be freed when object is freed. */
    };

struct starSpliceJunction *starSpliceJunctionLoad(char **row);
/* Load a starSpliceJunction from row fetched with select * from starSpliceJunction
 * from database.  Dispose of this with starSpliceJunctionFree(). */

struct starSpliceJunction *starSpliceJunctionLoadAll(char *fileName);
/* Load all starSpliceJunction from whitespace-separated file.
 * Dispose of this with starSpliceJunctionFreeList(). */

struct starSpliceJunction *starSpliceJunctionLoadAllByChar(char *fileName, char chopper);
/* Load all starSpliceJunction from chopper separated file.
 * Dispose of this with starSpliceJunctionFreeList(). */

#define starSpliceJunctionLoadAllByTab(a) starSpliceJunctionLoadAllByChar(a, '\t');
/* Load all starSpliceJunction from tab separated file.
 * Dispose of this with starSpliceJunctionFreeList(). */

struct starSpliceJunction *starSpliceJunctionCommaIn(char **pS, struct starSpliceJunction *ret);
/* Create a starSpliceJunction out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new starSpliceJunction */

void starSpliceJunctionFree(struct starSpliceJunction **pEl);
/* Free a single dynamically allocated starSpliceJunction such as created
 * with starSpliceJunctionLoad(). */

void starSpliceJunctionFreeList(struct starSpliceJunction **pList);
/* Free a list of dynamically allocated starSpliceJunction's */

void starSpliceJunctionOutput(struct starSpliceJunction *el, FILE *f, char sep, char lastSep);
/* Print out starSpliceJunction.  Separate fields with sep. Follow last field with lastSep. */

#define starSpliceJunctionTabOut(el,f) starSpliceJunctionOutput(el,f,'\t','\n');
/* Print out starSpliceJunction as a line in a tab-separated file. */

#define starSpliceJunctionCommaOut(el,f) starSpliceJunctionOutput(el,f,',',',');
/* Print out starSpliceJunction as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */
#endif /* STARSPLICEJUNCTION_H */

