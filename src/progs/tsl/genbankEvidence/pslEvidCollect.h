#ifndef pslEvidCollect_h
#define pslEvidCollect_h
struct psl;
#include <stdbool.h>

/* Object to collect evidence from a PSLs and write to external evidence file
 * for loading into databases.  By collecting the splice sites up front, this
 * allows detecting and adjusted reverse-complemented sequences.
 */
struct pslEvidCollect;

/* constructor */
struct pslEvidCollect* pslEvidCollectNew(char *twoBitFile, char *cdnaAlignFile, bool ignoreMatch);

/* destructor */
void pslEvidCollectFree(struct pslEvidCollect* self);

/* Analyze a psl and record splice information for possible output.
 * WARNING: psl maybe reverse-complemented. */
void pslEvidCollectAnalyze(struct pslEvidCollect* self, struct psl *psl, int representsCnt);

/* Write evidence for current psl */
void pslEvidCollectWrite(struct pslEvidCollect* self);

/* determine weighted direction of based on apparent intron  */
int pslEvidCollectWeightedDirection(struct pslEvidCollect* self);

/* reverse complement current psl and splice sites */
void pslEvidCollectReverseComplement(struct pslEvidCollect* self);

#endif
