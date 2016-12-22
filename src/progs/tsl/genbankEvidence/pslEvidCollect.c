#include "common.h"
#include "pslEvidCollect.h"
#include "dnautil.h"
#include "twoBit.h"
#include "psl.h"
#include "binRange.h"
#include "pslEvidence.h"
#include "dnautil.h"
#include <stdbool.h>

/*
 * Each output line consists an ALN record at the start of each alignment
 * followed by BLK records
 *
 * ALN cdnaAccver cdnaSize cdnaStrand cdnaStart cdnaEnd chrom chromSize chromStrand chromStart chromEnd chromBin ident representsCnt
 *   - ident can be empty if not computed.
 *   - representsCnt is used for ESTs when multiple ESTs are represented by a single best.
 *
 * BLK cdnaRelStart cdnaRelEnd chromRelStart chromRelEnd donor acceptor
 *   - BLK records are emitted for block or gaps, with the one of the start
 *     records empty for a gap.  The donor and acceptor are only written on
 *     target inserts of at least 4 bases.
 */

static const bool debug = false; 

static const int spliceSitesSizeIncr = 32;  // initial size and increment size

/* object that holds donor and acceptor bases */
struct spliceSite {
    char donor[3];
    char acceptor[3];
};

static struct spliceSite spliceSiteInit = {
    {'\0', '\0', '\0'}, 
    {'\0', '\0', '\0'}, 
};


/* Object to collect evidence from a PSLs and write to external evidence file
 * for loading into databases.  By collecting the splice sites up front, this
 * allows detecting and adjusted reverse-complemented sequences.
 */
struct pslEvidCollect {
    struct twoBitFile *genomeSeqs;
    FILE *outFh;
    bool ignoreMatch;
    int representsCnt;
    struct psl *psl; // current psl being analyzed, reversed-complement so Q
                     // is always positive unless sequence is reversed.
    struct spliceSite *spliceSites;  // indexed by following block
    int spliceSiteCapacity; // current capacity
    bool reverseComplemented;  // current has been reverse complemented
};

/* convert a splice site to a string. WARNING: static return */
static char *spliceSiteToString(struct spliceSite ss) {
    static char buf[32];
    safef(buf, sizeof(buf), "%s..%s", ss.donor, ss.acceptor);
    return buf;
}

/* reverse complement a splice site and return it */
static struct spliceSite spliceSiteReverseComplement(struct spliceSite ss) {
    struct spliceSite rss;
    strcpy(rss.donor, ss.acceptor);
    reverseComplement(rss.donor, strlen(rss.donor));
    strcpy(rss.acceptor, ss.donor);
    reverseComplement(rss.acceptor, strlen(rss.acceptor));
    return rss;
}

/* constructor */
struct pslEvidCollect* pslEvidCollectNew(char *twoBitFile, char *cdnaAlignFile, bool ignoreMatch) {
    struct pslEvidCollect* self;
    AllocVar(self);
    self->ignoreMatch = ignoreMatch;
    self->genomeSeqs = twoBitOpen(twoBitFile);
    self->outFh = mustOpen(cdnaAlignFile, "w");
    self->spliceSiteCapacity = spliceSitesSizeIncr;
    AllocArray(self->spliceSites, self->spliceSiteCapacity);
    return self;
}

/* destructor */
void pslEvidCollectFree(struct pslEvidCollect* self) {
    twoBitClose(&self->genomeSeqs);
    carefulClose(&self->outFh);
    freeMem(self);
}

/* Ensure that there is splice site array for the current psl, expanded if
 * needed. Zero all entries. */
static void spliceSitesSetup(struct pslEvidCollect* self) {
    if (self->psl->blockCount > self->spliceSiteCapacity) {
        int newSize = self->spliceSiteCapacity + self->psl->blockCount + spliceSitesSizeIncr;
        ExpandArray(self->spliceSites, self->spliceSiteCapacity, newSize);
        self->spliceSiteCapacity = newSize;
    }
    for (int iBlk = 0; iBlk < self->psl->blockCount; iBlk++) {
        self->spliceSites[iBlk] = spliceSiteInit;
    }
}

/* get the specified number of target bases, reverse complemented base on psl */
static void getTargetBases(struct twoBitFile *genomeSeqs, struct psl *psl, int tStartRel, int numBases, char *bases) {
    int tStart = tStartRel, tEnd = tStartRel + numBases;
    if (pslTStrand(psl) == '-') {
        reverseIntRange(&tStart, &tEnd, psl->tSize);
    }
    struct dnaSeq *dna = twoBitReadSeqFragExt(genomeSeqs, psl->tName, tStart, tEnd, FALSE, NULL);
    strncpy(bases, dna->dna, numBases);
    bases[numBases] = '\0';
    dnaSeqFree(&dna);
    if (pslTStrand(psl) == '-') {
        reverseComplement(bases, numBases);
    }
}

/* check for and record spliced sites BEFORE the specified block, if there is
 * at least a four bases target insert in the gap. */
static void saveGapSplice(struct pslEvidCollect* self, int iBlk) {
    if (pslTGapSize(self->psl, iBlk) >= 4) {
        getTargetBases(self->genomeSeqs, self->psl, pslTEnd(self->psl, iBlk-1), 2, self->spliceSites[iBlk].donor);
        getTargetBases(self->genomeSeqs, self->psl, pslTStart(self->psl, iBlk)-2, 2, self->spliceSites[iBlk].acceptor);
        if (debug) {
            fprintf(stderr, "%s[%d] gap: %s\n", self->psl->tName, iBlk, spliceSiteToString(self->spliceSites[iBlk]));
        }
    }
}

/* record splice sites from one alignment.
 * WARNING: should be called after setting reverse complementing for alignment
 * direction.  Doesn't do any reverse complement based on strand. 
 */
static void recordSpliceSites(struct pslEvidCollect* self) {
    for (int iBlk = 1; iBlk < self->psl->blockCount; iBlk++) {
        saveGapSplice(self, iBlk);
    }
}

/* Analyze a psl and record splice information for possible output.
 * WARNING: psl maybe reverse-complemented. */
void pslEvidCollectAnalyze(struct pslEvidCollect* self, struct psl *psl, int representsCnt) {
    self->psl = psl;
    self->reverseComplemented = false;
    self->representsCnt = representsCnt;
    // force blocks to transcription order *MUST* do first, before splice site collection
    if (pslQStrand(self->psl) == '-') {
        pslRc(self->psl);
    }
    spliceSitesSetup(self);
    recordSpliceSites(self);
}

/* write an ALN record */
static void alnWrite(struct pslEvidCollect* self) {
    fprintf(self->outFh, "ALN\t%s\t%d\t%c\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%d\t",
            self->psl->qName, self->psl->qSize, pslQStrand(self->psl), self->psl->qStart, self->psl->qEnd,
            self->psl->tName, self->psl->tSize, pslTStrand(self->psl), self->psl->tStart, self->psl->tEnd,
            binFromRange(self->psl->tStart, self->psl->tEnd));
    if (!self->ignoreMatch) {
        fprintf(self->outFh, "%0.4f", pslIdent(self->psl));
    }
    fprintf(self->outFh, "\t%d", self->representsCnt);
    fputc('\n', self->outFh);
}

/* write an BLK record, -1 for positions for gaps */
static void blkWrite(struct pslEvidCollect* self, int cdnaRelStart, int cdnaRelEnd, int chromRelStart, int chromRelEnd, char *donor, char *acceptor) {
    fprintf(self->outFh, "BLK");
    if (cdnaRelStart < 0) {
        fprintf(self->outFh, "\t\t");
    } else {
        fprintf(self->outFh, "\t%d\t%d", cdnaRelStart, cdnaRelEnd);
    }
    if (chromRelStart < 0) {
        fprintf(self->outFh, "\t\t");
    } else {
        fprintf(self->outFh, "\t%d\t%d", chromRelStart, chromRelEnd);
    }
    fprintf(self->outFh, "\t%s\t%s\n", donor, acceptor);
}

/* process a pairwise aligned block */
static void processBlock(struct pslEvidCollect* self, int iBlk) {
    blkWrite(self, 
             pslQStart(self->psl, iBlk), pslQEnd(self->psl, iBlk),
             pslTStart(self->psl, iBlk), pslTEnd(self->psl, iBlk),
             "", "");
}

/* process a gap BEFORE the specified block */
static void processGap(struct pslEvidCollect* self, int iBlk) {
    if (pslQEnd(self->psl, iBlk-1) < self->psl->qStarts[iBlk]) {
        // query insert
        blkWrite(self, pslQEnd(self->psl, iBlk-1), self->psl->qStarts[iBlk], -1, -1, "", "");
    }
    if (pslTEnd(self->psl, iBlk-1) < pslTStart(self->psl, iBlk)) {
        // target insert
        blkWrite(self,  -1, -1, pslTEnd(self->psl, iBlk-1), pslTStart(self->psl, iBlk), self->spliceSites[iBlk].donor, self->spliceSites[iBlk].acceptor);
    }
}

/* Write evidence for current psl */
void pslEvidCollectWrite(struct pslEvidCollect* self) {
    alnWrite(self);
    for (int iBlk = 0; iBlk < self->psl->blockCount; iBlk++) {
        if (iBlk > 0) {
            processGap(self, iBlk);
        }
        processBlock(self, iBlk);
    }
}

/* compare two splice site objects */
static bool spliceSiteEq(struct spliceSite *ss1, struct spliceSite *ss2) {
    return sameString(ss1->donor, ss2->donor) && sameString(ss1->acceptor, ss2->acceptor);
}

/* get direction of a spliceSite, Only consider gt..ag and gc..ag for now. */
static int spliceSiteDir(struct spliceSite *ss) {
    static const int numCmp = 2;
    static struct spliceSite fwdSplices[] = {{"gt", "ag"}, {"gc", "ag"}};
    static struct spliceSite revSplices[] = {{"ct", "ac"}, {"ct", "gc"}};
    for (int i = 0; i < numCmp; i++) {
        if (spliceSiteEq(&fwdSplices[i], ss)) {
            return 1;
        }
    }
    for (int i = 0; i < numCmp; i++) {
        if (spliceSiteEq(&revSplices[i], ss)) {
            return -1;
        }
    }
    return 0;
}

/* determine direction of an apparent intron  */
static int intronDir(struct pslEvidCollect* self, int iBlk) {
    if (pslTGapSize(self->psl, iBlk) >= pslEvidenceMinIntronSize) {
        return spliceSiteDir(&self->spliceSites[iBlk]);
    } else {
        return 0;
    }
}

/* determine weighted direction of based on apparent intron  */
int pslEvidCollectWeightedDirection(struct pslEvidCollect* self) {
    int dir = 0;
    // look at gaps between blocks
    for (int iBlk = 1; iBlk < self->psl->blockCount; iBlk++) {
        dir += intronDir(self, iBlk);
    }
    return dir;
}

/* reverse complement current psl and splice sites */
void pslEvidCollectReverseComplement(struct pslEvidCollect* self) {
    assert(pslQStrand(self->psl) == '+');
    self->reverseComplemented = true;
    pslRc(self->psl);
    // reverse splice site array and reverse-complement splice sites
    for (int iBlk = 0; iBlk < self->psl->blockCount; iBlk++) {
        self->spliceSites[iBlk] = spliceSiteReverseComplement(self->spliceSites[iBlk]);
        if (debug) {
            fprintf(stderr, "%s[%d] rc-gap: %s\n", self->psl->tName, iBlk, spliceSiteToString(self->spliceSites[iBlk]));
        }
    }
}

