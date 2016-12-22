#include "common.h"
#include "options.h"
#include "psl.h"
#include "verbose.h"
#include "pslEvidCollect.h"
#include "pslEvidence.h"

static struct optionSpec optionSpecs[] = {
    {"ignoreMatch", OPTION_BOOLEAN},
    {"allMultiExon", OPTION_BOOLEAN},
    {"inclSingleExon", OPTION_BOOLEAN},
    {"maxIgnoreTGapSize", OPTION_INT},
    {NULL, 0}
};
static bool ignoreMatch = FALSE;
static bool allMultiExon = FALSE;
static bool inclSingleExon = FALSE;
static int maxIgnoreTGapSize = 10;

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "estPslEvidence estPslFile twoBitFile cdnaAlignFile\n"
        "\n"
        "Pick the best of a set of consistent spliced ESTs, discarding others.\n"
        "Determine the direction of transcription and reverse-complement.\n"
        "ESTs without any minimum intron sized target gaps are discarded unless\n"
        "single-exon are requested.\n"
        "\n"
        "cdnaAlignFile are a tab-separated file describe in pslEvidCollect.c.\n"
        "\n"
        "Options:\n"
        "  -ignoreMatch - match/mismatch counts are incorrect and should be ignored\n"
        "  -inclSingleExon - include single-exon transcripts\n"
        "  -maxIgnoreTGapSize=10 - maximum size of a gap to ignore when picking ESTs that have\n"
        "   the same structure.\n"
        "  -allMultiExon - don't output representative ESTs, output all apparent multi-exon ESTs.\n"
        "  -verbose=n\n";
    errAbort("%s:\n%s", msg, usageMsg);
}

/* does a EST appear to have a minimum sized intron? */
static bool hasMinIntron(struct psl *psl) {
    if (inclSingleExon) {
        return true;
    }
    for (int iBlk = 1; iBlk < psl->blockCount; iBlk++) {
        if (pslTGapSize(psl, iBlk) >= pslEvidenceMinIntronSize) {
            return true;
        }
    }
    return false;
}

/* pop a set of EST with overlapping blocks, discarding ESTs with no target
 * gaps of sufficient size. Assumes list is sort by target lowest start and
 * highest end.  so longest overlapping is found first. completely contained
 * psls are returned to the list */
static struct psl *pslPopWithOverlapingTBlocks(struct psl **psls) {
    struct psl *longestPsl = NULL, *overPsls = NULL, *savedPsls = NULL, *psl;
    while ((psl = slPopHead(psls)) != NULL) {
        if (longestPsl == NULL) {
            // first psl to keep
            longestPsl = overPsls = psl;
        } else if (!pslTOverlap(psl, longestPsl)) {
            slAddHead(&savedPsls, psl);   // follows longest
            break;  // LEAVE LOOP!!
        } else if (pslTGapsSimilar(psl, longestPsl, maxIgnoreTGapSize)) {
            slAddHead(&overPsls, psl);
        } else {
            slAddHead(&savedPsls, psl);  // in gap, save for next pass
        }
    }

    // put back for later
    slReverse(&savedPsls);
    *psls = slCat(savedPsls, *psls);

    slReverse(&overPsls);
    return overPsls;
}

/* check splice sites and reverse complement as needed. */
static void correctForTranscriptionDir(struct pslEvidCollect* pslEvidCollect) {
    // FIXME: should indicate when direction can't be determined.
    // Only reverse if evidence shows sequence is reverse-complemented.
    // Note: that pslRc has already been called to put in direction of
    // transcription.
    int dir = pslEvidCollectWeightedDirection(pslEvidCollect);
    if (dir < 0) {
        pslEvidCollectReverseComplement(pslEvidCollect);
    }
}

static void verbosePrintOverset(struct psl *pslOverSet) {
    struct psl *psl0 = pslOverSet;
    fprintf(stderr, "overlapping set: %s %s:%d-%d\n", psl0->qName, psl0->tName, psl0->tStart, psl0->tEnd);
    for (struct psl *psl = psl0->next; psl != NULL; psl = psl->next) {
        fprintf(stderr, "\t%s %s:%d-%d\n", psl->qName, psl->tName, psl->tStart, psl->tEnd);
    }
}

/* select and write from set of overlapping ESTs.  FIXME: just write longest
 * one for starts, which is stupid.  */
static void processEstOverlapSet(struct pslEvidCollect* pslEvidCollect, struct psl *pslOverSet) {
    struct psl *psl0 = pslOverSet;
    pslEvidCollectAnalyze(pslEvidCollect, psl0, slCount(pslOverSet));
    correctForTranscriptionDir(pslEvidCollect);
    pslEvidCollectWrite(pslEvidCollect);
    if (verboseLevel() > 1) {
        verbosePrintOverset(pslOverSet);
    }
}

/* extract the evidence and summarize as a representative */
static void representativeEvidence(struct pslEvidCollect* pslEvidCollect, struct psl *psls) {
    while (psls != NULL) {
        struct psl *pslOverSet = pslPopWithOverlapingTBlocks(&psls);
        processEstOverlapSet(pslEvidCollect, pslOverSet);
        pslFreeList(&pslOverSet);
    }
}

/* extract the evidence without summarizing */
static void allEvidence(struct pslEvidCollect* pslEvidCollect, struct psl *psls) {
    struct psl *psl;
    while ((psl = slPopHead(&psls)) != NULL) {
        processEstOverlapSet(pslEvidCollect, psl);
        pslEvidCollectAnalyze(pslEvidCollect, psl, 1);
        correctForTranscriptionDir(pslEvidCollect);
        pslEvidCollectWrite(pslEvidCollect);
        pslFreeList(&psl);
    }
}

/* extract the evidence */
static void estPslEvidence(char *cdnaPslFile, char *twoBitFile, char *cdnaAlignFile) {
    struct pslEvidCollect* pslEvidCollect = pslEvidCollectNew(twoBitFile, cdnaAlignFile, ignoreMatch);
    struct psl *psls = pslEvidenceLoad(cdnaPslFile, hasMinIntron);
    if (allMultiExon) {
        allEvidence(pslEvidCollect, psls);
    } else {
        representativeEvidence(pslEvidCollect, psls);
    }
    pslEvidCollectFree(pslEvidCollect);
}

/* entry */
int main(int argc, char** argv) {
    dnaUtilOpen();
    optionInit(&argc, argv, optionSpecs);
    if (argc != 4)
        usage("wrong # args");
    ignoreMatch = optionExists("ignoreMatch");
    allMultiExon = optionExists("allMultiExon");
    inclSingleExon = optionExists("inclSingleExon");
    maxIgnoreTGapSize = optionInt("maxIgnoreTGapSize", maxIgnoreTGapSize);
    estPslEvidence(argv[1], argv[2], argv[3]);
    return 0;
}
