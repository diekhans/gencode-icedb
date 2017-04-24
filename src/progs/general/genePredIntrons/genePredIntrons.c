#include "common.h"
#include "intronMap.h"
#include "options.h"
#include "twoBit.h"
#include "dnaseq.h"
#include "dnautil.h"

static struct optionSpec optionSpecs[] = {
    {NULL, 0}
};

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "genePredIntrons genePred twoBit spliceTsv\n\n"
        "Collect introns splice information for a geneSet\n"
        "\n"
        "Options:\n"
        ;
    errAbort("%s:\n%s", msg, usageMsg);
}

/** donor/acceptor pair */
struct splicing {
    char donor[3];
    char accept[3];
};

/* construct a splicing object */
static struct splicing splicingMk(char *donor, char *accept) {
    struct splicing sp;
    safecpy(sp.donor, sizeof(sp.donor), donor);
    safecpy(sp.accept, sizeof(sp.accept), accept);
    return sp;
}

/* load splicing from genome. */
static struct splicing loadSplicing(char *chrom, int start, int end, char strand,
                                    struct twoBitFile *genomeSeqs) {
    struct dnaSeq *dna1 = twoBitReadSeqFragLower(genomeSeqs, chrom, start, start+2);
    struct dnaSeq *dna2 = twoBitReadSeqFragLower(genomeSeqs, chrom, end-2, end);
    strUpper(dna1->dna);
    strUpper(dna2->dna);
    struct splicing sp;
    if (strand == '+') {
        sp = splicingMk(dna1->dna, dna2->dna);
    } else {
        reverseComplement(dna1->dna, strlen(dna1->dna));
        reverseComplement(dna2->dna, strlen(dna2->dna));
        sp = splicingMk(dna2->dna, dna1->dna);
    }
    dnaSeqFree(&dna1);
    dnaSeqFree(&dna2);
    return sp;
}

/* get donor/acceptor */
static void getIntronSplicing(struct twoBitFile* genomeSeqs,
                              struct intronInfo* intronInfo) {
    struct splicing sp = loadSplicing(intronInfo->chrom, 
                                      intronInfo->chromStart,
                                      intronInfo->chromEnd,
                                      intronInfo->transStrand[0],
                                      genomeSeqs);
    safecpy(intronInfo->transDonor, sizeof(intronInfo->transDonor), sp.donor);
    safecpy(intronInfo->transAcceptor, sizeof(intronInfo->transAcceptor), sp.accept);
}

/* get splice sites for all introns in transcripts */
static void getTranscriptSpliceSites(char* twoBitFile,
                                     struct intronMap* intronMap) {
    struct twoBitFile* genomeSeqs = twoBitOpen(twoBitFile);
    struct intronInfo* intronInfos = intronMapGetSorted(intronMap);
    for (struct intronInfo* intronInfo = intronInfos; intronInfo != NULL; intronInfo = intronInfo->next) {
        getIntronSplicing(genomeSeqs, intronInfo);
    }
    twoBitClose(&genomeSeqs);
}

/* collect into information  */
static void genePredIntrons(char* genePredFile, char* twoBitFile,
                            char* spliceTsv) {
    struct intronMap* intronMap = intronMapNew();
    intronMapLoadTranscripts(intronMap, genePredFile);
    getTranscriptSpliceSites(twoBitFile, intronMap);
    intronMapSaveTranscriptSpliceSites(intronMap, spliceTsv);
    intronMapFree(&intronMap);
}

/* entry */
int main(int argc, char** argv) {
    optionInit(&argc, argv, optionSpecs);
    if (argc != 4) {
        usage("wrong # args");
    }

    genePredIntrons(argv[1], argv[2], argv[3]);
    return 0;
}

