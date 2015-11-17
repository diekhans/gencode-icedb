/* convert fasta with UCSC names to GENCODE names */

#include "common.h"
#include "options.h"
#include "fa.h"
#include "hash.h"
#include "linefile.h"

static struct optionSpec optionSpecs[] = {
    {NULL, 0}
};

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "faUcscToGencode ucscFa nameMapTsv gencodeFa\n"
        "\n"
        "Convert fasta with UCSC names to GENCODE names.  Records that don't map are dropped.\n"
        "\n"
        "Options:\n";
    errAbort("%s:\n%s", msg, usageMsg);
}

/* load mapping of UCSC name into GENCODE name */
static struct hash* loadNameMap(char* nameMapTsv) {
    struct hash* nameMap = hashNew(0);
    struct lineFile* lf = lineFileOpen(nameMapTsv, TRUE);
    char* row[2];
    if (!lineFileNextRowTab(lf, row, 2)) {
        errAbort("empty file: %s", nameMapTsv);
    }
    while (lineFileNextRowTab(lf, row, 2)) {
        hashAdd(nameMap, row[0], cloneString(row[1]));
    }
    lineFileClose(&lf);
    return nameMap;
}

/* process on fasta record, outputting or skipping */
static void processRecord(struct hash* nameMap, struct dnaSeq* dnaSeq, FILE* outFh) {
    char* gencodeName = hashFindVal(nameMap, dnaSeq->name);
    if (gencodeName != NULL) {
        faWriteNext(outFh, gencodeName, dnaSeq->dna, dnaSeq->size);
    }
}

static void faUcscToGencode(char* ucscFa, char* nameMapTsv, char* gencodeFa) {
    struct hash* nameMap = loadNameMap(nameMapTsv);
    FILE* inFh = mustOpen(ucscFa, "r");
    FILE* outFh = mustOpen(gencodeFa, "w");

    struct dnaSeq* dnaSeq;
    while (faReadMixedNext(inFh, TRUE, NULL, TRUE, NULL, &dnaSeq)) {
        processRecord(nameMap, dnaSeq, outFh);
    }
    carefulClose(&inFh);
    carefulClose(&outFh);
}

/* entry */
int main(int argc, char** argv) {
    dnaUtilOpen();
    optionInit(&argc, argv, optionSpecs);
    if (argc != 4) {
        usage("wrong # args");
    }
    faUcscToGencode(argv[1], argv[2], argv[3]);
    return 0;
}
