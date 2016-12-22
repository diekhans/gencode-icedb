#include "common.h"
#include "options.h"
#include "psl.h"
#include "pslEvidCollect.h"
#include "pslEvidence.h"

static struct optionSpec optionSpecs[] = {
    {"ignoreMatch", OPTION_BOOLEAN},
    {NULL, 0}
};
static bool ignoreMatch = FALSE;

/* usage message and abort */
static void usage(char *msg) {
    static char* usageMsg = "cDnaPslEvidence cdnaPslFile twoBitFile cdnaAlignFile\n"
        "\n"
        "cdnaAlignFile are a tab-separated file describe in pslEvidCollect.c.\n"
        "\n"
        "Options:\n"
        "  -ignoreMatch - match/mismatch counts are incorrect and should be ignored\n";
    errAbort("%s:\n%s", msg, usageMsg);
}

/* extract the evidence */
static void cDnaPslEvidence(char *cdnaPslFile, char *twoBitFile, char *cdnaAlignFile) {
    struct pslEvidCollect* pslEvidCollect = pslEvidCollectNew(twoBitFile, cdnaAlignFile, ignoreMatch);
    struct psl *psls = pslEvidenceLoad(cdnaPslFile, NULL);
    for (struct psl *psl = psls; psl != NULL; psl = psl->next) {
        pslEvidCollectAnalyze(pslEvidCollect, psl, 1);
        pslEvidCollectWrite(pslEvidCollect);
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
    cDnaPslEvidence(argv[1], argv[2], argv[3]);
    return 0;
}
