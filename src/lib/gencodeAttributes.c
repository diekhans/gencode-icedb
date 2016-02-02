#include "common.h"
#include "gencodeAttributes.h"
#include "hash.h"
#include "encode/wgEncodeGencodeAttrs.h"
#include "linefile.h"

/* load attributes into a hash by index */
struct hash* gencodeAttributesLoad(char* attributesTsv) {
    struct hash* attribMap = hashNew(0);
    struct lineFile* tsvLf = lineFileOpen(attributesTsv, TRUE);
    lineFileSkip(tsvLf, 1); // skip header
    char* row[WGENCODEGENCODEATTRS_NUM_COLS];
    while (lineFileNextRowTab(tsvLf, row, WGENCODEGENCODEATTRS_NUM_COLS)) {
        struct wgEncodeGencodeAttrs* attrs = wgEncodeGencodeAttrsLoad(row);
        hashAdd(attribMap, attrs->transcriptId, attrs);
    }
    lineFileClose(&tsvLf);
    return attribMap;
}

/* free gencode attributes */
void gencodeAttributesFree(struct hash* attribMap) {
    hashFreeWithVals(&attribMap, wgEncodeGencodeAttrsFree);
}

/* look up a transcript attribute */
struct wgEncodeGencodeAttrs* gencodeAttributesGet(struct hash* attribMap,
                                                  char* transcriptId) {
    // deal with PAR hack (ENSTR0000431238.7)
    if (transcriptId[4] == 'R') {
        char fixedName[32];
        safecpy(fixedName, sizeof(fixedName), transcriptId);
        fixedName[4] = '0';
        return hashMustFindVal(attribMap, fixedName);
    } else {
        return hashMustFindVal(attribMap, transcriptId);
    }
}
