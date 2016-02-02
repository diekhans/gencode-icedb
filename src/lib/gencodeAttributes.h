#ifndef gencodeAttributes_h
#define gencodeAttributes_h
struct hash;
struct wgEncodeGencodeAttrs;

/* load attributes into a hash by index */
struct hash* gencodeAttributesLoad(char* attributesTsv);

/* free gencode attributes */
void gencodeAttributesFree(struct hash* attribMap);

/* look up a transcript attribute */
struct wgEncodeGencodeAttrs* gencodeAttributesGet(struct hash* attribMap,
                                                  char* transcriptId);


#endif
