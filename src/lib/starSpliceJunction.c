/* starSpliceJunction.c was originally generated by the autoSql program, which also 
 * generated starSpliceJunction.h and starSpliceJunction.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "starSpliceJunction.h"
#include "rslAnalysisSet.h"


char *starSpliceJunctionCommaSepFieldNames = "chrom,chromStart,chromEnd,strand,intronMotif,annotated,numUniqueMapReads,numMultiMapReads,maxOverhang";

struct starSpliceJunction *starSpliceJunctionLoad(char **row)
/* Load a starSpliceJunction from row fetched with select * from starSpliceJunction
 * from database.  Dispose of this with starSpliceJunctionFree(). */
{
struct starSpliceJunction *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]) - 1;  /* CUSTOMIZED */
ret->chromEnd = sqlUnsigned(row[2]);
ret->strand = sqlUnsigned(row[3]);
ret->intronMotif = sqlUnsigned(row[4]);
ret->annotated = sqlUnsigned(row[5]);
ret->numUniqueMapReads = sqlUnsigned(row[6]);
ret->numMultiMapReads = sqlUnsigned(row[7]);
ret->maxOverhang = sqlUnsigned(row[8]);
return ret;
}

struct starSpliceJunction *starSpliceJunctionLoadAll(char *fileName) 
/* Load all starSpliceJunction from a whitespace-separated file.
 * Dispose of this with starSpliceJunctionFreeList(). */
{
struct starSpliceJunction *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[9];

while (lineFileRow(lf, row))
    {
    el = starSpliceJunctionLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct starSpliceJunction *starSpliceJunctionLoadAllByChar(char *fileName, char chopper) 
/* Load all starSpliceJunction from a chopper separated file.
 * Dispose of this with starSpliceJunctionFreeList(). */
{
struct starSpliceJunction *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[9];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = starSpliceJunctionLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct starSpliceJunction *starSpliceJunctionCommaIn(char **pS, struct starSpliceJunction *ret)
/* Create a starSpliceJunction out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new starSpliceJunction */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s) - 1;
ret->chromEnd = sqlUnsignedComma(&s);
ret->strand = sqlUnsignedComma(&s);
ret->intronMotif = sqlUnsignedComma(&s);
ret->annotated = sqlUnsignedComma(&s);
ret->numUniqueMapReads = sqlUnsignedComma(&s);
ret->numMultiMapReads = sqlUnsignedComma(&s);
ret->maxOverhang = sqlUnsignedComma(&s);
*pS = s;
return ret;
}

void starSpliceJunctionFree(struct starSpliceJunction **pEl)
/* Free a single dynamically allocated starSpliceJunction such as created
 * with starSpliceJunctionLoad(). */
{
struct starSpliceJunction *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
rslAnalysisLinkFreeList(&el->srcAnalyses);
freez(pEl);
}

void starSpliceJunctionFreeList(struct starSpliceJunction **pList)
/* Free a list of dynamically allocated starSpliceJunction's */
{
struct starSpliceJunction *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    starSpliceJunctionFree(&el);
    }
*pList = NULL;
}

void starSpliceJunctionOutput(struct starSpliceJunction *el, FILE *f, char sep, char lastSep) 
/* Print out starSpliceJunction.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->chromStart + 1);
fputc(sep,f);
fprintf(f, "%u", el->chromEnd);
fputc(sep,f);
fprintf(f, "%u", el->strand);
fputc(sep,f);
fprintf(f, "%u", el->intronMotif);
fputc(sep,f);
fprintf(f, "%u", el->annotated);
fputc(sep,f);
fprintf(f, "%u", el->numUniqueMapReads);
fputc(sep,f);
fprintf(f, "%u", el->numMultiMapReads);
fputc(sep,f);
fprintf(f, "%u", el->maxOverhang);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */
