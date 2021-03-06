#ifndef starOps_h
#define starOps_h
#include "common.h"
#include "starSpliceJunction.h"

/* convert STAR strand code to a character */
INLINE char starStrandCodeToChar(unsigned strandCode) {
    switch (strandCode) {
        case 0: return '?'; 
        case 1: return '+';
        case 2: return '-';;
        default:
            errAbort("unknown STAR strand code: %d", strandCode);
            return '\0';
    }
}

/* convert STAR intron motif code to a string.  non-canonical
 * is returned as ??/?? */
INLINE char* starMotifCodeToStr(unsigned intronMotifCode) {
    switch (intronMotifCode) {
        case 0: return "\?\?/\?\?"; 
        case 1: return "GT/AG";
        case 2: return "GT/AG"; // rev(CT/AC)
        case 3: return "GC/AG";
        case 4: return "GC/AG"; // rev(CT/GC)
        case 5: return "AT/AC";
        case 6: return "AT/AC"; // rev(GT/AT)
        default:
            errAbort("unknown STAR intron motif code: %d", intronMotifCode);
            return NULL;
    }
}

#endif
