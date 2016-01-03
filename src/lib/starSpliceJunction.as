object starSpliceJunction
"Splice junction output from STAR"
    (
    string chrom;             "chromsome"
    uint   chromStart;        "first base of the intron (1-based in file, converted in memory)"
    uint   chromEnd;          "last base of the intron"
    uint   strand;            "strand 0: undefined, 1: +, 2: -"
    uint   intronMotif;       "0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT"
    uint   annotated;         "0: unannotated, 1: annotated (only if splice junctions database is used)"
    uint   numUniqueMapReads; "number of uniquely mapping reads crossing the junction"
    uint   numMultiMapReads;  "number of multi-mapping reads crossing the junction"
    uint   maxOverhang;       "maximum spliced alignment overhang"
    )
