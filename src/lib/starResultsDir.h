#ifndef starResultsDirs_h
#define starResultsDirs_h

/* An star results directory TSV */
struct starResults {
    struct starResults *next;
    char *run_acc;
    char *mapping_param_symid;
    char *mapping_symid;
    char *sjout;   // this is converted to absolute path
};

/* load an STAR results directory TSV files */
struct starResults *starResultsDirLoad(char *starResultsDirTsv);

/* free STAR results dir list */
void starResultsDirFree(struct starResults *starResultsDir);


#endif
