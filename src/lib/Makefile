ROOT = ../..

LIB_MODS =  intronMap starResultsDir starSpliceJunction rslAnalysisSet sqliteEz

include ${ROOT}/config.mk

include ${ROOT}/rules.mk

all:: ${LIB_OBJS} ${ICEDBLIB}(${LIB_OBJS})

clean::
	rm -f ${ICEDBLIB}


# Local Variables:
# mode: makefile-gmake
# End:
