MACH = $(shell uname -m)
SYS = $(shell uname -s)

# edit to set to UCSC browser kent/src
KENTDIR = ${HOME}/kent/src

KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc
KENTLIBDIR = ${KENTDIR}/lib/${MACH}
KENTLIBS = ${KENTLIBDIR}/jkhgap.a ${KENTLIBDIR}/jkweb.a
LIBS = ${KENTLIBS} -lssl -lcrypto -lz -lpthread

# autodetect UCSC installation of samtabix
ifeq (${SAMTABIXDIR},)
    SAMTABIXDIR = /hive/data/outside/samtabix/${MACH}
    ifeq ($(wildcard ${SAMTABIXDIR}),)
        SAMTABIXDIR =
    endif
endif
ifneq (${SAMTABIXDIR},)
    LIBS += ${SAMTABIXDIR}/libsamtabix.a
endif


CFLAGS =  -Wall -Werror -Wno-sign-compare 
CDEBUG = -g -O0

CFLAGS += ${KENTINC} ${CXXDEBUG}

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs
