MACHTYPE = $(shell uname -m)
SYS = $(shell uname -s)

.SECONDARY:

# edit to set to UCSC browser kent/src
KENTDIR = ${HOME}/kent/src

KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc
KENTLIBDIR = ${KENTDIR}/lib/${MACHTYPE}
KENTLIBS = ${KENTLIBDIR}/jkhgap.a ${KENTLIBDIR}/jkweb.a
LIBS = ${KENTLIBS} -lssl -lcrypto -lz -lpthread

ifeq (${HTSDIR},)
    HTSDIR = /hive/data/outside/htslib/${MACHTYPE}
    ifneq ($(wildcard ${HTSDIR}),)
        ifeq (${USE_HTS},)
            USE_HTS=1
        endif
    endif
endif
ifeq (${USE_HTS},1)
    HG_DEFS+=-DUSE_HTS
    USE_SAMTABIX=1
    SAMTABIXDIR = ${HTSDIR}
    SAMTABIXLIB=${HTSDIR}/libhts.a
    LIBS += ${SAMTABIXLIB}
endif

MYSQLLIBS = $(shell mysql_config --libs)
LIBS += ${MYSQLLIBS}


CFLAGS =  -Wall -Werror -Wno-sign-compare -std=c99
CDEBUG = -g -O0

CFLAGS += ${KENTINC} ${CDEBUG}

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs
