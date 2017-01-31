# Before including, declare:
#  - PROGS - list of program names to compile, doesn't include BINDIR.
#  - progName_PROG_MODS = base name (less .o or .c) of modules to build for program
#  - LIB_MODS = base names files in current directory to compile to library


MACHTYPE = $(shell uname -m)
SYS = $(shell uname -s)

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs


.SECONDARY:

CFLAGS =  -Wall -Werror -Wno-sign-compare -std=c99
CDEBUG = -g -O0

ICEDBINC = -I${ROOT}/src/lib
ICEDBLIB = ${OBJDIR}/libicedb.a

CFLAGS += ${KENTINC} ${ICEDBINC} ${CDEBUG}
LIBS += ${ICEDBLIB}

# edit to set to UCSC browser kent/src
KENTDIR = ${HOME}/kent/src

KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc
KENTLIBDIR = ${KENTDIR}/lib/${MACHTYPE}
KENTLIBS = ${KENTLIBDIR}/jkhgap.a ${KENTLIBDIR}/jkweb.a ${KENTDIR}/htslib/libhts.a
LIBS += ${KENTLIBS} -lssl -lcrypto -lz -lpthread

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

ifneq (${MED_OPT},)
   CFLAGS += -I${MED_OPT}/include
   LIBS += -L${MED_OPT}/lib
endif

LIBS += -lsqlite3
