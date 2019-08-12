# Before including, declare:
#  - PROGS - list of program names to compile, doesn't include BINDIR.
#  - progName_PROG_MODS = base name (less .o or .c) of modules to build for program
#  - LIB_MODS = base names files in current directory to compile to library


PYTHON = python3

MACHTYPE = $(shell uname -m)
SYS = $(shell uname -s)

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs

PYLIBDIR = ${ROOT}/lib/gencode_icedb

# programs
faUcscToGencode = ${BINDIR}/faUcscToGencode
genePredIntrons = ${BINDIR}/genePredIntrons
icedbLoadPslAligns = ${BINDIR}/icedbLoadPslAligns
rslEncodeDccQuery = ${BINDIR}/rslEncodeDccQuery
rslGencodeCollectNovel = ${BINDIR}/rslGencodeCollectNovel
rslGencodeCollectNovelFinishJobs = ${BINDIR}/rslGencodeCollectNovelFinishJobs
rslGencodeCollectNovelMkJobs = ${BINDIR}/rslGencodeCollectNovelMkJobs
rslGencodeCollectSupport = ${BINDIR}/rslGencodeCollectSupport
rslGencodeCollectSupportFinishJobs = ${BINDIR}/rslGencodeCollectSupportFinishJobs
rslGencodeCollectSupportMkJobs = ${BINDIR}/rslGencodeCollectSupportMkJobs
rslMappingMetadataDbLoad = ${BINDIR}/rslMappingMetadataDbLoad
rslMkStarSjOutSplits = ${BINDIR}/rslMkStarSjOutSplits
rslMkStarSjSupMergeJobs = ${BINDIR}/rslMkStarSjSupMergeJobs
sraRunInfoDbLoad = ${BINDIR}/sraRunInfoDbLoad
sraRunInfoFilter = ${BINDIR}/sraRunInfoFilter
rslStarSjOutSplit = ${BINDIR}/rslStarSjOutSplit
rslArrayExpressSourceLoad = ${BINDIR}/rslArrayExpressSourceLoad
tslGbffGetProblemCases = ${BINDIR}/tslGbffGetProblemCases
tslGenbankProblemCasesLoad = ${BINDIR}/tslGenbankProblemCasesLoad
tslCollectSupport = ${BINDIR}/tslCollectSupport
tslCollectSupportFinishJobs = ${BINDIR}/tslCollectSupportFinishJobs
tslCollectSupportMkJobs = ${BINDIR}/tslCollectSupportMkJobs
tslGetEnsemblRnaAligns = ${BINDIR}/tslGetEnsemblRnaAligns
tslGetUcscRnaAligns = ${BINDIR}/tslGetUcscRnaAligns
tslLoadAlignEvid = ${BINDIR}/tslLoadAlignEvid
ucscGencodeDbLoad = ${BINDIR}/ucscGencodeDbLoad



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
LIBS += ${KENTLIBS} -lssl -lcrypto -lz -lpthread -lstdc++

# autodetect UCSC installation of hal:
ifeq (${HALDIR},)
    HALDIR = /hive/groups/browser/hal/halRelease
    ifneq ($(wildcard ${HALDIR}),)
        ifeq (${USE_HAL},)
          USE_HAL=1
        endif
    endif
endif

ifeq (${USE_HAL},1)
    HALLIBS=${HALDIR}/lib/halMaf.a ${HALDIR}/lib/halChain.a ${HALDIR}/lib/halMaf.a ${HALDIR}/lib/halLiftover.a ${HALDIR}/lib/halLod.a ${HALDIR}/lib/halLib.a ${HALDIR}/lib/sonLib.a ${HALDIR}/lib/libhdf5_cpp.a ${HALDIR}/lib/libhdf5.a ${HALDIR}/lib/libhdf5_hl.a -lstdc++
    LIBS += ${HALLIBS}
endif

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


diff = diff -u
# force sort to be consistent
export LC_ALL=C

