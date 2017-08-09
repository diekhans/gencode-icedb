ROOT = .
include ${ROOT}/config.mk

TSLCLTESTDIR = tests/tsl/classify

pylibs = $(wildcard ${PYLIBDIR}/*.py) $(wildcard ${PYLIBDIR}/rsl/*.py) $(wildcard ${PYLIBDIR}/tsl/*.py)
pytests = $(wildcard ${TSLCLTESTDIR}/*.py)
progs = icedbProgSetup.py \
	gencodeDbLoad \
	tslGbffGetProblemCases tslGenbankProblemCasesLoad tslGetEnsemblRnaAligns \
	rslEncodeDccQuery \
	rslSraRunInfoFilter rslSraRunInfoDbLoad rslMappingMetadataDbLoad \
	rslMkStarSjOutSplits \
	rslGencodeCollectSupport

# not done:  rslRnaSeqIntronEvidPlot rslRnaSeqIntronEvidStats

all::
	(cd src && ${MAKE})

test:: all
	(cd tests && ${MAKE} test)

lint:
	flake8 ${pylibs} ${progs:%=bin/%} ${pytests}

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -f ${pylibs:%=%c} ${BINDIR}/*.pyc
	rm -rf  ${OBJDIR}
