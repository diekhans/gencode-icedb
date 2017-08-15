ROOT = .
include ${ROOT}/config.mk

pylibs = $(wildcard ${PYLIBDIR}/*.py) \
	$(wildcard ${PYLIBDIR}/general/*.py) \
	$(wildcard ${PYLIBDIR}/rsl/*.py) \
	$(wildcard ${PYLIBDIR}/tsl/*.py)
pytests = $(wildcard tests/general/*.py)
progs = icedbProgSetup.py \
	gencodeDbLoad \
	tslGbffGetProblemCases tslGenbankProblemCasesLoad tslGetEnsemblRnaAligns \
	rslEncodeDccQuery \
	rslSraRunInfoFilter rslSraRunInfoDbLoad rslMappingMetadataDbLoad \
	rslMkStarSjOutSplits \
	rslGencodeCollectSupport rslGencodeCollectSupportMkJobs rslGencodeCollectSupportFinishJobs \
	rslGencodeCollectNovel

# not done:  rslRnaSeqIntronEvidPlot rslRnaSeqIntronEvidStats

all::
	(cd src && ${MAKE})

test:: all
	(cd tests && ${MAKE} test)

mondoTest::
	(cd tests && ${MAKE} mondoTest)

lint:
	flake8 ${pylibs} ${progs:%=bin/%} ${pytests}

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -f ${pylibs:%=%c} ${BINDIR}/*.pyc
	rm -rf  ${OBJDIR}
