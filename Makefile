ROOT = .
include ${ROOT}/config.mk

TSLCLTESTDIR = tests/tsl/classify

pylibs = $(wildcard ${PYLIBDIR}/*.py) $(wildcard ${PYLIBDIR}/rsl/*.py) $(wildcard ${PYLIBDIR}/tsl/*.py)
pytests = $(wildcard ${TSLCLTESTDIR}/*.py)
progs = encodeDccQuery estimateReadLength gbffGetProblemCases genbankProblemCasesLoad \
	getEnsemblRnaAligns icedbProgSetup.py mkRnaSeqSupportBatch rnaSeqDataRegister \
	rnaSeqIntronEvidBed \
	rnaSeqIntronSupport runRnaSeqSupport sjCollectEvidence starGenerateGenome \
	starSpliceJunctionMap

# not done:  rnaSeqIntronEvidPlot rnaSeqIntronEvidStats

all::
	(cd src && ${MAKE})

test::
	(cd tests && ${MAKE} test)

lint:
	flake8 ${pylibs} ${progs:%=bin/%} ${pytests}

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -f  ${PYLIBDIR}/*.pyc ${BINDIR}/*.pyc ${TSLCLTESTDIR}/*.pyc
	rm -rf  ${OBJDIR}
