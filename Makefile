ROOT = .
include ${ROOT}/config.mk

pylibs = $(wildcard ${PYLIBDIR}/*.py) $(wildcard ${PYLIBDIR}/rsl/*.py) $(wildcard ${PYLIBDIR}/tsl/*.py)

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
	flake8 ${pylibs} ${progs:%=bin/%}

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -f  ${PYLIBDIR}/*.pyc ${BINDIR}/*.pyc
	rm -rf  ${OBJDIR}
