ROOT = .
include ${ROOT}/config.mk

pylibs = $(wildcard ${PYLIBDIR}/*.py) $(wildcard ${PYLIBDIR}/rsl/*.py) $(wildcard ${PYLIBDIR}/tsl/*.py)

progs = encodeDccQuery estimateReadLength rnaSeqIntronEvidBed \
	rnaSeqIntronSupport \
	icedbProgSetup.py mkRnaSeqSupportBatch rnaSeqDataRegister runRnaSeqSupport \
	starGenerateGenome starSpliceJunctionMap \
	getEnsemblRnaAligns \
	gbffGetProblemCases genbankProblemCasesLoad
#gencodeIntronEvidStats gencodeIntronEvidPlot 

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
