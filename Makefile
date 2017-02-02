ROOT = .
include ${ROOT}/config.mk

libs = arrayExpress.py config.py dataOps.py gencodeIntronEvid.py \
	pipelineOps.py rnaSeqData.py rnaSeqSupportToil.py \
	genbankProblemCasesDb.py
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
	flake8 ${libs:%=lib/gencode_icedb/%} ${progs:%=bin/%}

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -rf ${BINDIR}/*.dSYM ${OBJS} lib/gencode_icedb/*.pyc objs ${BINDIR}/*.pyc
