ROOT = .
include ${ROOT}/config.mk

progs = icedbProgSetup.py \
	gencodeDbLoad \
	tslGbffGetProblemCases tslGenbankProblemCasesLoad tslGetEnsemblRnaAligns \
	rslEncodeDccQuery \
	rslSraRunInfoFilter rslSraRunInfoDbLoad rslMappingMetadataDbLoad \
	rslMkStarSjOutSplits \
	rslGencodeCollectSupport rslGencodeCollectSupportMkJobs rslGencodeCollectSupportFinishJobs \
	rslGencodeCollectNovel rslGencodeCollectNovelMkJobs rslGencodeCollectNovelFinishJobs

all::
	(cd src && ${MAKE})

test:: all
	(cd tests && ${MAKE} test)

mondoTest::
	(cd tests && ${MAKE} mondoTest)

lint:
	flake8-3 tests/general lib/gencode_icedb ${progs:%=bin/%}

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -rf  ${OBJDIR}
	find . -type f -name '*.pyc' -exec rm -f {} \;
	find . -type d -name __pycache__ -exec rmdir {} \;
