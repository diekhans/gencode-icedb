ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

# IMPORTANT: use bin/getTestCaseInfo to obtain data for new tests
testCases = input/v28-test-cases.tsv

allGenes = $(shell tawk 'NR>1 && !seen[$$4]{seen[$$4] = 1; print $$4}' ${testCases})

gencodeDb = output/db/gencode.db
evidDbDir = output/db/evidDb
createTestData = bin/createTestData

testDbDone = output/db/db.done

test:: classifyUnitTests supportCollectGenesTest supportCollectMkJobsTest

classifyUnitTests: ${testDbDone}
	${PYTHON} classifyUnitTests.py

supportCollectGenesTest: ${testDbDone} mkdirs
	${tslCollectSupport} ${gencodeDb} 02c995e3-372c-4cde-b216-5d3376c51988 ${evidDbDir}/RNA.psl.gz --details=output/$@.details.tsv output/$@.support.tsv ENSG00000177663.13
	${diff} expected/$@.support.tsv output/$@.support.tsv
	${diff} expected/$@.details.tsv output/$@.details.tsv

supportCollectMkJobsTest: ${testDbDone} mkdirs
	rm -rf output/$@.work output/$@.db
	${tslCollectSupportMkJobs} --genesPerJob=2 --allDetails ${gencodeDb} ${evidDbDir} output/$@.work
	sed -Ee 's#^.*(/bin/tslCollectSupport)#??\1#' output/$@.work/batch.jobs > output/$@.batch.jobs.cleaned
	${diff} expected/$@.batch.jobs.cleaned output/$@.batch.jobs.cleaned
	${jobsToMake} output/$@.work/batch.jobs output/$@.mk
	${MAKE} -f output/$@.mk
	(cd output/$@.work/results && find . -type f) | sort > output/$@.results.lst
	${diff} expected/$@.results.lst output/$@.results.lst
	${tslCollectSupportFinishJobs} --details=output/$@.support-eval-details.tsv output/$@.work output/$@.db
	$(call sqldumpdiff,gencode_support_eval)


${testDbDone}: ${createTestData} ${testCases}
	@mkdir -p $(dir $@)
	${createTestData} ${testCases} ${hsGencodeVer} ${hsUcscGenomeDb} ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${gencodeDb} ${evidDbDir}
	touch $@


mkdirs:
	@mkdir -p output

clean::
	rm -rf output
