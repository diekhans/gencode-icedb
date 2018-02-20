ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

# IMPORTANT: use bin/getTestCaseInfo to obtain data for new tests
testCases = input/v27-test-cases.tsv

allGenes = $(shell tawk 'NR>1 && !seen[$$4]{seen[$$4] = 1; print $$4}' ${testCases})

gencodeDb = output/gencode.db
evidDb = output/evid.db
createTestData = bin/createTestData

testInputDone = output/test-input.done

test:: classifyUnitTests supportCollectGenesTest supportCollectMkJobsTest

classifyUnitTests: ${testInputDone}
	${PYTHON} classifyUnitTests.py

supportCollectGenesTest: ${testInputDone} mkdirs
	${tslGencodeCollectSupport} ${gencodeDb} ${evidDb} --details=output/$@.details.tsv output/$@.tsl.tsv ENSG00000177663.13
	${diff} expected/$@.tsl.tsv output/$@.tsl.tsv
	${diff} expected/$@.details.tsv output/$@.details.tsv

supportCollectMkJobsTest: ${testInputDone} mkdirs
	rm -rf output/$@.work output/$@.db
	${tslGencodeCollectSupportMkJobs} --allDetails ${gencodeDb} ${evidDb} output/$@.work
	sed -Ee 's#^.*(/bin/tslGencodeCollectSupport)#??\1#' output/$@.work/batch.jobs > output/$@.batch.jobs.cleaned
	${diff} expected/$@.batch.jobs.cleaned output/$@.batch.jobs.cleaned
	${jobsToCmds} output/$@.work/batch.jobs | bash -e  
	ls -1 output/$@.work/results/ > output/$@.results.lst
	${diff} expected/$@.results.lst output/$@.results.lst
	${tslGencodeCollectSupportFinishJobs} --saveDetails output/$@.work output/$@.db
	$(call sqldumpdiff,gencode_transcript_support)
	$(call sqldumpdiff,gencode_transcript_support_details)


${testInputDone}: ${createTestData} ${testCases}
	@mkdir -p output
	${createTestData} ${testCases} ${hsGencodeVer} ${hsUcscGenomeDb} ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${evidDb} ${gencodeDb}
	touch $@


mkdirs:
	@mkdir -p output

clean::
	rm -rf output
