ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

# IMPORTANT: use bin/getTestCaseInfo to obtain data for new tests
testCases = input/v28-test-cases.tsv

allGenes = $(shell tawk 'NR>1 && !seen[$$4]{seen[$$4] = 1; print $$4}' ${testCases})

gencodeDb = output/db/gencode.db
genbankDbDir = output/db/evidDb

rnaName = GenBank-RNA
rnaUuid = 02c995e3-372c-4cde-b216-5d3376c51988
rnaPsl =  ${genbankDbDir}/GenBank-RNA.psl.gz

estName = GenBank-EST
estUuid = fcfc5121-7e07-42f6-93a2-331a513eeb2c
estPsl = ${genbankDbDir}/GenBank-EST.psl.gz

ucscDRnaName = ucsc-drna
ucscDRnaUuid = f108c86b-8e0b-4dbd-96ae-6e75d8566eac
ucscDRnaPsl = ${genbankDbDir}/${ucscDRnaName}.psl.gz
ucscDRnaBam = output/${ucscDRnaName}.bam

createTestData = bin/createTestData

testDbDone = output/db/db.done

test:: classifyUnitTests supportCollectGenesTest supportCollectMkJobsTest supportCollectMkJobsPrimaryTest ucscDRnaTest ucscDRnaBamTest

classifyUnitTests: ${testDbDone}
	${PYTHON} classifyUnitTests.py

supportCollectGenesTest: ${testDbDone} mkdirs
	${tslCollectSupport} ${gencodeDb} 02c995e3-372c-4cde-b216-5d3376c51988 ${genbankDbDir}/GenBank-RNA.psl.gz --details=output/$@.details.tsv output/$@.support.tsv ENSG00000177663.13
	${diff} expected/$@.support.tsv output/$@.support.tsv
	${diff} expected/$@.details.tsv output/$@.details.tsv

supportCollectMkJobsTest: ${testDbDone} mkdirs
	rm -rf output/$@.tmp output/$@.db
	${tslCollectSupportMkJobs} --genesPerJob=2 ${gencodeDb} ${rnaName} ${rnaUuid} ${rnaPsl} output/$@.tmp
	${jobsToMake} output/$@.tmp/${rnaName}/batch.jobs output/$@.RNA.mk
	${MAKE} -f output/$@.RNA.mk
	${tslCollectSupportMkJobs} --genesPerJob=2 ${gencodeDb} ${estName} ${estUuid} ${estPsl} output/$@.tmp
	${jobsToMake} output/$@.tmp/${estName}/batch.jobs output/$@.EST.mk
	${MAKE} -f output/$@.EST.mk
	(cd output/$@.tmp && find . -type f) | sort > output/$@.results.lst
	${diff} expected/$@.results.lst output/$@.results.lst
	${tslCollectSupportFinishJobs} ${rnaName} output/$@.tmp output/$@.db
	${tslCollectSupportFinishJobs} ${estName} output/$@.tmp output/$@.db
	$(call sqldumpdiff,gencode_support_eval)

# this doesn't do that much, since we only have primary data in path
supportCollectMkJobsPrimaryTest: ${testDbDone} mkdirs
	rm -rf output/$@.tmp output/$@.tmp output/$@.db
	${tslCollectSupportMkJobs} --primaryOnly --genesPerJob=2 ${gencodeDb} ${rnaName} ${rnaUuid} ${rnaPsl} output/$@.tmp

ucscDRnaTest: ${testDbDone} mkdirs
	rm -rf output/$@.tmp output/$@.db
	${tslCollectSupportMkJobs} --genesPerJob=2 --allDetails ${gencodeDb} ${ucscDRnaName} ${ucscDRnaUuid} ${ucscDRnaPsl} output/$@.tmp
	${jobsToMake} output/$@.tmp/${ucscDRnaName}/batch.jobs output/$@.mk
	${MAKE} -f output/$@.mk
	${tslCollectSupportFinishJobs} --details=output/$@.support-eval-details.tsv ${ucscDRnaName} output/$@.tmp output/$@.db
	$(call sqldumpdiff,gencode_support_eval)

ucscDRnaBamTest: ${testDbDone} ${ucscDRnaBam} mkdirs
	rm -rf output/$@.tmp output/$@.db
	${tslCollectSupportMkJobs} --genesPerJob=2 --allDetails ${gencodeDb} ${ucscDRnaName} ${ucscDRnaUuid} ${ucscDRnaBam} output/$@.tmp
	${jobsToMake} output/$@.tmp/${ucscDRnaName}/batch.jobs output/$@.mk
	${MAKE} -f output/$@.mk
	${tslCollectSupportFinishJobs} --details=output/$@.support-eval-details.tsv ${ucscDRnaName} output/$@.tmp output/$@.db
	$(call sqldumpdiff,gencode_support_eval)


${testDbDone}: ${createTestData} ${testCases}
	@mkdir -p $(dir $@)
	${createTestData} ${testCases} ${hsGencodeVer} ${hsUcscGenomeDb} ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${gencodeDb} ${genbankDbDir}
	${tslLoadAlignEvid} input/ucsc-drna.psl ${genbankDbDir} ${ucscDRnaName}
	touch $@

output/%.bam: input/%.sam
	samtools sort -O BAM -o $@ $<
	samtools index -b $@

mkdirs:
	@mkdir -p output

clean::
	rm -rf output
