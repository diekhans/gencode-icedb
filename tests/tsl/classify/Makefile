ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

hsEvidDb = output/hsEvid.db
hsGencodeDb = output/hsGencode.db

testInputDone = output/test-input.done

test:: ${testInputDone} classifyUnitTests

classifyUnitTests: ${testInputDone}
	${PYTHON} classifyUnitTests.py


supportCollectGenesTest: ${testInputDone} mkdirs
	${tslGencodeCollectSupport} ${genomeSeqsSpecs} output/$@/gencode.db output/$@/rnaseq.db output/$@.supp.tsv ENSG00000235169.7
	${diff} expected/$@.supp.tsv output/$@.supp.tsv


${testInputDone}: bin/createTestData input/v27-test-cases.tsv
	@mkdir -p output
	bin/createTestData input/v27-test-cases.tsv ${hsGencodeVer} ${hsUcscGenomeDb} ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${hsEvidDb} ${hsGencodeDb}
	touch $@


mkdirs:
	@mkdir -p output

clean::
	rm -rf output
