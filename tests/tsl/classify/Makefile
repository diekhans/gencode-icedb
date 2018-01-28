ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/testDefs.mk

hsGencodeVer = V27
hsEnsemblVer = 91_38
hgEnsemblCDnaDb = homo_sapiens_cdna_${hsEnsemblVer}
hsGrcRefAssembly = GRCh38
hsGrcRefAssemblyReport = ../genbankEvidence/input/GCF_000001405.36_GRCh38.p10_assembly_report.txt
hsUcscGenomeDb = hg38

hsEvidDb = output/hsEvid.db
hsGencodeDb = output/hsGencode.db

testInputDone = output/test-input.done

test:: ${testInputDone}


${testInputDone}: bin/createTestData input/v27-test-cases.tsv
	@mkdir -p output
	bin/createTestData input/v27-test-cases.tsv ${hsGencodeVer} ${hsUcscGenomeDb} ${hgEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${hsEvidDb} ${hsGencodeDb}
	touch $@


mkdirs:
	@mkdir -p output

clean::
	rm -rf output
