ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

test:: ucscTests ensemblTests problemCasesTests

##
# UCSC tsts
##
ucscTests: ucscRnaChrMTest ucscEstChrMTest

ucscRnaChrMTest: mkdirs
	rm -f output/$@.db
	${tslGetUcscRnaAligns} -chromSpec=chrM ${hsUcscGenomeDb} rna output/$@.db ${UCSC_RNA_ALN_TBL}
	sqlite3 --batch output/$@.db 'select * from ${UCSC_RNA_ALN_TBL};' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

ucscEstChrMTest: mkdirs
	rm -f output/$@.db
	${tslGetUcscRnaAligns} -chromSpec=chrM ${hsUcscGenomeDb} est output/$@.db ${UCSC_EST_ALN_TBL}
	sqlite3 --batch output/$@.db 'select * from ${UCSC_EST_ALN_TBL};' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

##
# Ensembl
##
ensemblTests: ensemblChrMTest ensemblIdListTest

ensemblChrMTest: mkdirs
	@rm -f output/$@.db
	${tslGetEnsemblRnaAligns} --chromSpec=chrM ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} output/$@.db ${ENSEMBL_RNA_ALN_TBL}
	sqlite3 --batch output/$@.db 'select * from ${ENSEMBL_RNA_ALN_TBL};' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

ensemblIdListTest: mkdirs
	@rm -f output/$@.db
	${tslGetEnsemblRnaAligns} --accverList=input/ensembl-cdna.acc ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} output/$@.db ${ENSEMBL_RNA_ALN_TBL}
	sqlite3 --batch output/$@.db 'select * from ${ENSEMBL_RNA_ALN_TBL};' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

##
# GBFF problem case parser
##
problemCasesTests: problemCasesParseTest problemCasesLoadTest

problemCasesParseTest: mkdirs
	${tslGbffGetProblemCases} --maxProcess=2 output/hs.$@.tab output/mm.$@.tab input/problem1.gbff input/problem2.gbff
	diff expected/hs.$@.tab output/hs.$@.tab
	diff expected/mm.$@.tab output/mm.$@.tab

problemCasesLoadTest: mkdirs
	${tslGenbankProblemCasesLoad} expected/hs.problemCasesParseTest.tab output/$@.db problemCases
	sqlite3 --batch output/$@.db 'select * from problemCases;' > output/$@.tab
	diff expected/hs.problemCasesParseTest.tab output/$@.tab

mkdirs:
	@mkdir -p output

clean::
	rm -rf output
