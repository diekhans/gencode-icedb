ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

test:: ucscTests ensemblTests problemCasesTests

##
# UCSC tsts
##
ucscTests: ucscRnaChrMTest ucscEstChrMTest

ucscRnaChrMTest: mkdirs
	@rm -f output/$@.db
	${tslGetUcscRnaAligns} -chromSpec=chrM ${hsUcscGenomeDb} rna output/$@.psl
	diff expected/$@.psl output/$@.psl

ucscEstChrMTest: mkdirs
	@rm -f output/$@.db
	${tslGetUcscRnaAligns} -chromSpec=chrM ${hsUcscGenomeDb} est output/$@.psl
	diff expected/$@.psl output/$@.psl

##
# Ensembl
##
ensemblTests: ensemblChrMTest ensemblIdListTest

ensemblChrMTest: mkdirs
	@rm -f output/$@.db
	${tslGetEnsemblRnaAligns} --chromSpec=chrM ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} output/$@.psl
	diff expected/$@.psl output/$@.psl

ensemblIdListTest: mkdirs
	@rm -f output/$@.db
	${tslGetEnsemblRnaAligns} --accverList=input/ensembl-cdna.acc ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} output/$@.psl
	diff expected/$@.psl output/$@.psl

##
# GBFF problem case parser
##
problemCasesTests: problemCasesParseTest problemCasesLoadTest

problemCasesParseTest: mkdirs
	${tslGbffGetProblemCases} output/hs.$@.tsv output/mm.$@.tsv input/problem1.gbff input/problem2.gbff
	diff expected/hs.$@.tsv output/hs.$@.tsv
	diff expected/mm.$@.tsv output/mm.$@.tsv

problemCasesLoadTest: mkdirs
	@rm -f output/$@.db
	${tslGenbankProblemCasesLoad} expected/hs.problemCasesParseTest.tsv output/$@.db
	sqlite3 --batch --header output/$@.db 'select organism, etype, startAcc, endAcc, reason from genbank_problem_case;' > output/$@.tsv
	diff expected/hs.problemCasesParseTest.tsv output/$@.tsv

mkdirs:
	@mkdir -p output

clean::
	rm -rf output
