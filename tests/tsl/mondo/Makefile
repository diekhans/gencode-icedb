##
# Semi-manual testing that builds a full database.  Used in development
# build not for regression, as it is slow.
##
ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

dbDir = db.d
hsEvidDb = ${dbDir}/hs.evid.db
hsAnnonDb = ${dbDir}/hs.annon.db
mmEvidDb = ${dbDir}/mm.evid.db
mmAnnonDb = ${dbDir}/mm.annon.db

tslDb = output/tsl.db
cmpTsv = output/cmp.tsv
detailsTsv = output/details.tsv

genesPerJob = 10
allDetails = --allDetails
# saveDetails = --saveDetails

tslCmp = bin/tslCmp
combineDetailsTsv = bin/combineDetailsTsv

test::
	@echo "run make setup or make finishup"

setup: ${hsAnnonDb} ${hsEvidDb}
	${tslGencodeCollectSupportMkJobs} --genesPerJob=${genesPerJob} ${allDetails} ${hsAnnonDb} ${hsEvidDb} output

finishup: ${cmpTsv} ${tslDb} ${detailsTsv}

tsldb: ${tslDb}
${tslDb}:
	@mkdir -p $(dir $@)
	${tslGencodeCollectSupportFinishJobs} ${saveDetails} output ${tslDb}.tmp
	mv -f ${tslDb}.tmp ${tslDb}


compare: ${cmpTsv}
${cmpTsv}: ${tslDb} ${tslCmp}
	@mkdir -p $(dir $@)
	${tslCmp} ${hsAnnonDb} ${tslDb} $@.tmp
	cut -f 4 $@.tmp |sort | uniq -c > $@.cnts
	mv -f  $@.tmp $@

details: ${detailsTsv}
${detailsTsv}:
	@mkdir -p $(dir $@)
	${combineDetailsTsv} output/results $@.tmp
	mv -f  $@.tmp $@

# database builds

dbbuild: ${hsAnnonDb} ${hsEvidDb} ${mmAnnonDb} ${mmEvidDb}

##
# Building of evidence db (~1.25 hours)
##
${hsEvidDb}:
	@mkdir -p $(dir $@)
	rm -f ${hsEvidDb} ${hsEvidDb}.tmp
	${tslGetEnsemblRnaAligns} ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${hsEvidDb}.tmp
	${tslGetUcscRnaAligns} ${hsUcscGenomeDb} rna ${hsEvidDb}.tmp
	${tslGetUcscRnaAligns} ${hsUcscGenomeDb} est ${hsEvidDb}.tmp
	${tslGenbankProblemCasesLoad} ${ROOT}/data/hs.genbank-problem-cases.tsv ${hsEvidDb}.tmp
	mv -f ${hsEvidDb}.tmp ${hsEvidDb}

${mmEvidDb}:
	@mkdir -p $(dir $@)
	rm -f ${mmEvidDb} ${mmEvidDb}.tmp
	${tslGetEnsemblRnaAligns} ${mmEnsemblCDnaDb} ${mmGrcRefAssemblyReport} ${mmEvidDb}.tmp
	${tslGetUcscRnaAligns} ${mmUcscGenomeDb} rna ${mmEvidDb}.tmp
	${tslGetUcscRnaAligns} ${mmUcscGenomeDb} est ${mmEvidDb}.tmp
	${tslGenbankProblemCasesLoad} ${ROOT}/data/mm.genbank-problem-cases.tsv ${mmEvidDb}.tmp
	mv -f ${mmEvidDb}.tmp ${mmEvidDb}

##
# building annotation db (~20 seconds)
##
${hsAnnonDb}:
	@mkdir -p $(dir $@)
	rm -f ${hsAnnonDb} ${hsAnnonDb}.tmp
	${ucscGencodeDbLoad} --hgdb=${hsUcscGenomeDb} --version=${hsGencodeVer} ${hsAnnonDb}.tmp
	mv -f ${hsAnnonDb}.tmp ${hsAnnonDb}

${mmAnnonDb}:
	@mkdir -p $(dir $@)
	rm -f ${mmAnnonDb} ${mmAnnonDb}.tmp
	${ucscGencodeDbLoad} --hgdb=${mmUcscGenomeDb} --version=${mmGencodeVer} ${mmAnnonDb}.tmp
	mv -f ${mmAnnonDb}.tmp ${mmAnnonDb}

clean::
	rm -rf output

realclean: clean
	rm -rf ${evidDir} ${annonDir}
