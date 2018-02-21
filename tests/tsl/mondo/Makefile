##
# Semi-manual testing that builds a full database.  Used in development
# build not for regression, as it is slow.
##
ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

dbDir = db.d
evidDb = ${dbDir}/hsEvid.db
annonDb = ${dbDir}/hsAnnon.db

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

setup: ${annonDb} ${evidDb}
	${tslGencodeCollectSupportMkJobs} --genesPerJob=${genesPerJob} ${allDetails} ${annonDb} ${evidDb} output

finishup: ${cmpTsv} ${tslDb} ${detailsTsv}

tsldb: ${tslDb}
${tslDb}:
	@mkdir -p $(dir $@)
	${tslGencodeCollectSupportFinishJobs} ${saveDetails} output ${tslDb}.tmp
	mv -f ${tslDb}.tmp ${tslDb}


compare: ${cmpTsv}
${cmpTsv}: ${tslDb} ${tslCmp}
	@mkdir -p $(dir $@)
	${tslCmp} ${annonDb} ${tslDb} $@.tmp
	cut -f 4 $@.tmp |sort | uniq -c > $@.cnts
	mv -f  $@.tmp $@

details: ${detailsTsv}
${detailsTsv}:
	@mkdir -p $(dir $@)
	${combineDetailsTsv} output/results $@.tmp
	mv -f  $@.tmp $@

##
# Building of evidence db (~1.25 hours)
##

${evidDb}:
	@mkdir -p $(dir $@)
	rm -f ${evidDb} ${evidDb}.tmp
	${tslGetEnsemblRnaAligns} ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${evidDb}.tmp
	${tslGetUcscRnaAligns} ${hsUcscGenomeDb} rna ${evidDb}.tmp
	${tslGetUcscRnaAligns} ${hsUcscGenomeDb} est ${evidDb}.tmp
	${tslGenbankProblemCasesLoad} ${ROOT}/data/hs.genbank-problem-cases.tsv ${evidDb}.tmp
	mv -f ${evidDb}.tmp ${evidDb}

##
# building annotation db (~20 seconds)
##
${annonDb}:
	@mkdir -p $(dir $@)
	rm -f ${annonDb} ${annonDb}.tmp
	${gencodeDbLoad} --hgdb=${hsUcscGenomeDb} --version=${hsGencodeVer} ${annonDb}.tmp
	mv -f ${annonDb}.tmp ${annonDb}

clean::
	rm -rf output

realclean: clean
	rm -rf ${evidDir} ${annonDir}
