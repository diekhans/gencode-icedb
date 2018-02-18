##
# Semi-manual testing that builds a full database.  Used in development
# build not for regression, as it is slow.
##
ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

evidDir = evid.d
evidDb = ${evidDir}/hsEvid.db

annonDir = annon.d
annonDb = ${annonDir}/hsAnnon.db

tslDb = output/tsl.db

genesPerJob = 10


test::
	@echo "run make setup or make finishup"

setup: ${annonDb} ${evidDb}
	${tslGencodeCollectSupportMkJobs} --genesPerJob=${genePerJob} ${annonDb} ${evidDb} output

finishup:
	@mkdir -p $(dir ${tslDb})
	${tslGencodeCollectSupportFinishJobs} output ${tslDb}

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
